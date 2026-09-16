"""Run the SAD algorithm (SADnm) on a SWOT reach and write `<reach_id>_sad.nc`.

Replaces the Julia module (`swot.jl` + `Sad.jl`) while preserving its output
specification exactly: same variables, dimensions, fill value and `valid` attribute.

SADnm does not identify Manning's `n` -- roughness and slope are absorbed jointly into
the Stage-1 level constant. `n` and `A0` are therefore back-derived via
`sadnm.effective_cross_section` and are *effective* quantities that inherit the monthly
prior's level; see the SADnm README. Both are written as FILL when they cannot be formed.
"""

import argparse
import datetime
import json
import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

import sadnm
from sadnm.uniform_flow import effective_cross_section

FILL = -999999999999
TIME_FILL = "no_data"
EPOCH = datetime.datetime(2000, 1, 1)
INDEX_SENTINEL = -256
MIN_SLOPE_OBS = 3       # overpasses with a usable slope needed to form n_eff

INPUT_DIR = Path("/mnt/data/input")
OUTPUT_DIR = Path("/mnt/data/output")
TMP_DIR = Path("/tmp")


def parse_commandline():
    p = argparse.ArgumentParser(description="Run SADnm on one SWOT reach.")
    p.add_argument("-i", "--index", type=int, default=0,
                   help="Index of reach to run on")
    p.add_argument("-r", "--reachfile", type=str, default="reaches.json",
                   help="Name of reaches JSON file")
    p.add_argument("-b", "--bucketkey", type=str, default="",
                   help="Bucket and key prefix to download SoS from")
    return p.parse_args()


def get_reach_files(indir, tmpdir, reachjson, index, sosbucket):
    """Resolve the input paths for one reach, downloading the SoS if requested."""
    with open(indir / reachjson) as jf:
        reach = json.load(jf)[index]

    if sosbucket:
        sosfile = tmpdir / reach["sos"]
        from sos_read.sos_read import download_sos
        download_sos(sosbucket, sosfile)
    else:
        sosfile = indir / "sos" / reach["sos"]

    return (int(reach["reach_id"]), indir / "swot" / reach["swot"], sosfile,
            indir / "sword" / reach["sword"])


def river_info(reach_id, swordfile):
    """Node IDs for a reach, ordered downstream to upstream by `dist_out`."""
    with Dataset(swordfile) as ds:
        g = ds.groups["nodes"]
        i = np.where(g["reach_id"][:] == reach_id)[0]
        nid = g["node_id"][i]
        x = g["dist_out"][i]
        keep = ~(np.ma.getmaskarray(nid) | np.ma.getmaskarray(x))
        nid, x = np.ma.getdata(nid)[keep], np.ma.getdata(x)[keep]
        order = np.argsort(x - x.min())
        return nid[order].astype(np.int64)


def _as_time_by_node(arr, n_node):
    """Orient a node variable to (n_time, n_node).

    The node axis is identified by length rather than assumed, so this is correct
    whichever way round the file stores it.
    """
    arr = np.ma.filled(arr.astype(float), np.nan)
    if arr.shape[0] == n_node and arr.shape[1] != n_node:
        return arr.T
    return arr


def read_swot_obs(ncfile, nids):
    """Load node WSE/width and reach slope/time, node-ordered to match `nids`."""
    with Dataset(ncfile) as ds:
        nodes, reaches = ds.groups["node"], ds.groups["reach"]

        file_nid = np.ma.getdata(nodes["node_id"][:]).astype(np.int64)
        dmap = {int(n): k for k, n in enumerate(file_nid)}
        sel = [dmap[int(n)] for n in nids if int(n) in dmap]

        wse = _as_time_by_node(nodes["wse"][:], len(file_nid))[:, sel]
        width = _as_time_by_node(nodes["width"][:], len(file_nid))[:, sel]
        slope = np.ma.filled(reaches["slope2"][:].astype(float), np.nan)

        tvar = reaches["time"]
        traw = np.ma.filled(tvar[:].astype(float), np.nan)
        # Some reaches carry overpasses with no timestamp. They cannot be assigned to a
        # month, so they are excluded from the inversion via `time_ok` but still occupy a
        # slot on the output time axis.
        time_ok = np.isfinite(traw)
        time_s = np.zeros(len(traw))
        time_str = [TIME_FILL] * len(traw)
        if hasattr(tvar, "units"):
            dates = num2date(traw[time_ok], tvar.units,
                             only_use_cftime_datetimes=False, only_use_python_datetimes=True)
            for k, d in zip(np.flatnonzero(time_ok), dates):
                time_s[k] = (d - EPOCH).total_seconds()
                time_str[k] = d.isoformat()
        else:
            time_s[time_ok] = traw[time_ok]
            for k in np.flatnonzero(time_ok):
                time_str[k] = str(traw[k])

    # SWOT marks missing data with both masks and NaNs; treat them alike.
    bad = ~np.isfinite(wse)
    wse[bad] = np.nan
    width[bad] = np.nan
    node_id = np.asarray(nids)[[i for i, n in enumerate(nids) if int(n) in dmap]]
    return wse, width, slope, time_s, time_str, time_ok, node_id


def read_prior(sosfile, reach_id):
    """12-month discharge climatology for a reach, or None if unavailable."""
    with Dataset(sosfile) as ds:
        i = np.where(ds.groups["reaches"]["reach_id"][:] == reach_id)[0]
        if len(i) == 0:
            return None
        mq = np.ma.filled(ds.groups["model"]["monthly_q"][i[0], :].astype(float), np.nan)
    return mq if np.any(np.isfinite(mq) & (mq > 0)) else None


def reach_slope(slope):
    """Median of the usable per-overpass slopes.

    SWOT reports genuinely negative reach slopes on some overpasses, so non-positive
    values are dropped rather than clamped -- `n_eff` scales as sqrt(S), and a floored
    slope would yield a plausible-looking but badly wrong roughness.
    """
    s = slope[np.isfinite(slope) & (slope > 0)]
    return float(np.median(s)) if len(s) >= MIN_SLOPE_OBS else np.nan


def write_output(reachid, valid, outdir, A0, n, Qa, Qu, nx, time_str):
    """Write the SAD output NetCDF (same specification as the Julia module)."""
    outfile = Path(outdir) / f"{reachid}_sad.nc"
    with Dataset(outfile, "w") as out:
        out.valid = valid
        out.createDimension("nx", nx)
        out.createDimension("nt", len(time_str))

        v = out.createVariable("reach_id", "i8", (), fill_value=FILL)
        v[...] = reachid
        v = out.createVariable("A0", "f8", (), fill_value=FILL)
        v[...] = FILL if not np.isfinite(A0) else A0
        v = out.createVariable("n", "f8", (), fill_value=FILL)
        v[...] = FILL if not np.isfinite(n) else n
        v = out.createVariable("Qa", "f8", ("nt",), fill_value=FILL)
        v[:] = np.where(np.isfinite(Qa), Qa, FILL)
        v = out.createVariable("Q_u", "f8", ("nt",), fill_value=FILL)
        v[:] = np.where(np.isfinite(Qu), Qu, FILL)
        v = out.createVariable("time_str", str, ("nt",), fill_value=TIME_FILL)
        v[:] = np.array([t if t else TIME_FILL for t in time_str], dtype=object)


def main():
    args = parse_commandline()
    if args.index == INDEX_SENTINEL:
        index = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX", 0))
    else:
        index = args.index

    print(f"Index: {index}")
    print(f"Reach File: {args.reachfile}")
    print(f"Bucket Key: {args.bucketkey}")

    reachid, swotfile, sosfile, swordfile = get_reach_files(
        INPUT_DIR, TMP_DIR, args.reachfile, index, args.bucketkey)
    print(f"Reach ID: {reachid}")
    print(f"SWOT: {swotfile}")
    print(f"SOS: {sosfile}")
    print(f"SWORD: {swordfile}")

    nids = river_info(reachid, swordfile)
    wse, width, slope, time_s, time_str, time_ok, node_id = read_swot_obs(swotfile, nids)
    nt, nx = wse.shape

    empty = np.full(nt, np.nan)
    if not np.any(np.isfinite(wse)) or not np.any(np.isfinite(width)):
        print(f"{reachid}: INVALID, no usable SWOT observations")
        write_output(reachid, 0, OUTPUT_DIR, np.nan, np.nan, empty, empty, nx, time_str)
        return

    monthly_q = read_prior(sosfile, reachid)
    if monthly_q is None:
        print(f"{reachid}: INVALID, missing mean discharge")
        write_output(reachid, 0, OUTPUT_DIR, np.nan, np.nan, empty, empty, nx, time_str)
        return

    node_mask = np.isfinite(wse) & np.isfinite(width)
    overpass_mask = node_mask.any(axis=1) & time_ok
    params, _, inv_cfg = sadnm.load_config()

    res = sadnm.run_reach(
        np.nan_to_num(wse), np.nan_to_num(width), node_mask, overpass_mask, node_id,
        time_s, monthly_q,
        # SADnm de-normalises as x*std + mean; the wrapper supplies raw metres.
        dict(wse_mean=0.0, wse_std=1.0, width_mean=0.0, width_std=1.0),
        params, cfg=inv_cfg,
    )
    if not res.ok:
        print(f"{reachid}: INVALID, {res.reason}")
        write_output(reachid, 0, OUTPUT_DIR, np.nan, np.nan, empty, empty, nx, time_str)
        return

    Qa, Qu = empty.copy(), empty.copy()
    Qa[res.overpass_idx] = res.q
    # lognormal SD implied by the Stage-2 predictive log-sigma
    Qu[res.overpass_idx] = res.q * np.sqrt(np.expm1(res.log_sigma ** 2))

    ecs = effective_cross_section(res, reach_slope(slope))
    print(f"{reachid}: VALID (r={res.r_shape:.3f}, d0={res.d0:.2f}, n_eff={ecs.n_eff:.4f})")
    write_output(reachid, 1, OUTPUT_DIR, ecs.A0, ecs.n_eff, Qa, Qu, nx, time_str)


if __name__ == "__main__":
    main()
