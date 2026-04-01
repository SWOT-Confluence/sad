"""
swot.py

Main driver script for the SAD algorithm with SWOT data.
Translated from swot.jl.
"""

import argparse
import json
import os
import sys
import traceback

import numpy as np
import netCDF4 as nc

import jax
jax.config.update("jax_enable_x64", True)

from preprocess import preprocess
from priors import priors_from_sos, priors_from_qwbm, River
from infer import infer
from utils import strip_mask

FILL = -999999999999


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def get_reach_files(indir, tmpdir, reachjson, index, sosbucket):
    """
    Get reach file names and optionally download SoS file from S3.
    """
    with open(os.path.join(indir, reachjson)) as f:
        reachlist = json.load(f)[index]

    if sosbucket:
        sosfile = os.path.join(tmpdir, reachlist["sos"])
        sys.path.insert(0, "/app/sos_read")
        from sos_read import download_sos
        download_sos(sosbucket, sosfile)
    else:
        sosfile = os.path.join(indir, "sos", reachlist["sos"])

    return (
        reachlist["reach_id"],
        os.path.join(indir, "swot",  reachlist["swot"]),
        sosfile,
        os.path.join(indir, "sword", reachlist["sword"]),
    )


# ---------------------------------------------------------------------------
# SWORD reader
# ---------------------------------------------------------------------------

def river_info(reach_id, swordfile):
    """
    Retrieve node IDs and chainage from SWORD NetCDF file.
    Returns (nids, x) sorted downstream -> upstream.
    """
    with nc.Dataset(swordfile) as fd:
        g       = fd.groups["nodes"]
        mask    = g.variables["reach_id"][:] == reach_id
        nid     = g.variables["node_id"][mask]
        x       = g.variables["dist_out"][mask]

        # Drop missing node IDs
        valid   = ~np.ma.getmaskarray(nid)
        nid     = strip_mask(nid[valid]).astype(int)
        x       = strip_mask(x[valid])

        # Reset chainage to start from 0 and sort downstream -> upstream
        x      -= x.min()
        order   = np.argsort(x)
        return nid[order], x[order]


# ---------------------------------------------------------------------------
# SWOT reader
# ---------------------------------------------------------------------------

def read_swot_obs(ncfile, nids):
    """
    Load SWOT node-level observations, ordered by nids (downstream -> upstream).

    Returns H, W, S as (nobs x nt) float64 arrays with NaN sentinels,
    and time_str as a list of strings.
    """
    with nc.Dataset(ncfile) as ds:
        nodes   = ds.groups["node"]
        reaches = ds.groups["reach"]

        # Node-level arrays: stored as (nt x nnodes) -> transpose to (nnodes x nt)
        S_raw   = strip_mask(nodes.variables["slope2"][:]).T
        H_raw   = strip_mask(nodes.variables["wse"][:]).T
        W_raw   = strip_mask(nodes.variables["width"][:]).T

        # Map node IDs to row indices
        node_ids = strip_mask(nodes.variables["node_id"][:]).astype(int)
        id_to_idx = {nid: k for k, nid in enumerate(node_ids)}
        rows = [id_to_idx[n] for n in nids if n in id_to_idx]

        H = H_raw[rows]
        W = W_raw[rows]
        S = S_raw[rows]

        # NaN-out W and S where H is NaN (consistency)
        nan_h       = np.isnan(H)
        W[nan_h]    = np.nan
        S[nan_h]    = np.nan

        # Reach-level time strings
        time_var    = reaches.variables["time_str"]
        time_raw    = np.array(time_var[:])
        if time_raw.ndim == 2:
            time_str = ["".join(row).strip() for row in time_raw]
        else:
            time_str = [str(t).strip() for t in time_raw]

        return H, W, S, time_str


# ---------------------------------------------------------------------------
# Output writer
# ---------------------------------------------------------------------------

def write_output(reachid, valid, outdir, A0, n, Qa, Qu, W, time_str):
    """Write SAD output to NetCDF."""
    outfile = os.path.join(outdir, f"{reachid}_sad.nc")
    nt      = W.shape[1] if W.ndim == 2 else len(time_str)
    nx      = W.shape[0] if W.ndim == 2 else 1

    with nc.Dataset(outfile, "w") as out:
        out.valid = int(valid)

        out.createDimension("nx", nx)
        out.createDimension("nt", nt)

        rid_v           = out.createVariable("reach_id", "i8", (), fill_value=FILL)
        rid_v[:]        = reachid

        A0_v            = out.createVariable("A0",  "f8", (), fill_value=FILL)
        A0_v[:]         = float(A0) if A0 is not None else FILL

        n_v             = out.createVariable("n",   "f8", (), fill_value=FILL)
        n_v[:]          = float(n)  if n  is not None else FILL

        Qa_v            = out.createVariable("Qa",  "f8", ("nt",), fill_value=FILL)
        Qa_arr          = np.asarray(Qa, dtype=np.float64)
        Qa_arr[np.isnan(Qa_arr)] = FILL
        Qa_v[:]         = Qa_arr

        Qu_v            = out.createVariable("Q_u", "f8", ("nt",), fill_value=FILL)
        Qu_arr          = np.asarray(Qu, dtype=np.float64)
        Qu_arr[np.isnan(Qu_arr)] = FILL
        Qu_v[:]         = Qu_arr

        ts_v            = out.createVariable("time_str", str, ("nt",))
        ts_v[:]         = np.array(time_str[:nt], dtype=object)


# ---------------------------------------------------------------------------
# Posterior uncertainty
# ---------------------------------------------------------------------------

def posterior_uncertainty(res, nt):
    """
    Compute per-timestep posterior Q uncertainty (std dev) from the
    ensemble. Returns NaN where the timestep was skipped.
    """
    Qu = np.full(nt, np.nan)
    for t, A in enumerate(res["A_post"]):
        if A is not None:
            Qu[t] = float(A[0].std())
    return Qu


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="SAD algorithm with SWOT data")
    p.add_argument("--index", "-i", type=int, default=0,
                   help="Index of reach to run on")
    p.add_argument("--reachfile", "-r", type=str, default="reaches.json",
                   help="Name of reaches JSON file")
    p.add_argument("--bucketkey", "-b", type=str, default="",
                   help="Bucket and key prefix to download SoS from S3")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    indir  = os.path.join("/mnt", "data", "input")
    outdir = os.path.join("/mnt", "data", "output")
    tmpdir = "/tmp"

    args = parse_args()

    # AWS Batch array job index handling
    if args.index == -256:
        try:
            index = int(os.environ["AWS_BATCH_JOB_ARRAY_INDEX"])
        except KeyError:
            index = 0
    else:
        index = args.index

    reachfile = args.reachfile
    bucketkey = args.bucketkey
    print(f"Index: {index}")
    print(f"Reach File: {reachfile}")
    print(f"Bucket Key: {bucketkey}")

    reachid, swotfile, sosfile, swordfile = get_reach_files(
        indir, tmpdir, reachfile, index, bucketkey,
    )
    print(f"Reach ID: {reachid}")
    print(f"SWOT:     {swotfile}")
    print(f"SOS:      {sosfile}")
    print(f"SWORD:    {swordfile}")

    nids, x = river_info(reachid, swordfile)
    H, W, S, time_str = read_swot_obs(swotfile, nids)
    nt = H.shape[1]

    # Default invalid outputs
    A0 = None
    n  = None
    Qa = np.full(nt, np.nan)
    Qu = np.full(nt, np.nan)

    if np.all(np.isnan(H)) or np.all(np.isnan(W)) or np.all(np.isnan(S)):
        print(f"{reachid}: INVALID")
        write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)
        return

    try:
        reach  = preprocess(x, H, W, S)
        hmin   = reach.hmin
        priors = priors_from_sos(sosfile, hmin, reachid)

        if priors is None:
            print(f"{reachid}: INVALID, missing mean discharge in SoS")
            write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)
            return

        res = infer(priors, reach)

        Qa = res["Q_post"]
        Qu = posterior_uncertainty(res, nt)

        # Posterior reach parameter means from the reach ensemble
        reach_ens = res["reach_ensemble"]
        n  = float(reach_ens[0].mean())   # Manning n
        r  = float(reach_ens[1].mean())   # Dingman r
        z0 = float(reach_ens[2].mean())   # bed elevation
        # A0 = bankfull cross-sectional area at downstream node
        # Using Dingman formula at bankfull depth
        hbf_ds = float(reach.hbf(reach.x[0]))
        Yb_ds  = max(hbf_ds - z0, 0.01)
        Wb_ds  = float(reach.wbf(reach.x[0]))
        Ym_ds  = (r + 1) / r * Yb_ds
        A0     = Wb_ds * (Ym_ds / Yb_ds) ** (1.0 / r) * Yb_ds

        print(f"{reachid}: VALID")
        write_output(reachid, 1, outdir, A0, n, Qa, Qu, W, time_str)

    except Exception:
        traceback.print_exc()
        print(f"{reachid}: INVALID")
        write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)


if __name__ == "__main__":
    main()
