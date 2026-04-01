"""
rejection.py

Prior refinement for the SAD inference pipeline.

Workflow
--------
1. select_representative_timesteps  — stratify valid timesteps by downstream
                                      WSE quantiles, pick best-covered per bin
2. rejection_sample                 — refine (n, r, z0, α) prior once using
                                      representative timesteps
3. q_ensemble                       — narrow Q prior per timestep using the
                                      accepted reach parameter ensemble and
                                      observed downstream WSE

All arrays are float64 NumPy.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from preprocess import SWOTReach
from priors import SWOTPriors
from gvf import gvf_solve, gvf_solve_ensemble

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Representative timestep selection
# ---------------------------------------------------------------------------

def _coverage(reach: SWOTReach, t: int) -> int:
    """
    Count non-NaN H observations at timestep t across all SWOT nodes.
    More observations → better spatial constraint on reach parameters.
    """
    return int(np.sum(np.isfinite(reach.obs.H[:, t])))


def select_representative_timesteps(
    reach:  SWOTReach,
    n_bins: int = 5,
) -> list[int]:
    """
    Select timesteps that cover the full range of flow conditions by
    stratifying the downstream WSE distribution into quantile bins and
    picking the timestep with the best spatial observation coverage in
    each bin.

    Quantile binning ensures equal representation of low, mid, and high
    flow conditions regardless of the shape of the WSE distribution —
    equal-width bins would underrepresent high flow where WSE varies
    slowly with Q.

    Parameters
    ----------
    reach : SWOTReach
    n_bins : int
        Number of quantile bins. Default 5. Automatically capped at the
        number of valid timesteps.

    Returns
    -------
    List of selected timestep indices (into reach.H columns).
    """
    valid_ts = np.where(reach.valid)[0]
    if valid_ts.size == 0:
        return []

    h_ds = reach.H[0, valid_ts]   # downstream WSE at each valid timestep

    # Degenerate case: all downstream WSE values identical
    if np.allclose(h_ds, h_ds[0]):
        best = valid_ts[np.argmax([_coverage(reach, t) for t in valid_ts])]
        logger.warning(
            "All downstream WSE values identical — returning single timestep"
        )
        return [int(best)]

    n_bins = min(n_bins, valid_ts.size)

    # Quantile edges: equal number of timesteps per bin
    quantiles = np.linspace(0.0, 1.0, n_bins + 1)
    edges     = np.quantile(h_ds, quantiles)

    selected = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        if i == n_bins - 1:
            in_bin = np.where((h_ds >= lo) & (h_ds <= hi))[0]
        else:
            in_bin = np.where((h_ds >= lo) & (h_ds < hi))[0]

        if in_bin.size == 0:
            continue

        # Pick timestep with most non-NaN spatial observations
        coverages = [_coverage(reach, valid_ts[j]) for j in in_bin]
        best      = in_bin[np.argmax(coverages)]
        selected.append(int(valid_ts[best]))

    return selected


# ---------------------------------------------------------------------------
# Rejection sampling  →  (n, r, z0, α) ensemble
# ---------------------------------------------------------------------------

def rejection_sample(
    priors:       SWOTPriors,
    reach:        SWOTReach,
    N:            int   = 500,
    eps_rel:      float = 0.5,
    max_attempts: int   = 20_000,
    n_bins:       int   = 5,
) -> np.ndarray:
    """
    Sample from the prior and reject parameter sets whose GVF-predicted WSE
    deviates from observations by more than eps_rel (mean relative RMSE
    across representative timesteps).

    Time-invariant reach parameters (n, r, z0, α) are drawn once per
    attempt and evaluated jointly across all representative timesteps.
    Q is drawn independently per timestep since it varies in time, but
    is not stored in the output — Q is handled per-timestep by q_ensemble.

    Parameters
    ----------
    priors : SWOTPriors
    reach : SWOTReach
    N : int
        Target number of accepted samples.
    eps_rel : float
        Relative RMSE threshold for acceptance. Loosely tuned — the goal
        is to prune implausible parameter space, not do full inference.
        Start at 0.5; tighten if acceptance rate > 20%, loosen if < 1%.
    max_attempts : int
        Hard cap on total draws to prevent infinite loops.
    n_bins : int
        Number of quantile bins for representative timestep selection.

    Returns
    -------
    (4 × N_accepted) float64 array with rows [n, r, z0, α].
    Returns empty (4 × 0) array if no samples are accepted.
    """
    rep_ts = select_representative_timesteps(reach, n_bins=n_bins)
    if not rep_ts:
        logger.warning("No valid timesteps for rejection sampling")
        return np.empty((4, 0))

    # Reference WSE scale: median downstream WSE across valid timesteps
    H_ref = float(np.median(reach.H[0, reach.valid]))

    # Pre-extract observations for each representative timestep
    obs_per_t = []
    for t in rep_ts:
        good = np.where(np.isfinite(reach.obs.H[:, t]))[0]
        obs_per_t.append({
            "x":    reach.obs.x[good],
            "H":    reach.obs.H[good, t].astype(np.float64),
            "H_bc": float(reach.H[0, t]),
        })

    accepted = np.empty((4, N))
    n_acc    = 0
    n_try    = 0

    while n_acc < N and n_try < max_attempts:
        n_try += 1

        # Draw time-invariant reach parameters
        n_draw = float(priors.np.rvs())
        r_draw = float(priors.rp.rvs())
        z0     = float(priors.zp.rvs())
        alpha  = float(priors.ap.rvs())

        rmse_all = []
        failed   = False

        for obs in obs_per_t:
            Q_t  = float(priors.Qp.rvs())
            pred = gvf_solve(Q_t, n_draw, r_draw, alpha, z0,
                             obs["H_bc"], reach, saveat=obs["x"])
            if pred is None:
                failed = True
                break
            rmse = float(np.sqrt(np.mean((pred - obs["H"]) ** 2)))
            rmse_all.append(rmse)

        if failed:
            continue
        if np.mean(rmse_all) / H_ref > eps_rel:
            continue

        accepted[:, n_acc] = [n_draw, r_draw, z0, alpha]
        n_acc += 1

    rate = 100.0 * n_acc / n_try if n_try > 0 else 0.0
    logger.info(
        "Rejection sampling: accepted %d / %d (%.1f%%) across %d timesteps",
        n_acc, n_try, rate, len(rep_ts),
    )
    if rate > 20.0:
        logger.warning("Acceptance rate > 20%% — consider tightening eps_rel")
    if rate < 1.0 and n_acc < N:
        logger.warning("Acceptance rate < 1%% — consider loosening eps_rel")

    return accepted[:, :n_acc]


# ---------------------------------------------------------------------------
# Manning lower bound on Q
# ---------------------------------------------------------------------------

def _manning_q_lower(
    reach:  SWOTReach,
    t:      int,
    priors: SWOTPriors,
) -> float:
    """
    Estimate a physically grounded lower bound on Q at timestep t using
    Manning's equation applied to the observed WSE profile.

    The observed depth and slope must be consistent with at least Q_lower
    discharge. Uses the most conservative (highest n) prior estimate to
    give the lowest plausible Q.

    The bound is asymmetric by design:
    - Lower bound: tight, from physics (observed depth requires minimum Q)
    - Upper bound: loose, kept at ppf(0.999) — WSE alone cannot rule out
      high Q since n or width could be underestimated.

    Parameters
    ----------
    reach : SWOTReach
    t : int
        Timestep index.
    priors : SWOTPriors
        Used to get the 95th percentile of n (most conservative estimate).

    Returns
    -------
    float : lower bound on Q [m3/s], clipped to prior support.
    """
    h = reach.obs.H[:, t]
    good = np.where(np.isfinite(h))[0]
    if good.size < 2:
        return float(priors.Qp.ppf(0.001))

    x_good = reach.obs.x[good]
    h_good = h[good].astype(np.float64)

    # Median slope — robust to local WSE noise
    S_obs = float(np.median(np.diff(h_good) / np.diff(x_good)))
    S_obs = max(S_obs, 1e-5)

    # Conservative depth: minimum observed depth across valid nodes
    z_good = reach.z(x_good) + reach.hmin
    d_obs  = np.maximum(h_good - z_good, 0.1)
    d_min  = float(np.min(d_obs))

    # Conservative width: minimum observed width at valid nodes
    w_obs = reach.obs.W[good, t]
    w_obs = w_obs[np.isfinite(w_obs)]
    w_min = float(np.min(w_obs)) if w_obs.size > 0 else float(np.nanmin(reach.wbf(x_good)))
    w_min = max(w_min, 10.0)

    # Highest plausible n gives the lowest Q estimate
    n_hi = float(priors.np.ppf(0.95))

    # Manning Q at minimum depth, width and maximum n
    Q_lower = (1.0 / n_hi) * w_min * d_min ** (5.0 / 3.0) * S_obs ** 0.5

    # 50% safety factor — Manning at a single cross-section underestimates
    # true Q due to channel irregularity and 3D effects
    Q_lower *= 0.5

    # Clip to lower half of prior support — Manning cannot be used to rule
    # out Q above the median, only to establish a minimum
    return float(np.clip(Q_lower, priors.Qp.ppf(0.001), priors.Qp.ppf(0.5)))


# ---------------------------------------------------------------------------
# Q ensemble  ->  per-timestep Q draws conditioned on upstream WSE
# ---------------------------------------------------------------------------

def q_ensemble(
    priors:          SWOTPriors,
    reach_ensemble:  np.ndarray,
    reach:           SWOTReach,
    t:               int,
    eps_abs:         float = 1.0,
    n_tries:         int   = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Find Q values consistent with the observed upstream WSE, conditioned
    on the accepted reach parameter ensemble.

    Design choices:

    1. Upstream node conditioning: the downstream WSE is the ODE boundary
       condition and carries no Q information. The upstream node WSE varies
       strongly with Q and provides a meaningful acceptance criterion.

    2. Manning lower bound: the observed WSE depth and slope constrain the
       minimum plausible Q. This removes the systematic low-Q bias from
       purely uniform sampling at low-flow timesteps where the prior mean
       would otherwise sit well above the true Q. The upper bound remains
       at ppf(0.999) since WSE alone cannot rule out high Q.

    3. Uniform sampling between [Q_lo, Q_hi] ensures the ensemble spans
       the full plausible range for the LETKF to update against.

    Uses a single vmapped GVF call over all (N x n_tries) candidates.

    Parameters
    ----------
    priors : SWOTPriors
    reach_ensemble : (4 x N) float64
        Accepted reach parameters [n, r, z0, alpha] from rejection_sample.
    reach : SWOTReach
    t : int
        Timestep index.
    eps_abs : float
        Absolute WSE tolerance [m] for acceptance at the upstream node.
        Default 1.0 m.
    n_tries : int
        Number of Q draws per ensemble member. Default 10.

    Returns
    -------
    Q_accepted : (N_accepted,) float64
        Accepted Q values.
    member_idx : (N_accepted,) int
        Index into reach_ensemble columns for each accepted Q.
    Returns (empty, empty) if reach.valid[t] is False or upstream WSE
    is NaN.
    """
    empty = (np.empty(0, dtype=np.float64), np.empty(0, dtype=np.intp))

    if not reach.valid[t]:
        return empty

    H_col       = reach.obs.H[:, t]
    valid_nodes = np.where(np.isfinite(H_col))[0]
    if valid_nodes.size < 2:
        return empty

    upstream_node = valid_nodes[-1]
    H_obs_up      = float(H_col[upstream_node])
    x_up          = np.array([reach.obs.x[upstream_node]])
    H_bc          = float(reach.H[0, t])
    N             = reach_ensemble.shape[1]

    # Manning lower bound — physically grounded minimum Q
    # Upper bound from prior — keeps the full high-Q tail available
    Q_lo = _manning_q_lower(reach, t, priors)
    Q_hi = float(priors.Qp.ppf(0.999))

    logger.debug(
        "Timestep %d: Q sampling bounds [%.0f, %.0f] m3/s",
        t, Q_lo, Q_hi,
    )

    n_total  = N * n_tries
    Q_all    = np.random.uniform(Q_lo, Q_hi, n_total).astype(np.float64)

    tile_idx = np.tile(np.arange(N), n_tries)
    n_ens    = reach_ensemble[0, tile_idx]
    r_ens    = reach_ensemble[1, tile_idx]
    z0_ens   = reach_ensemble[2, tile_idx]
    a_ens    = reach_ensemble[3, tile_idx]

    HA = gvf_solve_ensemble(
        Q_all, n_ens, r_ens, a_ens, z0_ens,
        H_bc, reach, saveat=x_up,
    )   # (n_total, 1)

    within_tol = np.isfinite(HA[:, 0]) & (np.abs(HA[:, 0] - H_obs_up) < eps_abs)

    accepted_Q: dict[int, float] = {}
    for i in np.where(within_tol)[0]:
        member = int(tile_idx[i])
        if member not in accepted_Q:
            accepted_Q[member] = float(Q_all[i])

    if len(accepted_Q) == 0:
        logger.warning("q_ensemble: no Q accepted at timestep %d", t)
        return empty

    members = np.array(list(accepted_Q.keys()),   dtype=np.intp)
    Q_out   = np.array(list(accepted_Q.values()), dtype=np.float64)
    return Q_out, members
