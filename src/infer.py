"""
inference.py

LETKF state estimation for the SAD inference pipeline.

Workflow
--------
1. rejection_sample (rejection.py)  — refine (n, r, z0, α) prior once
2. q_ensemble (rejection.py)        — narrow Q prior per timestep
3. letkf                            — full state update (Q, n, r, z0, α)
                                      per timestep using spatial WSE profile
4. infer                            — orchestrates the full pipeline

Observation error
-----------------
R is fixed per-run via sigma_obs. Optionally, per-node errors are inflated
by the inverse square root of temporal data completeness — nodes that are
frequently missing are less reliable when they do appear.

State vector rows: [Q, n, r, z0, α]  (5 × N_ens)
All arrays are float64 NumPy unless inside a JAX call.
"""

from __future__ import annotations

import logging
from typing import Optional

from tqdm import tqdm

import numpy as np
from scipy.linalg import sqrtm

from preprocess import SWOTReach
from priors import SWOTPriors
from gvf import gvf_solve_ensemble
from rejection import (
    rejection_sample,
    q_ensemble,
    select_representative_timesteps,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Observation completeness weights
# ---------------------------------------------------------------------------

def _obs_completeness(reach: SWOTReach) -> np.ndarray:
    """
    Compute temporal data completeness for each SWOT node.

    Completeness is the fraction of timesteps with a valid (non-NaN) H
    observation at that node. Used to inflate observation error for nodes
    that are frequently missing — they are less reliable when they do appear.

    Parameters
    ----------
    reach : SWOTReach

    Returns
    -------
    (nobs,) float64 array of completeness values in (0, 1].
    Clipped to a minimum of 1/nt to avoid division by zero.
    """
    nobs, nt   = reach.obs.H.shape
    valid_frac = np.array([
        np.sum(np.isfinite(reach.obs.H[j, :])) / nt
        for j in range(nobs)
    ])
    return np.clip(valid_frac, 1.0 / nt, 1.0)


# ---------------------------------------------------------------------------
# LETKF
# ---------------------------------------------------------------------------

def letkf(
    A:      np.ndarray,
    d:      np.ndarray,
    HA:     np.ndarray,
    R:      np.ndarray,
    xs:     Optional[list[list[int]]] = None,
    ys:     Optional[list[list[int]]] = None,
    rho:    float = 1.05,
) -> np.ndarray:
    """
    Local Ensemble Transform Kalman Filter.

    Translated from kalman.jl.

    Parameters
    ----------
    A : (ndim, nens) float64
        State ensemble matrix.
    d : (nobs,) float64
        Observation vector.
    HA : (nobs, nens) float64
        Model-predicted observation ensemble.
    R : (nobs,) or (nobs, nobs) float64
        Observation error covariance. If 1-D, treated as diagonal.
    xs : list of lists of int, optional
        State indices for each local patch. Defaults to global update.
    ys : list of lists of int, optional
        Observation indices for each local patch.
    rho : float
        Covariance inflation factor. Default 1.05.

    Returns
    -------
    (ndim, nens) float64 analysis ensemble.
    """
    ndim, nens = A.shape
    nobs       = d.shape[0]
    Aa         = np.zeros_like(A)

    # Build R matrix
    if R.ndim == 1:
        R_mat = np.diag(R)
    else:
        R_mat = R

    # Ensemble perturbations
    Y = HA - HA.mean(axis=1, keepdims=True)   # (nobs, nens)
    X = A  - A.mean(axis=1,  keepdims=True)   # (ndim, nens)

    # Default: single global patch
    if xs is None or ys is None:
        xs = [list(range(ndim))]
        ys = [list(range(nobs))]

    for lx, ly in zip(xs, ys):
        lx = list(lx)
        ly = list(ly)

        Xl = X[np.ix_(lx, range(nens))]               # (|lx|, nens)
        Yl = Y[np.ix_(ly, range(nens))]               # (|ly|, nens)
        Rl = R_mat[np.ix_(ly, ly)]                    # (|ly|, |ly|)

        C  = Yl.T @ np.linalg.pinv(Rl)                                 # (nens, |ly|)
        P  = np.linalg.pinv((nens - 1) / rho * np.eye(nens) + C @ Yl) # (nens, nens)
        W  = np.real(sqrtm((nens - 1) * P))                            # (nens, nens)
        w  = P @ C @ (d[ly] - HA[ly, :].mean(axis=1))                  # (nens,)
        W  = W + w[:, np.newaxis]

        Aa[np.ix_(lx, range(nens))] = (
            Xl @ W + A[lx, :].mean(axis=1, keepdims=True)
        )

    return Aa


# ---------------------------------------------------------------------------
# Full inference loop
# ---------------------------------------------------------------------------

def infer(
    priors:              SWOTPriors,
    reach:               SWOTReach,
    sigma_obs:           float = 0.1,
    N:                   int   = 500,
    eps_rel:             float = 0.5,
    eps_abs_q:           float = 1.0,
    n_bins:              int   = 5,
    rho:                 float = 1.05,
    completeness_weight: bool  = True,
) -> dict:
    """
    Full SAD inference pipeline for a single SWOT reach.

    Parameters
    ----------
    priors : SWOTPriors
    reach : SWOTReach
    sigma_obs : float
        Observation error standard deviation [m]. Default 0.1 m.
        For Pepsi data (no measurement error) this reflects GVF structural
        error only. For real SWOT data combine measurement and model error:
        sigma_obs = sqrt(sigma_meas^2 + sigma_model^2).
    N : int
        Target rejection sample size (reach parameter ensemble).
    eps_rel : float
        Relative RMSE threshold for rejection sampling.
    eps_abs_q : float
        Absolute upstream WSE tolerance for q_ensemble [m].
    n_bins : int
        Quantile bins for representative timestep selection.
    rho : float
        LETKF covariance inflation factor.
    completeness_weight : bool
        If True, inflate per-node R by the inverse of temporal data
        completeness. Nodes with frequent missing data receive higher
        observation error. Default True.

    Returns
    -------
    dict with keys:
        "reach_ensemble" : (4 × N_acc) array  [n, r, z0, α]
        "Q_post"         : (nt,) array of posterior mean Q per timestep
                           (NaN for invalid or skipped timesteps)
        "A_post"         : list of (5 × N) posterior ensemble per timestep
                           (None for invalid or skipped timesteps)
        "rep_ts"         : list of representative timestep indices
        "completeness"   : (nobs,) array of per-node data completeness
    """
    # --- Stage 1: refine reach parameter prior ---
    rep_ts         = select_representative_timesteps(reach, n_bins=n_bins)
    reach_ensemble = rejection_sample(
        priors, reach, N=N, eps_rel=eps_rel, n_bins=n_bins,
    )
    if reach_ensemble.shape[1] == 0:
        raise RuntimeError("Rejection sampling produced no accepted samples.")

    N_acc = reach_ensemble.shape[1]
    logger.info("Reach ensemble size: %d", N_acc)

    # Pre-compute per-node completeness weights once
    completeness = _obs_completeness(reach)
    logger.info(
        "Node completeness: min=%.2f  median=%.2f  max=%.2f",
        completeness.min(), np.median(completeness), completeness.max(),
    )

    nt     = reach.nt
    Q_post = np.full(nt, np.nan)
    A_post = [None] * nt

    valid_ts = np.where(reach.valid)[0]
    for t in tqdm(valid_ts, desc="Timesteps", unit="t"):

        # --- Stage 2: Q ensemble conditioned on upstream WSE ---
        Q_ens, member_idx = q_ensemble(
            priors, reach_ensemble, reach, t, eps_abs=eps_abs_q,
        )
        if Q_ens.size == 0:
            logger.warning("Skipping timestep %d: empty Q ensemble", t)
            continue

        # Use original member indices — preserves (Q, n, r, z0, α) pairing
        n_ens  = reach_ensemble[0, member_idx]
        r_ens  = reach_ensemble[1, member_idx]
        z0_ens = reach_ensemble[2, member_idx]
        a_ens  = reach_ensemble[3, member_idx]

        # Full state ensemble: rows [log(Q), n, r, z0, α]
        # log-space enforces Q > 0 after the LETKF update and improves
        # Gaussianity of the Q marginal distribution
        A    = np.vstack([np.log(Q_ens), n_ens, r_ens, z0_ens, a_ens])  # (5, n_q)
        H_bc = float(reach.H[0, t])

        # --- Stage 3: forward model ---
        good = np.where(np.isfinite(reach.obs.H[:, t]))[0]
        if good.size < 2:
            logger.warning(
                "Skipping timestep %d: fewer than 2 observations", t
            )
            continue

        d     = reach.obs.H[good, t].astype(np.float64)
        x_obs = reach.obs.x[good]

        # A[0] is log(Q) — exponentiate for the GVF forward model
        HA_full = gvf_solve_ensemble(
            np.exp(A[0]), n_ens, r_ens, a_ens, z0_ens,
            H_bc, reach, saveat=x_obs,
        )   # (n_q, nobs)

        # Drop ensemble members where GVF failed
        good_ens = np.where(np.all(np.isfinite(HA_full), axis=1))[0]
        if good_ens.size < 2:
            logger.warning(
                "Skipping timestep %d: fewer than 2 valid ensemble members", t
            )
            continue

        A_filt  = A[:, good_ens]
        HA_filt = HA_full[good_ens].T   # (nobs, n_valid)

        n_q     = A_filt.shape[1]
        n_valid = good_ens.size
        if n_valid < n_q * 0.5:
            logger.warning(
                "Timestep %d: only %d / %d ensemble members valid after GVF",
                t, n_valid, n_q,
            )

        # --- Stage 4: build R ---
        if completeness_weight:
            sigma_node = sigma_obs / np.sqrt(completeness[good])
        else:
            sigma_node = np.full(d.size, sigma_obs)
        R_diag = sigma_node ** 2   # variance

        # --- Stage 5: LETKF update ---
        A_analysis = letkf(A_filt, d, HA_filt, R_diag, rho=rho)

        # Convert log(Q) back to natural space
        log_Q_post = A_analysis[0]
        Q_post_ens = np.exp(log_Q_post)

        if not np.all(np.isfinite(Q_post_ens)):
            logger.warning(
                "Timestep %d: non-finite Q after exponentiation — skipping", t
            )
            continue

        Q_post[t] = float(Q_post_ens.mean())

        # Store posterior in natural space [Q, n, r, z0, α]
        A_post[t] = np.vstack([
            Q_post_ens,
            A_analysis[1],
            A_analysis[2],
            A_analysis[3],
            A_analysis[4],
        ])

        logger.info(
            "Timestep %d: Q = %.1f m³/s (spread = %.1f)",
            t, Q_post[t], float(Q_post_ens.std()),
        )

    return {
        "reach_ensemble": reach_ensemble,
        "Q_post":         Q_post,
        "A_post":         A_post,
        "rep_ts":         rep_ts,
        "completeness":   completeness,
    }
