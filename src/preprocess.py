"""
preprocess.py

Preprocessing for SWOT river reach observations.

Translated from preprocess.jl. Uses NaN sentinels throughout rather than
Union{Missing, Float64} — NaN means the satellite did not observe that
node at that timestep.

All arrays follow the convention: axis 0 = nodes (downstream → upstream),
axis 1 = timesteps.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator

from utils import strip_mask


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SWOTObs:
    """
    Raw SWOT observations at node locations.

    Attributes
    ----------
    x : (nobs,) float64
        Chainage of SWOT nodes, downstream → upstream [m].
    H : (nobs, nt) float64
        Water surface elevation observations [m]. NaN where not observed.
    W : (nobs, nt) float64
        Width observations [m]. NaN where not observed.
    """
    x: np.ndarray
    H: np.ndarray
    W: np.ndarray


@dataclass
class SWOTReach:
    """
    Preprocessed SWOT reach ready for GVF solving.

    Attributes
    ----------
    obs : SWOTObs
        Raw observations (NaN-sentineled).
    x : (nx,) float64
        Uniform computational chainage, downstream → upstream [m].
    H : (nx, nt) float64
        WSE interpolated to computational chainage. Zero for invalid timesteps.
    W : (nx, nt) float64
        Width interpolated to computational chainage. Zero for invalid timesteps.
    valid : (nt,) bool
        True where timestep has >= 2 valid H observations.
    S0 : PchipInterpolator
        Bed slope interpolant S0(x) [m/m].
    wbf : PchipInterpolator
        Bankfull width interpolant wbf(x) [m].
    hbf : PchipInterpolator
        Bankfull WSE interpolant hbf(x) [m].
    z : PchipInterpolator
        Cumulative bed elevation interpolant z(x) [m], z(0) = 0.
    hmin : float
        Minimum observed WSE at the downstream node [m].
    nx : int
        Number of computational nodes.
    nobs : int
        Number of SWOT nodes with at least one valid observation.
    nt : int
        Number of timesteps.
    """
    obs:   SWOTObs
    x:     np.ndarray
    H:     np.ndarray
    W:     np.ndarray
    valid: np.ndarray
    S0:    PchipInterpolator
    wbf:   PchipInterpolator
    hbf:   PchipInterpolator
    z:     PchipInterpolator
    hmin:  float
    nx:    int
    nobs:  int
    nt:    int


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _pchip(x: np.ndarray, y: np.ndarray) -> PchipInterpolator:
    """
    Build a PCHIP interpolant with linear extrapolation beyond the data range.

    Parameters
    ----------
    x : (n,) float64
        Monotonically increasing independent variable (chainage).
    y : (n,) float64
        Dependent variable values at x.
    """
    return PchipInterpolator(x, y, extrapolate=True)


def _drop_unobserved(
    x: np.ndarray,
    H: np.ndarray,
    W: np.ndarray,
    S: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Remove nodes that have no valid (non-NaN) observations in either H or W
    across all timesteps. Resets chainage to start from zero.

    If all nodes are invalid the original arrays are returned unchanged.
    """
    has_h = ~np.all(np.isnan(H), axis=1)
    has_w = ~np.all(np.isnan(W), axis=1)
    keep  = np.where(has_h & has_w)[0]

    if keep.size == 0:
        return x, H, W, S

    x_out = x[keep] - x[keep].min()
    return x_out, H[keep], W[keep], S[keep]


def _build_chainage(xobs: np.ndarray, dx: float) -> np.ndarray:
    """
    Build a uniform computational chainage starting at 0 and ending at or
    beyond the most upstream SWOT node.

    Parameters
    ----------
    xobs : (nobs,) float64
        SWOT node chainage values [m].
    dx : float
        Computational node spacing [m].
    """
    xmax = xobs.max()
    n    = math.ceil(xmax / dx) + 1
    return np.linspace(0.0, (n - 1) * dx, n)


def _estimate_bed_slope(S: np.ndarray, min_slope: float) -> np.ndarray:
    """
    Estimate the bed slope profile from the time-averaged SWOT water surface
    slope. Negative and zero slopes are excluded from the mean (they indicate
    noise or adverse-slope artefacts). Rows with no positive slopes fall back
    to the reach-mean positive slope, or min_slope if none exists.

    Parameters
    ----------
    S : (nobs, nt) float64
        Water surface slope observations. NaN where not observed.
    min_slope : float
        Minimum allowable slope [m/m].
    """
    nobs = S.shape[0]

    # Reach-mean from all positive, non-NaN values
    all_valid = S[np.isfinite(S) & (S > 0)]
    S_reach   = min_slope if all_valid.size == 0 else all_valid.mean()

    S0 = np.empty(nobs)
    for k in range(nobs):
        row_valid = S[k][np.isfinite(S[k]) & (S[k] > 0)]
        S0[k]     = S_reach if row_valid.size == 0 else row_valid.mean()

    return S0


def _obs_chainage(
    xobs:      np.ndarray,
    H:         np.ndarray,
    W:         np.ndarray,
    min_slope: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fill spatial gaps in WSE and width via PCHIP interpolation and enforce
    upstream-increasing monotonicity in WSE.

    Timesteps with fewer than 2 valid H observations are marked invalid.

    Parameters
    ----------
    xobs : (nobs,) float64
        Node chainage [m].
    H : (nobs, nt) float64
        WSE observations, NaN-sentineled.
    W : (nobs, nt) float64
        Width observations, NaN-sentineled.
    min_slope : float
        Minimum slope used to enforce WSE monotonicity [m/m].

    Returns
    -------
    Hout : (nobs, nt) float64
    Wout : (nobs, nt) float64
    valid : (nt,) bool
    """
    nobs, nt = H.shape
    Hout  = np.zeros((nobs, nt))
    Wout  = np.zeros((nobs, nt))
    valid = np.zeros(nt, dtype=bool)

    for t in range(nt):
        h = H[:, t]
        w = W[:, t]

        good_h = np.where(np.isfinite(h))[0]
        good_w = np.where(np.isfinite(w))[0]

        if good_h.size < 2:
            continue

        # Interpolate W (fall back to nearest if only one valid node)
        if good_w.size >= 2:
            wout = _pchip(xobs[good_w], w[good_w])(xobs)
        elif good_w.size == 1:
            wout = np.full(nobs, w[good_w[0]])
        else:
            wout = np.zeros(nobs)

        # Interpolate H
        hout = _pchip(xobs[good_h], h[good_h])(xobs)

        # Enforce upstream-increasing monotonicity
        for k in range(1, nobs):
            dx       = xobs[k] - xobs[k - 1]
            required = hout[k - 1] + min_slope * dx
            if hout[k] < required:
                hout[k] = required

        Hout[:, t] = hout
        Wout[:, t] = wout
        valid[t]   = True

    return Hout, Wout, valid


def _interpolate_to_chainage(
    xobs: np.ndarray,
    y:    np.ndarray,
    x:    np.ndarray,
) -> np.ndarray:
    """Interpolate a complete (no NaN) node-space vector to computational chainage."""
    return _pchip(xobs, y)(x)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def preprocess(
    xobs:      np.ndarray,
    H:         np.ndarray,
    W:         np.ndarray,
    S:         np.ndarray,
    dx:        float = 200.0,
    min_slope: float = 1e-5,
) -> SWOTReach:
    """
    Preprocess raw SWOT observations into a SWOTReach ready for GVF solving.

    Parameters
    ----------
    xobs : (nobs,) float64
        Node chainage, downstream → upstream [m]. Need not start at zero.
    H : (nobs, nt) float64
        WSE observations [m]. NaN where not observed.
    W : (nobs, nt) float64
        Width observations [m]. NaN where not observed.
    S : (nobs, nt) float64
        Water surface slope observations [m/m]. NaN where not observed.
    dx : float
        Computational node spacing [m]. Default 200 m.
    min_slope : float
        Minimum bed slope [m/m]. Default 1e-5.

    Returns
    -------
    SWOTReach
    """
    xobs = strip_mask(xobs)
    H    = strip_mask(H)
    W    = strip_mask(W)
    S    = strip_mask(S)
    xobs, H, W, S = _drop_unobserved(xobs, H, W, S)
    nobs, nt      = H.shape

    # Uniform computational chainage
    x  = _build_chainage(xobs, dx)
    nx = x.size

    # Bed slope interpolant
    S0obs = _estimate_bed_slope(S, min_slope)
    S0    = _interpolate_to_chainage(xobs, S0obs, x)
    S0    = np.maximum(S0, min_slope)
    S0_itp = _pchip(x, S0)

    # Cumulative bed elevation, anchored at z(0) = 0
    dz    = S0[:-1] * np.diff(x)
    z     = np.concatenate([[0.0], np.cumsum(dz)])
    z_itp = _pchip(x, z)

    # Bankfull width and WSE interpolants
    # nanmax along axis=1: maximum observed value at each node across time
    wbf_obs = np.nanmedian(W, axis=1)
    hbf_obs = np.nanmedian(H, axis=1)
    wbf_itp = _pchip(x, _interpolate_to_chainage(xobs, wbf_obs, x))
    hbf_itp = _pchip(x, _interpolate_to_chainage(xobs, hbf_obs, x))

    # Fill and enforce monotonicity in observed profiles
    H_f, W_f, valid = _obs_chainage(xobs, H, W, min_slope)

    # Interpolate valid timesteps to computational chainage
    H_comp = np.zeros((nx, nt))
    W_comp = np.zeros((nx, nt))
    for t in np.where(valid)[0]:
        H_comp[:, t] = _interpolate_to_chainage(xobs, H_f[:, t], x)
        W_comp[:, t] = _interpolate_to_chainage(xobs, W_f[:, t], x)

    # Minimum downstream WSE across all valid observations
    ds_h  = H[0]
    hmin  = float(np.nanmin(ds_h))

    return SWOTReach(
        obs   = SWOTObs(x=xobs, H=H, W=W),
        x     = x,
        H     = H_comp,
        W     = W_comp,
        valid = valid,
        S0    = S0_itp,
        wbf   = wbf_itp,
        hbf   = hbf_itp,
        z     = z_itp,
        hmin  = hmin,
        nx    = nx,
        nobs  = nobs,
        nt    = nt,
    )
