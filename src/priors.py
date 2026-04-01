"""
priors.py

Prior distributions for the SAD inference algorithm.

Translated from priors.jl. Uses scipy.stats distributions throughout.
All distributions expose a consistent .rvs() / .ppf() / .pdf() interface
via scipy.stats frozen distribution objects.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Optional

import numpy as np
from scipy.stats import (
    loguniform,
    lognorm,
    truncnorm,
    uniform,
)
from netCDF4 import Dataset

from preprocess import SWOTReach
from utils import strip_mask


# ---------------------------------------------------------------------------
# River planform classification
# ---------------------------------------------------------------------------

class River(IntEnum):
    """
    River planform classification used to set uninformative prior bounds on
    the Dingman shape exponent r.
    """
    braided      = 1
    sinuous      = 2
    more_sinuous = 3
    straight     = 4


# r prior bounds indexed by River enum value
_R_BOUNDS = {
    River.braided:      (0.5,  1.0),
    River.sinuous:      (1.0,  5.0),
    River.more_sinuous: (5.0,  10.0),
    River.straight:     (10.0, 20.0),
}


# ---------------------------------------------------------------------------
# Data structure
# ---------------------------------------------------------------------------

@dataclass
class SWOTPriors:
    """
    Prior distributions for all inferred parameters in the SAD algorithm.

    All fields are frozen scipy.stats distribution objects, exposing
    .rvs(), .ppf(), .pdf(), .logpdf() etc.

    Attributes
    ----------
    Qp : frozen distribution
        Discharge prior [m³/s]. Truncated LogNormal.
    np : frozen distribution
        Manning roughness coefficient prior. Uniform.
    rp : frozen distribution
        Dingman shape exponent prior. Truncated Normal or Uniform.
    zp : frozen distribution
        Downstream bed elevation prior [m]. Uniform.
    ap : frozen distribution
        Slope correction factor prior (centred on 1). LogNormal.
    """
    Qp: object
    np: object
    rp: object
    zp: object
    ap: object


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _truncated_lognormal(
    mu_log: float,
    sigma_log: float,
    lo: float,
    hi: float,
) -> object:
    """
    Truncated LogNormal distribution.

    Parameters
    ----------
    mu_log : float
        Mean of the underlying Normal (log-space).
    sigma_log : float
        Std dev of the underlying Normal (log-space).
    lo, hi : float
        Lower and upper truncation bounds in natural space.
    """
    # scipy parameterises lognorm as lognorm(s, scale=exp(mu))
    # truncation in log-space maps to Normal truncation
    a = (np.log(lo) - mu_log) / sigma_log
    b = (np.log(hi) - mu_log) / sigma_log
    # Use truncnorm on log-space then exponentiate via lognorm
    # Easier: use lognorm with explicit ppf-based rejection via
    # scipy's truncated approach
    # scipy has no built-in truncated lognormal, so we build it from
    # truncnorm on the log scale and wrap in a thin adapter
    return _TruncatedLogNormal(mu_log, sigma_log, lo, hi)


class _TruncatedLogNormal:
    """
    Minimal frozen-distribution adapter for a truncated LogNormal,
    mirroring the scipy frozen distribution interface (.rvs, .logpdf, .ppf).

    Internally uses a truncated Normal on the log scale.
    """
    def __init__(self, mu: float, sigma: float, lo: float, hi: float):
        self._mu    = mu
        self._sigma = sigma
        self._lo    = lo
        self._hi    = hi
        a = (np.log(lo) - mu) / sigma
        b = (np.log(hi) - mu) / sigma
        self._tn = truncnorm(a, b, loc=mu, scale=sigma)

    def rvs(self, size=None, random_state=None):
        return np.exp(self._tn.rvs(size=size, random_state=random_state))

    def logpdf(self, x):
        log_x = np.log(np.clip(x, 1e-300, None))
        # log p(x) = log p_normal(log x) - log x
        return self._tn.logpdf(log_x) - log_x

    def pdf(self, x):
        return np.exp(self.logpdf(x))

    def ppf(self, q):
        return np.exp(self._tn.ppf(q))

    @property
    def support(self):
        return (self._lo, self._hi)


def _truncated_normal(
    mu: float,
    sigma: float,
    lo: float,
    hi: float,
) -> object:
    """Truncated Normal distribution via scipy.stats.truncnorm."""
    a = (lo - mu) / sigma
    b = (hi - mu) / sigma
    return truncnorm(a, b, loc=mu, scale=sigma)


def _uniform(lo: float, hi: float) -> object:
    """Uniform distribution via scipy.stats.uniform."""
    return uniform(loc=lo, scale=hi - lo)


def _lognormal(mu_log: float, sigma_log: float) -> object:
    """
    Untruncated LogNormal via scipy.stats.lognorm.
    scipy parameterisation: lognorm(s=sigma, scale=exp(mu))
    """
    return lognorm(s=sigma_log, scale=np.exp(mu_log))


def _z0_prior(
    qwbm:  float,
    hmin:  float,
    reach: Optional[SWOTReach] = None,
) -> object:
    """
    Estimate the downstream bed elevation prior.

    If reach is provided, estimates depth from Manning scaling at mid-reach
    geometry. Otherwise falls back to a fixed depth estimate based on qwbm.

    Parameters
    ----------
    qwbm : float
        Mean discharge estimate [m³/s].
    hmin : float
        Minimum observed downstream WSE [m].
    reach : SWOTReach, optional
        Preprocessed reach. If provided, mid-reach slope and bankfull width
        are used to estimate depth via Manning's equation.

    Returns
    -------
    Uniform distribution over [z0_est - 3, z0_est + 3].
    """
    if reach is not None:
        x_mid  = reach.x[-1] / 2.0
        S_med  = float(reach.S0(x_mid))
        W_med  = float(reach.wbf(x_mid))
        n_est  = 0.035 if qwbm > 500.0 else 0.030
        S_use  = max(S_med, 1e-5)
        depth_est = (n_est * qwbm / (W_med * np.sqrt(S_use))) ** 0.6
        depth_est = float(np.clip(depth_est, 2.0, 20.0))
    else:
        if qwbm > 500.0:
            depth_est = 7.0
        elif qwbm > 100.0:
            depth_est = 5.0
        else:
            depth_est = 3.0

    z0_est = hmin - depth_est
    return _uniform(z0_est - 3.0, z0_est + 3.0)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def priors_from_sos(
    sosfile: str,
    hmin:    float,
    reachid: int,
) -> Optional[SWOTPriors]:
    """
    Derive prior distributions from the SoS (SWORD of Science) database.

    Parameters
    ----------
    sosfile : str
        Path to SoS NetCDF file.
    hmin : float
        Minimum observed downstream WSE [m]. Used as upper bound for z0 prior.
    reachid : int
        SWORD reach ID.

    Returns
    -------
    SWOTPriors, or None if the discharge prior cannot be constructed
    (missing mean_q in SoS).
    """
    with Dataset(sosfile) as f:
        gr = f.groups["reaches"]
        reach_ids = gr.variables["reach_id"][:]
        idx = int(np.where(reach_ids == reachid)[0][0])

        # --- Roughness ---
        g   = f.groups["gbpriors"].groups["reach"]
        n_l = float(np.exp(strip_mask(g.variables["lowerbound_logn"][idx])))
        n_u = float(np.exp(strip_mask(g.variables["upperbound_logn"][idx])))
        try:
            assert n_l < n_u
            np_ = _uniform(n_l, n_u)
        except Exception:
            np_ = _uniform(0.01, 0.10)

        # --- Channel shape ---
        r_m = float(np.exp(strip_mask(g.variables["logr_hat"][idx])))
        r_s = float(np.exp(strip_mask(g.variables["logr_sd"][idx])))
        r_l = float(np.exp(strip_mask(g.variables["lowerbound_logr"][idx])))
        r_u = float(np.exp(strip_mask(g.variables["upperbound_logr"][idx])))
        try:
            assert r_l < r_u and r_s > 0
            rp = _truncated_normal(r_m, r_s, r_l, r_u)
        except Exception:
            rp = _uniform(1.0, 10.0)

        # --- Discharge ---
        m   = f.groups["model"]
        q_m = m.variables["mean_q"][idx]
        q_u = float(strip_mask(m.variables["max_q"][idx]))
        q_l = float(strip_mask(m.variables["min_q"][idx]))

        if np.ma.is_masked(q_m) or q_m is None:
            return None

        q_m = float(strip_mask(q_m))
        qm  = np.log(q_m) - 2.0 ** 2 / 2.0
        if not np.isfinite(qm):
            qm = (q_u + q_l) / 2.0
        try:
            assert q_l < q_u
            Qp = _truncated_lognormal(qm, 2.0, q_l, q_u)
        except Exception:
            Qp = _truncated_lognormal(qm, 2.0, 0.1 * q_m, 20.0 * q_m)

        # --- Bed elevation ---
        z0_est = hmin - 5.0
        zp = _uniform(z0_est - 3.0, z0_est + 3.0)

        # --- Slope correction ---
        ap = _lognormal(0.0, 0.2)

        return SWOTPriors(Qp=Qp, np=np_, rp=rp, zp=zp, ap=ap)


def priors_from_qwbm(
    qwbm:  float,
    hmin:  float,
    river: River,
    reach: Optional[SWOTReach] = None,
) -> SWOTPriors:
    """
    Construct uninformative priors when SoS data are unavailable.

    Parameters
    ----------
    qwbm : float
        Mean discharge estimate [m³/s].
    hmin : float
        Minimum observed downstream WSE [m].
    river : River
        Planform classification — controls r prior bounds.
    reach : SWOTReach, optional
        If provided, z0 is estimated from reach geometry via Manning scaling
        rather than a fixed depth offset from hmin.

    Returns
    -------
    SWOTPriors
    """
    # Discharge
    qm = np.log(qwbm) - 2.0 ** 2 / 2.0
    Qp = _truncated_lognormal(qm, 2.0, 0.1 * qwbm, 20.0 * qwbm)

    # Roughness
    if qwbm > 500.0:
        n_lo = 0.025
    elif qwbm > 100.0:
        n_lo = 0.020
    else:
        n_lo = 0.015
    np_ = _uniform(n_lo, 0.07)

    # Shape exponent
    r_lo, r_hi = _R_BOUNDS[river]
    rp = _uniform(r_lo, r_hi)

    # Bed elevation
    zp = _z0_prior(qwbm, hmin, reach)

    # Slope correction
    ap = _lognormal(0.0, 0.2)

    return SWOTPriors(Qp=Qp, np=np_, rp=rp, zp=zp, ap=ap)
