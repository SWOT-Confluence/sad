"""
gvf.py

Gradually-Varied Flow (GVF) solver for SWOT river reaches.

Translated from gvf.jl. Uses JAX + Diffrax for the ODE solve so that
gvf_solve can be vmapped over ensemble members for efficient batch inference.

Design notes
------------
- The SciPy PCHIP interpolants in SWOTReach cannot be traced by JAX. They
  are evaluated on a fixed grid *before* the ODE solve and passed in as
  plain JAX arrays. The ODE then uses jnp.interp (linear) for in-solver
  lookups. For large low-gradient rivers the difference between PCHIP and
  linear interpolation within the solver is negligible.
- gvf_rhs is a pure JAX function and can be JIT-compiled and vmapped.
- gvf_solve returns NaN-filled arrays on solver failure rather than None,
  so it is vmap-safe.
"""

from __future__ import annotations

from typing import Optional

import jax
import jax.numpy as jnp
import diffrax
import numpy as np

from preprocess import SWOTReach

jax.config.update("jax_enable_x64", True)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Regularisation term in GVF denominator (1 - Fr² + ε).
# Prevents instability near critical flow. Has no physical effect for
# large low-gradient rivers where Fr << 1.
_GVF_EPS = 1e-6

_G = 9.806   # gravitational acceleration [m/s²]


# ---------------------------------------------------------------------------
# Cross-section geometry
# ---------------------------------------------------------------------------

def _area(y: jax.Array, Wb: jax.Array, Yb: jax.Array, r: jax.Array) -> jax.Array:
    """
    Flow area for a Dingman power-law cross section.

        A = Wb * (Ym / Yb)^(1/r) * y

    where Ym = (r+1)/r * y is the thalweg depth.

    Parameters
    ----------
    y  : mean flow depth [m]
    Wb : bankfull width [m]
    Yb : bankfull mean depth [m]
    r  : Dingman shape exponent
    """
    Ym = (r + 1.0) / r * y
    return Wb * (Ym / Yb) ** (1.0 / r) * y


# ---------------------------------------------------------------------------
# GVF ODE right-hand side
# ---------------------------------------------------------------------------

def _make_gvf_rhs(
    x_grid:   jax.Array,
    S0_grid:  jax.Array,
    wbf_grid: jax.Array,
    hbf_grid: jax.Array,
    z_grid:   jax.Array,
):
    """
    Close over pre-evaluated reach grids and return a Diffrax-compatible
    ODE term function.

    The returned function has signature  f(x, y, args)  where:
      x    : chainage (independent variable) [m]
      y    : mean flow depth [m]  (scalar)
      args : (Q, n, r, α, z0)

    Integration runs upstream (x increases from x[0] to x[-1]).

        dy/dx = -(S0(x)·α - Sf(y)) / (1 - Fr(y)² + ε)

    Manning friction slope:  Sf = (n·Q / (A · y^(2/3)))²
    Froude number:           Fr = Q / (A · √(g·y))

    Both use the wide-channel approximation R ≈ y, standard for large
    low-gradient rivers observed by SWOT.
    """
    def rhs(x: jax.Array, y: jax.Array, args) -> jax.Array:
        Q, n, r, α, z0 = args
        y  = jnp.maximum(y, 0.01)

        # Linear interpolation on pre-evaluated grids (JAX-traceable)
        Wb  = jnp.interp(x, x_grid, wbf_grid)
        hbf = jnp.interp(x, x_grid, hbf_grid)
        z_r = jnp.interp(x, x_grid, z_grid)

        # Bed elevation: reference profile scaled by α, anchored at z0
        zx  = z0 + z_r * α
        Yb  = jnp.maximum(hbf - zx, 0.01)
        S0  = jnp.interp(x, x_grid, S0_grid) * α

        A   = _area(y, Wb, Yb, r)
        Sf  = (n * Q / (A * y ** (2.0 / 3.0))) ** 2
        Fr  = Q / (A * jnp.sqrt(_G * y))

        return -(S0 - Sf) / (1.0 - Fr ** 2 + _GVF_EPS)

    return rhs


# ---------------------------------------------------------------------------
# Single solve
# ---------------------------------------------------------------------------

def gvf_solve(
    Q:     float,
    n:     float,
    r:     float,
    alpha: float,
    z0:    float,
    H_bc:  float,
    reach: SWOTReach,
    saveat: Optional[np.ndarray] = None,
) -> Optional[np.ndarray]:
    """
    Solve the GVF equation for a single parameter set.

    Parameters
    ----------
    Q     : discharge [m³/s]
    n     : Manning roughness coefficient
    r     : Dingman channel shape exponent
    alpha : slope correction factor (S_actual = S0 · alpha)
    z0    : downstream bed elevation [m]
    H_bc  : downstream boundary WSE [m]  (reach.H[0, t])
    reach : SWOTReach from preprocess.py
    saveat : chainage locations at which to return WSE [m].
             Defaults to reach.obs.x (aligns output with SWOT observations).

    Returns
    -------
    (nout,) float64 array of predicted WSE at saveat locations,
    or None if the solve fails or produces non-finite values.
    """
    if saveat is None:
        saveat = reach.obs.x

    try:
        result = _gvf_solve_jax(
            jnp.float64(Q),
            jnp.float64(n),
            jnp.float64(r),
            jnp.float64(alpha),
            jnp.float64(z0),
            jnp.float64(H_bc),
            jnp.asarray(reach.x,           dtype=jnp.float64),
            jnp.asarray(reach.S0(reach.x),  dtype=jnp.float64),
            jnp.asarray(reach.wbf(reach.x), dtype=jnp.float64),
            jnp.asarray(reach.hbf(reach.x), dtype=jnp.float64),
            jnp.asarray(reach.z(reach.x),   dtype=jnp.float64),
            jnp.asarray(saveat,             dtype=jnp.float64),
        )
    except Exception:
        return None

    result_np = np.asarray(result)
    if not np.all(np.isfinite(result_np)):
        return None
    return result_np


@jax.jit
def _gvf_solve_jax(
    Q:        jax.Array,
    n:        jax.Array,
    r:        jax.Array,
    alpha:    jax.Array,
    z0:       jax.Array,
    H_bc:     jax.Array,
    x_grid:   jax.Array,
    S0_grid:  jax.Array,
    wbf_grid: jax.Array,
    hbf_grid: jax.Array,
    z_grid:   jax.Array,
    saveat:   jax.Array,
) -> jax.Array:
    """
    JAX-traceable GVF solve. Returns WSE at saveat locations.
    Returns NaN-filled array on solver failure (safe for vmap).
    """
    rhs  = _make_gvf_rhs(x_grid, S0_grid, wbf_grid, hbf_grid, z_grid)
    term = diffrax.ODETerm(rhs)

    # Downstream boundary condition: convert WSE to mean depth
    y_bc = jnp.maximum((H_bc - z0) * r / (r + 1.0), 0.01)

    sol = diffrax.diffeqsolve(
        term,
        diffrax.Tsit5(),
        t0      = x_grid[0],
        t1      = x_grid[-1],
        dt0     = (x_grid[-1] - x_grid[0]) / 100.0,
        y0      = y_bc,
        args    = (Q, n, r, alpha, z0),
        saveat  = diffrax.SaveAt(ts=saveat),
        stepsize_controller = diffrax.PIDController(
            rtol=1e-3, atol=1e-3, jump_ts=None,
        ),
        max_steps = 10_000,
        throw     = False,   # return NaN rather than raising on failure
    )

    # Convert depth solution back to WSE
    # WSE = thalweg_depth + bed_elevation
    #     = y * (r+1)/r  +  z0 + z_ref(x) * alpha
    y_saveat = sol.ys                                      # mean depth at saveat
    z_saveat = z0 + jnp.interp(saveat, x_grid, z_grid) * alpha
    wse      = y_saveat * ((r + 1.0) / r) + z_saveat

    # NaN-fill if solver did not succeed
    success = sol.result == diffrax.RESULTS.successful
    return jnp.where(success, wse, jnp.full_like(wse, jnp.nan))


# ---------------------------------------------------------------------------
# Batched solve over an ensemble
# ---------------------------------------------------------------------------

def gvf_solve_ensemble(
    Q_ens:     np.ndarray,
    n_ens:     np.ndarray,
    r_ens:     np.ndarray,
    alpha_ens: np.ndarray,
    z0_ens:    np.ndarray,
    H_bc:      float,
    reach:     SWOTReach,
    saveat:    Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Solve GVF for a full ensemble in one vmapped call.

    Parameters
    ----------
    Q_ens, n_ens, r_ens, alpha_ens, z0_ens : (N,) float64
        Ensemble of parameter draws.
    H_bc : float
        Downstream boundary WSE [m], shared across ensemble.
    reach : SWOTReach
    saveat : (nout,) float64, optional
        Defaults to reach.obs.x.

    Returns
    -------
    (N, nout) float64 array of predicted WSE.
    NaN rows indicate solver failure for that ensemble member.
    """
    if saveat is None:
        saveat = reach.obs.x

    # Pre-evaluate SciPy interpolants on the computational grid once
    x_grid   = jnp.asarray(reach.x,            dtype=jnp.float64)
    S0_grid  = jnp.asarray(reach.S0(reach.x),  dtype=jnp.float64)
    wbf_grid = jnp.asarray(reach.wbf(reach.x), dtype=jnp.float64)
    hbf_grid = jnp.asarray(reach.hbf(reach.x), dtype=jnp.float64)
    z_grid   = jnp.asarray(reach.z(reach.x),   dtype=jnp.float64)
    saveat_j = jnp.asarray(saveat,             dtype=jnp.float64)
    H_bc_j   = jnp.float64(H_bc)

    Q_j     = jnp.asarray(Q_ens,     dtype=jnp.float64)
    n_j     = jnp.asarray(n_ens,     dtype=jnp.float64)
    r_j     = jnp.asarray(r_ens,     dtype=jnp.float64)
    alpha_j = jnp.asarray(alpha_ens, dtype=jnp.float64)
    z0_j    = jnp.asarray(z0_ens,    dtype=jnp.float64)

    result = _gvf_solve_ensemble_jax(
        Q_j, n_j, r_j, alpha_j, z0_j,
        H_bc_j,
        x_grid, S0_grid, wbf_grid, hbf_grid, z_grid,
        saveat_j,
    )
    return np.asarray(result)


@jax.jit
def _gvf_solve_ensemble_jax(
    Q_ens:     jax.Array,
    n_ens:     jax.Array,
    r_ens:     jax.Array,
    alpha_ens: jax.Array,
    z0_ens:    jax.Array,
    H_bc:      jax.Array,
    x_grid:    jax.Array,
    S0_grid:   jax.Array,
    wbf_grid:  jax.Array,
    hbf_grid:  jax.Array,
    z_grid:    jax.Array,
    saveat:    jax.Array,
) -> jax.Array:
    """vmap _gvf_solve_jax over the ensemble dimension."""
    def single(Q, n, r, alpha, z0):
        return _gvf_solve_jax(
            Q, n, r, alpha, z0, H_bc,
            x_grid, S0_grid, wbf_grid, hbf_grid, z_grid,
            saveat,
        )
    return jax.vmap(single)(Q_ens, n_ens, r_ens, alpha_ens, z0_ens)
