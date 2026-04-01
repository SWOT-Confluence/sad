"""
utils.py

Shared utilities for the SAD inference pipeline.
"""

from __future__ import annotations

import numpy as np


def strip_mask(arr) -> np.ndarray:
    """
    Convert a masked array to a plain float64 array with NaN sentinels.

    netCDF4 returns masked arrays by default. Calling this at every
    NetCDF boundary ensures masked arrays never propagate into the
    pipeline where they can cause silent shape or broadcast errors.

    Parameters
    ----------
    arr : array-like or np.ma.MaskedArray

    Returns
    -------
    Plain float64 ndarray with masked values replaced by NaN.
    """
    if isinstance(arr, np.ma.MaskedArray):
        return np.where(arr.mask, np.nan, arr.data).astype(np.float64)
    return np.asarray(arr, dtype=np.float64)
