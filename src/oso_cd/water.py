"""Water masks: where a DTM shows water, not ground.

A ground DTM has no ground under water: lidar returns from a water surface are classed as water (ASPRS class 9) or
not at all, and the gap-filled DTM bridges the lake or the river from its banks. The DoD then counts the water level
as ground change. Each epoch gets a water mask on the DTM grid, from the first source available:

1. the points of the water classes of that epoch (the 2014 delivery classifies water), counted per cell, closed to
   bridge the gaps between sparse water returns, small holes filled and a buffer added for the shore cells;
2. otherwise a polygon drawn for that epoch (the 2013 delivery has no water class: its river is drawn).

The DoD then keeps the two cases apart (``dod_water``): where t1 is water the cell leaves the budget, and the change
it would have added is reported separately; where only t0 was water (a channel since buried by the deposit) the cell
stays in the budget and is flagged, because the change there is underestimated by the t0 water depth.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

DRY, T1_ONLY, T0_ONLY, BOTH = 0, 1, 2, 3   # classes of the DoD water raster


def _disk(radius_px: int) -> np.ndarray:
    yy, xx = np.mgrid[-radius_px:radius_px + 1, -radius_px:radius_px + 1]
    return xx**2 + yy**2 <= radius_px**2


def count_on_grid(x: np.ndarray, y: np.ndarray, grid) -> np.ndarray:
    """Number of points per cell of ``grid`` (``dtm.Grid``); points outside the grid are ignored."""
    c = np.floor((np.asarray(x, dtype=np.float64) - grid.x0) / grid.res).astype(np.int64)
    r = np.floor((grid.y_top - np.asarray(y, dtype=np.float64)) / grid.res).astype(np.int64)
    ok = (c >= 0) & (c < grid.nx) & (r >= 0) & (r < grid.ny)
    return np.bincount(r[ok] * grid.nx + c[ok], minlength=grid.nx * grid.ny).reshape(grid.shape)


def buffer_mask(water: np.ndarray, res: float, buffer_m: float, valid: np.ndarray | None = None) -> np.ndarray:
    """Dilates ``water`` by ``buffer_m`` and keeps the ``valid`` cells."""
    r = int(round(buffer_m / res))
    if r > 0 and water.any():
        water = ndimage.binary_dilation(water, structure=_disk(r))
    return water & valid if valid is not None else water


def mask_from_counts(counts: np.ndarray, res: float, valid: np.ndarray | None = None, closing_m: float = 2.0,
                     max_hole_m2: float = 50.0, buffer_m: float = 1.0) -> np.ndarray:
    """Water mask from water points per cell: closing, holes up to ``max_hole_m2`` filled (islands kept), buffer."""
    water = np.asarray(counts) > 0
    r = int(round(closing_m / res))
    if r > 0 and water.any():
        padded = np.pad(water, r, mode="edge")   # so the closing does not erode water at the grid edge
        water = ndimage.binary_closing(padded, structure=_disk(r))[r:-r, r:-r]
    if max_hole_m2 > 0 and water.any():
        labels, n = ndimage.label(ndimage.binary_fill_holes(water) & ~water)
        if n:
            small = np.flatnonzero(np.bincount(labels.ravel())[1:] * res**2 <= max_hole_m2) + 1
            water |= np.isin(labels, small)
    return buffer_mask(water, res, buffer_m, valid)


def dod_water(water0: np.ndarray, water1: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Water class of each DoD cell: 0 dry, 1 water in t1 only, 2 water in t0 only, 3 water in both."""
    out = np.full(valid.shape, DRY, dtype=np.uint8)
    out[water1 & ~water0] = T1_ONLY
    out[water0 & ~water1] = T0_ONLY
    out[water0 & water1] = BOTH
    out[~valid] = DRY
    return out


def leave_out_t1_water(significant: np.ndarray, water: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Splits the significant dDEM: (kept, left out where t1 is water)."""
    t1 = np.isin(water, (T1_ONLY, BOTH))
    return np.where(t1, np.nan, significant), np.where(t1, significant, np.nan)
