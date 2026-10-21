"""Error model of a DEM of Difference (DoD) built from stable terrain.

Inside stable areas the true change is zero, so every raw DoD value there is error. Spread is measured with the
normalised median absolute deviation (NMAD), robust to the outliers stable areas always contain (Höhle and Höhle,
2009). Spatially varying error is modelled from stable terrain as a function of slope and of the distance to the
nearest measured cell (after Hugonnet et al., 2022):

- slope:          σ_DoD(θ) = sqrt(a² + (b · tan θ)²), fitted on cells measured in both epochs;
- interpolation:  σ_interp(d) = sqrt(NMAD_bin(d)² − NMAD_measured²), added only where an epoch is interpolated;
- systematic:     σ_sys = NMAD of the per-area medians, a fully correlated error for volumes over a zone.

All values are 1σ in metres, slopes in degrees.

References: hohle2009robust, hugonnet2022stable (see references.bib).
"""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np
from scipy.optimize import least_squares

NMAD_FACTOR = 1.4826
SLOPE_BINS_DEG = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 60, 90)
DISTANCE_BINS_M = (0, 2, 4, 8, 16, 32, 64)


def nmad(values: np.ndarray) -> float:
    """1.4826 · median(|x − median(x)|), ignoring NaN; NaN for an empty input."""
    x = np.asarray(values, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(NMAD_FACTOR * np.median(np.abs(x - np.median(x))))


def slope_model(slope_deg: np.ndarray, a_m: float, b_m: float) -> np.ndarray:
    """σ_DoD(θ) = sqrt(a² + (b · tan θ)²): vertical error a plus a horizontal error b projected by the slope."""
    t = np.tan(np.radians(np.asarray(slope_deg, dtype=np.float64)))
    return np.sqrt(a_m**2 + (b_m * t) ** 2)


def fit_slope_model(dz: np.ndarray, slope_deg: np.ndarray, measured: np.ndarray,
                    bins: Sequence[float] = SLOPE_BINS_DEG, max_slope_deg: float = 45.0,
                    min_cells_per_bin: int = 200) -> dict[str, Any]:
    """Fits a and b to the NMAD of measured stable cells per slope bin (weighted least squares, weights √n).

    Bins steeper than ``max_slope_deg`` or with fewer than ``min_cells_per_bin`` cells are reported but not fitted.
    With a single usable bin (flat site), b = 0 and a = NMAD of that bin.
    """
    dz, slope_deg, measured = (np.asarray(v).ravel() for v in (dz, slope_deg, measured))
    keep = measured.astype(bool) & np.isfinite(dz) & np.isfinite(slope_deg)
    dz, slope_deg = dz[keep], slope_deg[keep]

    rows = []
    for lo, hi in zip(bins[:-1], bins[1:], strict=True):
        sel = (slope_deg >= lo) & (slope_deg < hi)
        n = int(sel.sum())
        row = {"bin": f"{lo:g}-{hi:g} deg", "lo_deg": float(lo), "hi_deg": float(hi), "n": n,
               "median_m": float(np.median(dz[sel])) if n else None,
               "nmad_m": nmad(dz[sel]) if n else None,
               "tan_slope": float(np.median(np.tan(np.radians(slope_deg[sel])))) if n else None,
               "used": hi <= max_slope_deg and n >= min_cells_per_bin}
        rows.append(row)
    used = [r for r in rows if r["used"]]
    if not used:
        raise ValueError("No slope bin can be fitted: add stable areas or lower min_cells_per_bin.")

    if len(used) == 1:
        a_m, b_m = float(used[0]["nmad_m"]), 0.0
    else:
        t = np.array([r["tan_slope"] for r in used])
        s = np.array([r["nmad_m"] for r in used])
        w = np.sqrt(np.array([r["n"] for r in used], dtype=np.float64))
        # Start from the linear fit of σ² = a² + b² t², then fit σ itself
        a2, b2 = np.linalg.lstsq(np.column_stack((np.ones_like(t), t**2)) * w[:, None], s**2 * w, rcond=None)[0]
        fit = least_squares(lambda p: w * (np.hypot(p[0], p[1] * t) - s),
                            x0=[np.sqrt(max(a2, 1e-12)), np.sqrt(max(b2, 1e-12))], bounds=([0, 0], [np.inf, np.inf]))
        a_m, b_m = (float(v) for v in fit.x)

    for r in rows:
        r["model_nmad_m"] = float(np.hypot(a_m, b_m * r["tan_slope"])) if r["n"] else None
    residuals = np.array([r["nmad_m"] - r["model_nmad_m"] for r in used])
    return {"a_m": a_m, "b_m": b_m, "max_slope_deg": float(max_slope_deg), "bins": rows,
            "rms_residual_m": float(np.sqrt(np.mean(residuals**2))), "cells_used": int(sum(r["n"] for r in used))}


def sigma_dod_raster(slope_deg: np.ndarray, a_m: float, b_m: float, max_slope_deg: float = 45.0):
    """Per-cell σ_DoD from the slope model; NaN above ``max_slope_deg``, where the model has no data support.

    Returns (sigma, unreliable): a NaN σ gives a NaN Level of Detection, so those cells never count as significant.
    """
    slope_deg = np.asarray(slope_deg, dtype=np.float64)
    sigma = slope_model(slope_deg, a_m, b_m)
    unreliable = np.isfinite(slope_deg) & (slope_deg > max_slope_deg)
    sigma[unreliable | ~np.isfinite(slope_deg)] = np.nan
    return sigma, unreliable


def fit_interpolation_model(dz: np.ndarray, distance_m: np.ndarray, measured: np.ndarray,
                            bins: Sequence[float] = DISTANCE_BINS_M, min_cells_per_bin: int = 200) -> dict[str, Any]:
    """Extra σ of interpolated stable cells per distance bin: sqrt(NMAD_bin² − NMAD_measured²).

    ``distance_m`` is the distance to the nearest measured cell (the larger of the two epochs). Short bins, or bins
    no noisier than the measured cells, use the value fitted on all interpolated stable cells (``constant_m``).
    """
    dz, dist, measured = (np.asarray(v).ravel() for v in (dz, distance_m, measured))
    meas = measured.astype(bool) & np.isfinite(dz)
    interp = ~measured.astype(bool) & np.isfinite(dz) & np.isfinite(dist)
    if not meas.any():
        raise ValueError("No measured stable cell.")
    nmad_measured = nmad(dz[meas])

    def extra(values: np.ndarray) -> float:
        total = nmad(values)
        return float(np.sqrt(total**2 - nmad_measured**2)) if total > nmad_measured else 0.0

    constant = extra(dz[interp]) if interp.any() else 0.0
    edges = list(bins) + [np.inf]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        sel = interp & (dist >= lo) & (dist < hi)
        n = int(sel.sum())
        own = n >= min_cells_per_bin and nmad(dz[sel]) > nmad_measured
        rows.append({"lo_m": float(lo), "hi_m": None if np.isinf(hi) else float(hi), "n": n,
                     "sigma_interp_m": extra(dz[sel]) if own else constant, "source": "bin" if own else "constant"})
    return {"nmad_measured_m": nmad_measured, "constant_m": constant, "bins": rows,
            "cells_interpolated": int(interp.sum())}


def interpolation_sigma(distance_m: np.ndarray, model: dict[str, Any]) -> np.ndarray:
    """σ_interp at the given distances (NaN where the distance is NaN)."""
    distance_m = np.asarray(distance_m, dtype=np.float64)
    out = np.zeros_like(distance_m)
    for row in model["bins"]:
        hi = np.inf if row["hi_m"] is None else row["hi_m"]
        out[(distance_m >= row["lo_m"]) & (distance_m < hi)] = row["sigma_interp_m"]
    out[~np.isfinite(distance_m)] = np.nan
    return out


def systematic_error(dz: np.ndarray, area_ids: np.ndarray, measured: np.ndarray) -> dict[str, Any]:
    """σ_sys = NMAD of the per-area medians of measured stable cells; global median bias of all of them."""
    dz, area_ids, measured = (np.asarray(v).ravel() for v in (dz, area_ids, measured))
    keep = measured.astype(bool) & np.isfinite(dz) & (area_ids > 0)
    areas = [{"area_id": int(a), "n": int((keep & (area_ids == a)).sum()),
              "median_m": float(np.median(dz[keep & (area_ids == a)]))} for a in np.unique(area_ids[keep])]
    if not areas:
        raise ValueError("No stable area with measured cells.")
    medians = np.array([a["median_m"] for a in areas])
    return {"sigma_sys_m": nmad(medians), "global_median_bias_m": float(np.median(dz[keep])), "areas": areas}


def epoch_sigma(sigma_dod: np.ndarray, interpolated: np.ndarray, distance_m: np.ndarray,
                interp_model: dict[str, Any] | None) -> np.ndarray:
    """σ of one epoch: σ_DoD / √2, plus σ_interp in quadrature where this epoch is interpolated.

    σ_DoD of the stable areas combines both surfaces and is split evenly, so on measured cells
    LoD = k · sqrt(σ0² + σ1²) = k · σ_DoD.
    """
    sigma = np.asarray(sigma_dod, dtype=np.float64) / np.sqrt(2.0)
    if interp_model is None:
        return sigma
    extra = interpolation_sigma(distance_m, interp_model)
    sel = np.asarray(interpolated, dtype=bool) & np.isfinite(extra)
    out = sigma.copy()
    out[sel] = np.hypot(sigma[sel], extra[sel])
    return out
