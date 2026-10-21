"""DEM of Difference, Level of Detection and zonal sediment budget.

dDEM = Z_t1 − Z_t0 (positive = deposition). A cell is significant when |dDEM| > LoD = k · sqrt(σ0² + σ1²), with
per-cell σ of each epoch from ``uncertainty.epoch_sigma`` (Brasington et al., 2003; Lane et al., 2003;
Wheaton et al., 2010). Volumes are summed over significant cells only and split by DTM source: measured in both
epochs, or interpolated in at least one.

References: brasington2003sediment, lane2003braided, wheaton2010dod (see references.bib).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

SOURCE_NODATA, SOURCE_MEASURED, SOURCE_INTERPOLATED = 0, 1, 2


def compute_dod(z0: np.ndarray, z1: np.ndarray, sigma0: np.ndarray, sigma1: np.ndarray,
                k: float = 1.96) -> dict[str, np.ndarray]:
    """Raw dDEM, LoD and significant dDEM (NaN where not significant or where σ is NaN)."""
    raw = np.asarray(z1, dtype=np.float64) - np.asarray(z0, dtype=np.float64)
    lod = k * np.hypot(sigma0, sigma1)
    with np.errstate(invalid="ignore"):
        significant = np.where(np.isfinite(raw) & (np.abs(raw) > lod), raw, np.nan)
    return {"raw": raw, "lod": lod, "significant": significant}


def dod_source(source0: np.ndarray, source1: np.ndarray, raw: np.ndarray) -> np.ndarray:
    """1 = measured in both epochs, 2 = interpolated in at least one, 0 = nodata."""
    out = np.full(raw.shape, SOURCE_NODATA, dtype=np.uint8)
    valid = np.isfinite(raw)
    out[valid] = SOURCE_INTERPOLATED
    out[valid & (source0 == SOURCE_MEASURED) & (source1 == SOURCE_MEASURED)] = SOURCE_MEASURED
    return out


def _volumes(values: np.ndarray, cell_area_m2: float) -> dict[str, float]:
    values = values[np.isfinite(values)]
    erosion = float(-values[values < 0].sum() * cell_area_m2)
    deposition = float(values[values > 0].sum() * cell_area_m2)
    return {"erosion_m3": erosion, "deposition_m3": deposition, "net_m3": deposition - erosion,
            "erosion_area_m2": float((values < 0).sum() * cell_area_m2),
            "deposition_area_m2": float((values > 0).sum() * cell_area_m2)}


def zonal_volumes(significant: np.ndarray, source: np.ndarray, zones: dict[str, np.ndarray],
                  cell_area_m2: float, zone_areas_m2: dict[str, float] | None = None) -> pd.DataFrame:
    """One row per zone and part (total, measured, interpolated) of the significant dDEM."""
    rows = []
    for name, mask in zones.items():
        parts = {"total": mask, "measured": mask & (source == SOURCE_MEASURED),
                 "interpolated": mask & (source == SOURCE_INTERPOLATED)}
        for part, sel in parts.items():
            rows.append({"zone": name, "part": part,
                         "zone_area_m2": (zone_areas_m2 or {}).get(name, float(mask.sum() * cell_area_m2)),
                         **_volumes(significant[sel], cell_area_m2)})
    return pd.DataFrame(rows)


def volume_uncertainty(table: pd.DataFrame, sigma_sys_m: float, interp_sigma_m: float,
                       conservative_bias_m: float = 1.0) -> pd.DataFrame:
    """Uncertainty terms of each zone's net volume (m³), stated next to the volumes rather than folded into one σ.

    - systematic: σ_sys × zone area (fully correlated within a zone);
    - interpolation, stable-area level: σ_interp (NMAD of interpolated stable cells) × interpolated active area;
    - interpolation, conservative: a bias of ``conservative_bias_m`` of one sign on every interpolated active cell.
    Random error is negligible over 10⁵–10⁶ uncorrelated cells and is not listed.
    """
    rows = []
    for zone, group in table.groupby("zone", sort=False):
        total = group[group.part == "total"].iloc[0]
        interp = group[group.part == "interpolated"].iloc[0]
        interp_area = interp.erosion_area_m2 + interp.deposition_area_m2
        systematic = sigma_sys_m * total.zone_area_m2
        conservative = conservative_bias_m * interp_area
        rows.append({"zone": zone, "net_m3": total.net_m3, "systematic_m3": systematic,
                     "interpolation_stable_m3": interp_sigma_m * interp_area,
                     "interpolation_conservative_m3": conservative,
                     "stated_uncertainty_m3": systematic + conservative,
                     "interpolated_share_of_net_pct": 100.0 * interp.net_m3 / total.net_m3 if total.net_m3 else np.nan})
    return pd.DataFrame(rows)
