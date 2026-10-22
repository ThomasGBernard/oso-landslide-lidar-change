"""Figures of the case study: change map, error model, budget against published values."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .uncertainty import slope_model

# Published values the budget is compared with (Mm³), see README "Comparison with published results"
PUBLISHED = {
    "source": {"value": 3.8, "low": 2.8, "high": 4.7,
               "label": "Iverson et al. (2015): 8.3 mobilised − ~4.5 stranded in the source"},
    "runout": {"value": 3.9, "low": None, "high": None,
               "label": "Anderson et al. (2017): valley-floor deposit (3.4 in July 2015 + 0.506 eroded)"},
}


def _save(fig, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def dod_map(significant: np.ndarray, shade: np.ndarray, extent, zones_gdf, path: str | Path,
            vlim_m: float = 30.0, title: str = "", water: np.ndarray | None = None,
            labels: tuple[str, str] = ("2013", "2014")) -> Path:
    """Significant change over a hillshade, zones outlined; colour range ±``vlim_m``.

    With ``water`` (DoD water classes): water in t1 in grey (left out of the budget), water in t0 only outlined in
    blue (kept, flagged).
    """
    fig, ax = plt.subplots(figsize=(8, 9))
    ax.imshow(shade, cmap="gray", extent=extent, interpolation="nearest")
    im = ax.imshow(significant, cmap="RdBu", norm=TwoSlopeNorm(0, -vlim_m, vlim_m), extent=extent,
                   interpolation="nearest", alpha=0.85)
    handles = []
    if water is not None:
        t1 = np.isin(water, (1, 3))
        if t1.any():
            ax.imshow(np.where(t1, 1.0, np.nan), cmap=ListedColormap(["#8c8c8c"]), extent=extent,
                      interpolation="nearest")
            handles.append(Patch(facecolor="#8c8c8c", label=f"Water in {labels[1]}: left out"))
        t0_only = water == 2
        if t0_only.any():
            ax.contour(t0_only.astype(float), levels=[0.5], colors="tab:blue", linewidths=0.8, extent=extent,
                       origin="upper")
            handles.append(Line2D([], [], color="tab:blue", label=f"Water in {labels[0]} only: kept, flagged"))
    if zones_gdf is not None:
        zones_gdf.boundary.plot(ax=ax, color="k", linewidth=1.0)
        for _, z in zones_gdf.iterrows():
            p = z.geometry.representative_point()
            ax.annotate(str(z["polygon_id"]), (p.x, p.y), ha="center", fontsize=9, weight="bold")
    cb = fig.colorbar(im, ax=ax, shrink=0.7, extend="both")
    cb.set_label("Elevation change 2014 − 2013 (m)\nred = erosion, blue = deposition")
    ax.set_xlabel("Easting (m), EPSG:32148")
    ax.set_ylabel("Northing (m)")
    ax.set_title(title)
    ax.ticklabel_format(style="plain", useOffset=False)
    if handles:
        ax.legend(handles=handles, loc="lower right", fontsize=8)
    return _save(fig, path)


def error_model(fit: dict, path: str | Path) -> Path:
    """NMAD of measured stable cells per slope bin with the fitted σ_DoD(θ) curve."""
    rows = [r for r in fit["bins"] if r["n"]]
    fig, ax = plt.subplots(figsize=(6, 4))
    used = [r for r in rows if r["used"]]
    other = [r for r in rows if not r["used"]]
    mid = lambda r: (r["lo_deg"] + r["hi_deg"]) / 2  # noqa: E731
    ax.plot([mid(r) for r in used], [r["nmad_m"] for r in used], "o", label="stable cells, fitted bins")
    if other:
        ax.plot([mid(r) for r in other], [r["nmad_m"] for r in other], "x", color="grey", label="not fitted")
    theta = np.linspace(0, fit["max_slope_deg"], 100)
    ax.plot(theta, slope_model(theta, fit["a_m"], fit["b_m"]), "-",
            label=f"σ = √(a² + (b·tanθ)²), a = {fit['a_m']:.3f} m, b = {fit['b_m']:.3f} m")
    ax.set_xlabel("Slope (°)")
    ax.set_ylabel("NMAD of raw DoD (m)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    return _save(fig, path)


def budget_vs_published(volumes, uncertainty, path: str | Path) -> Path:
    """Source loss and runout gain with their stated uncertainty, next to the published values."""
    zones = ["source", "runout"]
    totals = volumes[volumes.part == "total"].set_index("zone")
    unc = uncertainty.set_index("zone")
    ours = [abs(totals.loc[z, "net_m3"]) / 1e6 for z in zones]
    err = [unc.loc[z, "stated_uncertainty_m3"] / 1e6 for z in zones]
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(zones))
    ax.bar(x - 0.18, ours, 0.36, yerr=err, capsize=4, label="This study (significant DoD)")
    pub = [PUBLISHED[z]["value"] for z in zones]
    lo = [PUBLISHED[z]["value"] - PUBLISHED[z]["low"] if PUBLISHED[z]["low"] else 0 for z in zones]
    hi = [PUBLISHED[z]["high"] - PUBLISHED[z]["value"] if PUBLISHED[z]["high"] else 0 for z in zones]
    ax.bar(x + 0.18, pub, 0.36, yerr=[lo, hi], capsize=4, color="grey", label="Published")
    ax.set_xticks(x, ["Volume that left the source\n(net loss)", "Valley-floor deposit\n(net gain)"])
    ax.set_ylabel("Volume (Mm³)")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    return _save(fig, path)
