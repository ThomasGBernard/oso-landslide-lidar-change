"""The workflow, step by step, driven by ``config/oso.toml``.

Each step writes its outputs and is skipped when they already exist (pass ``force=True`` to redo it). Large
intermediate products go to ``data/work/`` (git-ignored); the small results that the repository publishes go to
``results/`` and ``figures/``.

    python -m oso_cd.pipeline config/oso.toml          # whole workflow
"""

from __future__ import annotations

import json
import sys
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

from . import dod as dodm
from . import dtm as dtmm
from . import plots, prepare, registration
from . import uncertainty as unc


def load_config(path: str | Path) -> dict:
    """Reads the TOML file; every path in [paths] is resolved against the repository root (parent of config/)."""
    path = Path(path).resolve()
    cfg = tomllib.loads(path.read_text(encoding="utf-8"))
    root = path.parent.parent
    cfg["root"] = root
    for key, value in cfg["paths"].items():
        cfg["paths"][key] = [root / v for v in value] if isinstance(value, list) else root / value
    for key in ("work", "results", "figures"):
        cfg["paths"][key].mkdir(parents=True, exist_ok=True)
    return cfg


def _labels(cfg: dict) -> tuple[str, str]:
    return cfg["labels"]["t0"], cfg["labels"]["t1"]


def _expand(patterns: list[Path]) -> list[Path]:
    files: list[Path] = []
    for p in patterns:
        files += sorted(p.parent.glob(p.name)) if any(ch in p.name for ch in "*?[") else [p]
    missing = [f for f in files if not f.exists()]
    if not files or missing:
        raise FileNotFoundError(f"Input lidar not found: {missing or patterns}. See data/README.md.")
    return files


def _write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, default=float) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------------------------- step 0
def step_prepare(cfg: dict, force: bool = False) -> dict[str, Path]:
    """Feet -> metres (CRS override + reprojection), ground class only, one LAZ per epoch."""
    work, out = cfg["paths"]["work"], {}
    for label, key in zip(_labels(cfg), ("raw_t0", "raw_t1"), strict=True):
        target = work / f"{label}_ground_m.laz"
        if force or not target.exists():
            pipe = prepare.build_pipeline(_expand(cfg["paths"][key]), target, cfg["crs"]["source"],
                                          cfg["crs"]["target"], classes=cfg["registration"]["classes"])
            prepare.run_pipeline(pipe, record=work / f"{label}_prepare_pipeline.json")
        out[label] = target
    return out


# --------------------------------------------------------------------------------------------- step 1
def step_register(cfg: dict, ground: dict[str, Path], force: bool = False) -> Path:
    t0, t1 = _labels(cfg)
    target = cfg["paths"]["work"] / f"{t1}_ground_m_registered.laz"
    record_path = cfg["paths"]["results"] / "registration.json"
    if force or not target.exists():
        p = cfg["registration"]
        params = {k: p[k] for k in ("voxel_size", "initial_distance_threshold", "level_of_detection",
                                    "supervoxel_resolution", "min_svp_num")}
        _write_json(record_path, registration.register_file(ground[t0], ground[t1], target, params))
    return target


# --------------------------------------------------------------------------------------------- step 2
def _bounds(path: Path) -> np.ndarray:
    import laspy

    with laspy.open(str(path)) as f:
        return np.r_[f.header.mins[:2], f.header.maxs[:2]]


def step_dtm(cfg: dict, clouds: dict[str, Path], force: bool = False) -> dict:
    """1 m IDW DTM per epoch on a common grid, gap-filled inside the common footprint."""
    import laspy
    import rasterio
    from pyproj import CRS

    d, work = cfg["dtm"], cfg["paths"]["work"]
    labels = _labels(cfg)
    outputs = {lab: {k: work / f"dtm_{lab}{s}.tif" for k, s in
                     (("z", ""), ("source", "_source"), ("distance", "_interp_distance_m"))} for lab in labels}
    if not force and all(p.exists() for o in outputs.values() for p in o.values()):
        with rasterio.open(outputs[labels[0]]["z"]) as src:
            t = src.transform
            grid = dtmm.Grid(t.c, t.f, t.a, src.width, src.height)
        return {"grid": grid, "outputs": outputs}

    b0, b1 = _bounds(clouds[labels[0]]), _bounds(clouds[labels[1]])
    common = np.r_[np.maximum(b0[:2], b1[:2]), np.minimum(b0[2:], b1[2:])]
    grid = dtmm.Grid.from_bounds(*common, res=d["resolution_m"])
    crs_wkt = CRS.from_user_input(cfg["crs"]["target"]).to_wkt()

    gridded, footprints = {}, []
    for lab in labels:
        acc = dtmm.IdwAccumulator(grid, radius=d.get("idw_radius_m"), power=d.get("idw_power", 1.0))
        xs, ys = [], []
        with laspy.open(str(clouds[lab])) as f:
            for chunk in f.chunk_iterator(5_000_000):
                x, y, z = np.asarray(chunk.x), np.asarray(chunk.y), np.asarray(chunk.z)
                acc.add(x, y, z)
                keep = slice(None, None, 10)   # every 10th point is plenty for a 5 m occupancy footprint
                xs.append(x[keep])
                ys.append(y[keep])
        gridded[lab] = acc.result()
        footprints.append(dtmm.occupancy_footprint(np.concatenate(xs), np.concatenate(ys), grid,
                                                   cell_m=d["footprint_cell_m"],
                                                   closing_cells=d["footprint_closing_cells"]))
    footprint = footprints[0] & footprints[1]

    stats = {}
    for lab in labels:
        z, count = gridded[lab]
        filled, source, distance = dtmm.fill_dtm(z, footprint, d["fill_max_search_px"])
        dtmm.write_raster(outputs[lab]["z"], filled, grid, crs_wkt)
        dtmm.write_raster(outputs[lab]["source"], source, grid, crs_wkt)
        dtmm.write_raster(outputs[lab]["distance"], distance * grid.res, grid, crs_wkt)
        inside = source > 0
        stats[lab] = {"cells": int(inside.sum()), "measured_pct": float(100 * (source == 1).sum() / inside.sum()),
                      "mean_points_per_measured_cell": float(count[source == 1].mean()),
                      "max_interp_distance_m": float(np.nanmax(distance) * grid.res)}
    _write_json(cfg["paths"]["results"] / "dtm_summary.json",
                {"grid": {"x0": grid.x0, "y_top": grid.y_top, "res": grid.res, "nx": grid.nx, "ny": grid.ny},
                 "epochs": stats})
    return {"grid": grid, "outputs": outputs}


# --------------------------------------------------------------------------------------------- steps 3-4
def _rasterize(gdf, column: str | None, grid: dtmm.Grid) -> np.ndarray:
    from rasterio.features import rasterize

    values = gdf[column] if column and column in gdf else range(1, len(gdf) + 1)
    return rasterize(list(zip(gdf.geometry, values, strict=True)), out_shape=grid.shape,
                     transform=grid.transform, fill=0, dtype="int32")


def step_change(cfg: dict, dtms: dict) -> dict:
    """Error model on stable areas, LoD-thresholded DoD, zonal budget with its uncertainty."""
    import geopandas as gpd
    from pyproj import CRS

    u, labels, grid = cfg["uncertainty"], _labels(cfg), dtms["grid"]
    res_dir, work = cfg["paths"]["results"], cfg["paths"]["work"]
    o = dtms["outputs"]
    z0, z1 = (dtmm.read_raster(o[lab]["z"]) for lab in labels)
    s0, s1 = (dtmm.read_raster(o[lab]["source"]) for lab in labels)
    d0, d1 = (dtmm.read_raster(o[lab]["distance"]) for lab in labels)
    raw = z1 - z0
    slope = dtmm.slope_deg((z0 + z1) / 2.0, grid.res)   # slope of the mean of the two DTMs
    measured = (s0 == dtmm.SOURCE_MEASURED) & (s1 == dtmm.SOURCE_MEASURED)
    dist = np.fmax(np.where(s0 == dtmm.SOURCE_INTERPOLATED, d0, np.nan),
                   np.where(s1 == dtmm.SOURCE_INTERPOLATED, d1, np.nan))

    horizontal = CRS.from_user_input(cfg["crs"]["target"]).sub_crs_list[0]
    stable = gpd.read_file(cfg["paths"]["stable_areas"]).to_crs(horizontal)
    area_id = _rasterize(stable, "polygon_id", grid)
    st = (area_id > 0) & np.isfinite(raw) & np.isfinite(slope)

    slope_fit = unc.fit_slope_model(raw[st], slope[st], measured[st], bins=u["slope_bins_deg"],
                                    max_slope_deg=u["max_slope_deg"], min_cells_per_bin=u["min_cells_per_bin"])
    interp_fit = unc.fit_interpolation_model(raw[st], dist[st], measured[st], bins=u["interp_distance_bins_m"],
                                             min_cells_per_bin=u["min_cells_per_bin"])
    sys_err = unc.systematic_error(raw[st], area_id[st], measured[st])
    sigma_dod, unreliable = unc.sigma_dod_raster(slope, slope_fit["a_m"], slope_fit["b_m"], u["max_slope_deg"])
    sigma_dod[~np.isfinite(raw)] = np.nan
    sig0 = unc.epoch_sigma(sigma_dod, s0 == dtmm.SOURCE_INTERPOLATED, d0, interp_fit)
    sig1 = unc.epoch_sigma(sigma_dod, s1 == dtmm.SOURCE_INTERPOLATED, d1, interp_fit)
    change = dodm.compute_dod(z0, z1, sig0, sig1, k=u["lod_k"])
    source = dodm.dod_source(s0, s1, change["raw"])

    stable_measured = raw[st & measured]
    stable_interp = raw[st & ~measured]
    error_model = {"slope_model": slope_fit, "interpolation_model": interp_fit, "systematic": sys_err,
                   "stable_measured_nmad_m": unc.nmad(stable_measured),
                   "stable_measured_median_m": float(np.median(stable_measured)),
                   "stable_interpolated_nmad_m": unc.nmad(stable_interp),
                   "unreliable_area_share_pct": float(100 * unreliable[np.isfinite(raw)].mean()),
                   "lod_median_m": float(np.nanmedian(change["lod"]))}
    _write_json(res_dir / "error_model.json", error_model)

    zones_gdf = gpd.read_file(cfg["paths"]["zones"]).to_crs(horizontal)
    zone_masks = {str(z.polygon_id): _rasterize(zones_gdf.iloc[[i]], None, grid) > 0
                  for i, z in enumerate(zones_gdf.itertuples())}
    zone_areas = {str(z.polygon_id): float(z.geometry.area) for z in zones_gdf.itertuples()}
    volumes = dodm.zonal_volumes(change["significant"], source, zone_masks, grid.res**2, zone_areas)
    volumes.to_csv(res_dir / "zonal_volumes.csv", index=False, float_format="%.1f")
    uncertainty = dodm.volume_uncertainty(volumes, sys_err["sigma_sys_m"], error_model["stable_interpolated_nmad_m"],
                                          u.get("conservative_interp_bias_m", 1.0))
    uncertainty.to_csv(res_dir / "volume_uncertainty.csv", index=False, float_format="%.1f")

    crs_wkt = CRS.from_user_input(cfg["crs"]["target"]).to_wkt()
    dtmm.write_raster(work / "dod_significant.tif", change["significant"], grid, crs_wkt)
    dtmm.write_raster(work / "dod_lod.tif", change["lod"], grid, crs_wkt)
    return {"grid": grid, "z1": z1, "significant": change["significant"], "zones": zones_gdf,
            "volumes": volumes, "uncertainty": uncertainty, "error_model": error_model}


# --------------------------------------------------------------------------------------------- step 5
def step_validate(cfg: dict, volumes: pd.DataFrame) -> dict:
    """Compares the zonal net volumes with the reference run of the in-house toolbox.

    ``rel_tol`` is one number, or a table with one tolerance per zone; ``notes`` explains a zone-specific tolerance
    and is copied into ``results/validation.json``.
    """
    v = cfg["validation"]
    ref = pd.read_csv(cfg["paths"]["reference_volumes"]).set_index("zone_id")
    ours = volumes[volumes.part == "total"].set_index("zone")
    rows, ok = [], True
    for zone in ref.index:
        r, o = float(ref.loc[zone, "net_volume_m3"]), float(ours.loc[str(zone), "net_m3"])
        rel = abs(o - r) / abs(r)
        tol = v["rel_tol"][str(zone)] if isinstance(v["rel_tol"], dict) else v["rel_tol"]
        ok &= rel <= tol
        rows.append({"zone": zone, "reference_net_m3": r, "this_net_m3": o, "rel_diff": rel, "rel_tol": tol,
                     "note": v.get("notes", {}).get(str(zone), "")})
    result = {"passed": bool(ok), "zones": rows}
    _write_json(cfg["paths"]["results"] / "validation.json", result)
    return result


def step_figures(cfg: dict, change: dict) -> list[Path]:
    fig, grid = cfg["paths"]["figures"], change["grid"]
    extent = (grid.x0, grid.x0 + grid.nx * grid.res, grid.y_top - grid.ny * grid.res, grid.y_top)
    shade = dtmm.hillshade(change["z1"], grid.res)
    t0, t1 = _labels(cfg)
    return [
        plots.dod_map(change["significant"], shade, extent, change["zones"], fig / "dod_map.png",
                      title=f"Significant elevation change {t0}–{t1} (95 % LoD)"),
        plots.error_model(change["error_model"]["slope_model"], fig / "error_model.png"),
        plots.budget_vs_published(change["volumes"], change["uncertainty"], fig / "budget_vs_published.png"),
    ]


def run(config: str | Path, force: bool = False) -> dict:
    cfg = load_config(config)
    t0, t1 = _labels(cfg)
    ground = step_prepare(cfg, force)
    registered = step_register(cfg, ground, force)
    dtms = step_dtm(cfg, {t0: ground[t0], t1: registered}, force)
    change = step_change(cfg, dtms)
    validation = step_validate(cfg, change["volumes"])
    step_figures(cfg, change)
    print(change["volumes"][change["volumes"].part == "total"].to_string(index=False))
    print("Validation against the toolbox reference run:", "PASSED" if validation["passed"] else "FAILED")
    return {"change": change, "validation": validation}


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "config/oso.toml", force="--force" in sys.argv)
