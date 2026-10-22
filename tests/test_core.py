"""Tests of the numerical core on synthetic data (numpy / scipy only, no lidar needed)."""

import numpy as np
import pytest

from oso_cd import dod, water
from oso_cd import uncertainty as unc
from oso_cd.dtm import Grid, IdwAccumulator, occupancy_footprint, slope_deg
from oso_cd.prepare import build_pipeline
from oso_cd.registration import apply_transform, voxel_reduce

RNG = np.random.default_rng(42)


def test_nmad_equals_std_for_normal_errors():
    x = RNG.normal(0.0, 0.1, 200_000)
    assert unc.nmad(x) == pytest.approx(0.1, rel=0.02)
    assert np.isnan(unc.nmad(np.array([np.nan])))


def test_slope_model_fit_recovers_a_and_b():
    a, b = 0.06, 0.25
    slope = RNG.uniform(0, 44.9, 300_000)
    dz = RNG.normal(0, unc.slope_model(slope, a, b))
    fit = unc.fit_slope_model(dz, slope, np.ones_like(dz, bool))
    assert fit["a_m"] == pytest.approx(a, rel=0.05)
    assert fit["b_m"] == pytest.approx(b, rel=0.05)
    assert all(not r["used"] for r in fit["bins"] if r["hi_deg"] > 45)


def test_sigma_raster_is_nan_above_max_slope():
    sigma, unreliable = unc.sigma_dod_raster(np.array([10.0, 50.0, np.nan]), 0.06, 0.25, 45.0)
    assert np.isfinite(sigma[0]) and np.isnan(sigma[1]) and np.isnan(sigma[2])
    assert unreliable.tolist() == [False, True, False]


def test_interpolation_model_extra_sigma():
    n = 50_000
    measured = np.r_[np.ones(n, bool), np.zeros(n, bool)]
    dz = np.r_[RNG.normal(0, 0.09, n), RNG.normal(0, np.hypot(0.09, 0.15), n)]
    dist = np.r_[np.full(n, np.nan), RNG.uniform(0, 1.9, n)]
    model = unc.fit_interpolation_model(dz, dist, measured)
    assert model["bins"][0]["source"] == "bin"
    assert model["bins"][0]["sigma_interp_m"] == pytest.approx(0.15, rel=0.05)
    assert unc.interpolation_sigma(np.array([1.0, np.nan]), model)[0] == pytest.approx(0.15, rel=0.05)


def test_systematic_error_is_spread_of_area_medians():
    ids = np.repeat([1, 2, 3, 4], 1000)
    offsets = np.array([0.0, 0.02, -0.02, 0.04])
    dz = offsets[ids - 1] + RNG.normal(0, 0.05, ids.size)
    out = unc.systematic_error(dz, ids, np.ones_like(dz, bool))
    assert len(out["areas"]) == 4
    assert out["sigma_sys_m"] == pytest.approx(unc.nmad([a["median_m"] for a in out["areas"]]))


def test_epoch_sigma_splits_dod_sigma_and_adds_interpolation():
    model = {"bins": [{"lo_m": 0.0, "hi_m": None, "sigma_interp_m": 0.2}]}
    s = unc.epoch_sigma(np.array([0.1, 0.1]), np.array([False, True]), np.array([np.nan, 3.0]), model)
    assert s[0] == pytest.approx(0.1 / np.sqrt(2))
    assert s[1] == pytest.approx(np.hypot(0.1 / np.sqrt(2), 0.2))


def test_dod_lod_and_zonal_budget():
    z0 = np.zeros((4, 4))
    z1 = np.zeros((4, 4))
    z1[0, :2] = -5.0          # erosion, 2 cells
    z1[3, 2:] = 3.0           # deposition, 2 cells
    z1[2, 0] = 0.05           # below the LoD
    sig = np.full((4, 4), 0.1)
    out = dod.compute_dod(z0, z1, sig, sig, k=1.96)
    assert np.isnan(out["significant"][2, 0])
    src = np.ones((4, 4), np.uint8)
    src[3, 3] = 2
    source = dod.dod_source(src, np.ones((4, 4), np.uint8), out["raw"])
    zones = {"all": np.ones((4, 4), bool)}
    table = dod.zonal_volumes(out["significant"], source, zones, cell_area_m2=1.0)
    total = table[(table.part == "total")].iloc[0]
    assert total.erosion_m3 == pytest.approx(10.0)
    assert total.deposition_m3 == pytest.approx(6.0)
    assert total.net_m3 == pytest.approx(-4.0)
    interp = table[(table.part == "interpolated")].iloc[0]
    assert interp.deposition_m3 == pytest.approx(3.0)
    u = dod.volume_uncertainty(table, sigma_sys_m=0.03, interp_sigma_m=0.19)
    assert u.iloc[0].systematic_m3 == pytest.approx(0.03 * 16)
    assert u.iloc[0].interpolation_conservative_m3 == pytest.approx(1.0)


def test_idw_grid_constant_surface_and_counts():
    grid = Grid.from_bounds(0, 0, 10, 10, 1.0)
    assert grid.shape == (10, 10)
    x, y = RNG.uniform(0, 10, 20_000), RNG.uniform(0, 10, 20_000)
    acc = IdwAccumulator(grid)
    half = x.size // 2
    acc.add(x[:half], y[:half], np.full(half, 7.0))      # two chunks give the same result as one
    acc.add(x[half:], y[half:], np.full(x.size - half, 7.0))
    z, count = acc.result()
    assert np.nanmax(np.abs(z - 7.0)) < 1e-9
    assert count.sum() > 0


def test_idw_point_on_cell_centre_dominates():
    grid = Grid.from_bounds(0, 0, 2, 2, 1.0)
    acc = IdwAccumulator(grid)
    acc.add(np.array([0.5, 0.9]), np.array([1.5, 1.1]), np.array([10.0, 0.0]))   # top-left cell centre is (0.5, 1.5)
    z, _ = acc.result()
    assert z[0, 0] == pytest.approx(10.0, abs=1e-3)


def test_footprint_fills_interior_holes():
    grid = Grid.from_bounds(0, 0, 100, 100, 1.0)
    xx, yy = np.meshgrid(np.arange(0.5, 100), np.arange(0.5, 100))
    ring = (np.hypot(xx - 50, yy - 50) > 15) & (np.hypot(xx - 50, yy - 50) < 40)
    fp = occupancy_footprint(xx[ring], yy[ring], grid, cell_m=5.0, closing_cells=1)
    assert fp[50, 50]                 # the hole in the middle is filled
    assert not fp[2, 2]               # the corners outside the ring stay outside


def test_slope_of_a_30_degree_plane():
    x = np.arange(50.0)
    z = np.tile(x * np.tan(np.radians(30)), (50, 1))
    assert np.allclose(slope_deg(z, 1.0), 30.0)


def test_voxel_reduce_keeps_one_point_per_voxel():
    xyz = RNG.uniform(0, 10, (5000, 3))
    idx = voxel_reduce(xyz, 2.0)
    keys = np.floor((xyz[idx] - xyz.min(axis=0)) / 2.0).astype(int)
    assert len({tuple(k) for k in keys}) == len(idx) <= 125


def test_apply_transform_about_reduction_point():
    r = np.array([425_000.0, 142_000.0, 200.0])
    angle = np.radians(0.3)
    affine = np.eye(4)
    affine[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    affine[:3, 3] = [0.5, -0.2, 0.1]
    xyz = r + RNG.uniform(-100, 100, (10, 3))
    moved = apply_transform(xyz, affine, r)
    assert np.allclose(apply_transform(r[None], affine, r), r + affine[:3, 3])
    back = (moved - affine[:3, 3] - r) @ affine[:3, :3] + r
    assert np.allclose(back, xyz)


def test_pdal_pipeline_overrides_crs_and_keeps_ground():
    p = build_pipeline(["a.laz", "b.laz"], "out.laz", "EPSG:2285+6360", "EPSG:32148+5703", classes=[2])
    types = [s["type"] for s in p["pipeline"]]
    assert types == ["readers.las", "readers.las", "filters.merge", "filters.reprojection", "filters.range",
                     "writers.las"]
    assert all(s["override_srs"] == "EPSG:2285+6360" for s in p["pipeline"][:2])
    assert p["pipeline"][4]["limits"] == "Classification[2:2]"


def test_water_points_counted_on_grid():
    grid = Grid.from_bounds(0, 0, 4, 3, 1.0)
    counts = water.count_on_grid(np.array([0.2, 0.8, 1.5, -1.0, 4.5]), np.array([2.5, 2.9, 1.5, 2.5, 2.5]), grid)
    assert counts[0, 0] == 2 and counts[1, 1] == 1 and counts.sum() == 3   # points outside the grid ignored


def test_water_mask_closes_gaps_fills_small_holes_keeps_islands():
    sparse = np.zeros((60, 60), int)
    sparse[10:50:2, 10:50:2] = 1          # a return every other cell
    assert water.mask_from_counts(sparse, 1.0, closing_m=2.0, max_hole_m2=0, buffer_m=0)[11:48, 11:48].all()
    lake = np.zeros((60, 60), int)
    lake[10:50, 10:50] = 1
    lake[20:22, 20:22] = 0                # 4 m2 dropout: filled
    lake[30:40, 30:40] = 0                # 100 m2 island: kept
    m = water.mask_from_counts(lake, 1.0, closing_m=0, max_hole_m2=50, buffer_m=0)
    assert m[20:22, 20:22].all() and not m[32:38, 32:38].any()
    b = water.mask_from_counts(lake, 1.0, closing_m=0, max_hole_m2=50, buffer_m=2.0)
    assert b[8, 20] and not b[7, 20]      # 2 m buffer


def test_t1_water_leaves_the_budget_and_t0_water_is_flagged():
    z0 = np.zeros((4, 4))
    z1 = np.zeros((4, 4))
    z1[0, :] = 1.0                        # water level rise in t1 (row 0 is water in 2014)
    z1[3, :2] = 2.0                       # deposit in a former channel (water in t0 only)
    sig = np.full((4, 4), 0.1)
    change = dod.compute_dod(z0, z1, sig, sig)
    w0, w1 = np.zeros((4, 4), bool), np.zeros((4, 4), bool)
    w1[0, :] = True
    w0[3, :2] = True
    classes = water.dod_water(w0, w1, np.isfinite(change["raw"]))
    kept, left_out = water.leave_out_t1_water(change["significant"], classes)
    table = dod.zonal_volumes(kept, np.ones((4, 4), np.uint8), {"all": np.ones((4, 4), bool)}, 1.0,
                              water=classes, left_out=left_out).set_index("part")
    assert table.loc["total", "deposition_m3"] == pytest.approx(4.0)            # the t1 water is not counted
    assert table.loc["water_t1_left_out", "deposition_m3"] == pytest.approx(4.0)
    assert table.loc["water_t0_only", "deposition_m3"] == pytest.approx(4.0)    # kept, inside the total


def test_pdal_pipeline_for_water_points():
    p = build_pipeline(["a.laz"], "w.laz", "EPSG:2285+6360", "EPSG:32148+5703", classes=[9])
    assert p["pipeline"][2]["limits"] == "Classification[9:9]"
