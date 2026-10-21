"""Ground points to a gap-filled 1 m DTM, with a record of which cells are measured and which are interpolated.

Gridding is inverse-distance weighting of the ground points within ``radius`` of each cell centre (power 1), the
same rule as PDAL's ``writers.gdal`` with ``output_type="idw"`` (Butler et al.); the default radius res/√2 is the
circle through the cell corners. Empty cells are then filled with GDAL's inverse-distance fill (rasterio
``fillnodata``) up to ``max_search_px`` cells away, and cells outside the survey footprint are masked.

The source mask (1 measured, 2 interpolated, 0 nodata) and the distance to the nearest measured cell feed the
interpolation term of the error model.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

SOURCE_NODATA, SOURCE_MEASURED, SOURCE_INTERPOLATED = 0, 1, 2


@dataclass(frozen=True)
class Grid:
    """North-up raster grid: ``x0`` / ``y_top`` are the outer edges of the top-left cell."""

    x0: float
    y_top: float
    res: float
    nx: int
    ny: int

    @classmethod
    def from_bounds(cls, xmin: float, ymin: float, xmax: float, ymax: float, res: float) -> "Grid":
        """Grid snapped outward to multiples of ``res``."""
        x0, y0 = np.floor(xmin / res) * res, np.floor(ymin / res) * res
        x1, y1 = np.ceil(xmax / res) * res, np.ceil(ymax / res) * res
        return cls(float(x0), float(y1), float(res), int(round((x1 - x0) / res)), int(round((y1 - y0) / res)))

    @property
    def shape(self) -> tuple[int, int]:
        return self.ny, self.nx

    @property
    def transform(self):
        from rasterio.transform import from_origin

        return from_origin(self.x0, self.y_top, self.res, self.res)


class IdwAccumulator:
    """Accumulates IDW sums chunk by chunk, so a cloud of 10⁸ points never has to sit in memory at once."""

    def __init__(self, grid: Grid, radius: float | None = None, power: float = 1.0):
        self.grid = grid
        self.radius = grid.res / np.sqrt(2.0) if radius is None else float(radius)
        if self.radius > grid.res:
            raise ValueError("radius > res is not supported (only the 4 nearest cell centres are visited).")
        self.power = power
        size = grid.nx * grid.ny
        self._wz = np.zeros(size)
        self._w = np.zeros(size)
        self._count = np.zeros(size, dtype=np.int64)

    def add(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> None:
        g = self.grid
        # Continuous column / row of each point, cell centres at integers
        cf = (np.asarray(x, dtype=np.float64) - g.x0) / g.res - 0.5
        rf = (g.y_top - np.asarray(y, dtype=np.float64)) / g.res - 0.5
        c0, r0 = np.floor(cf).astype(np.int64), np.floor(rf).astype(np.int64)
        z = np.asarray(z, dtype=np.float64)
        for dc in (0, 1):
            for dr in (0, 1):
                c, r = c0 + dc, r0 + dr
                d = np.hypot(cf - c, rf - r) * g.res
                ok = (d <= self.radius) & (c >= 0) & (c < g.nx) & (r >= 0) & (r < g.ny)
                idx = r[ok] * g.nx + c[ok]
                w = 1.0 / np.maximum(d[ok], 1e-6) ** self.power
                self._wz += np.bincount(idx, weights=w * z[ok], minlength=self._wz.size)
                self._w += np.bincount(idx, weights=w, minlength=self._w.size)
                self._count += np.bincount(idx, minlength=self._count.size)

    def result(self) -> tuple[np.ndarray, np.ndarray]:
        """(z, count): z is NaN where no point fell within the radius."""
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(self._w > 0, self._wz / self._w, np.nan)
        return z.reshape(self.grid.shape), self._count.reshape(self.grid.shape)


def occupancy_footprint(x: np.ndarray, y: np.ndarray, grid: Grid, cell_m: float = 5.0,
                        closing_cells: int = 3) -> np.ndarray:
    """Survey footprint on ``grid``: occupied coarse cells, closed and hole-filled, as a boolean mask."""
    f = int(round(cell_m / grid.res))
    nx, ny = -(-grid.nx // f), -(-grid.ny // f)
    c = ((np.asarray(x) - grid.x0) / (grid.res * f)).astype(np.int64)
    r = ((grid.y_top - np.asarray(y)) / (grid.res * f)).astype(np.int64)
    ok = (c >= 0) & (c < nx) & (r >= 0) & (r < ny)
    coarse = np.zeros((ny, nx), dtype=bool)
    coarse[r[ok], c[ok]] = True
    if closing_cells:
        coarse = ndimage.binary_closing(coarse, iterations=closing_cells, border_value=0)
    coarse = ndimage.binary_fill_holes(coarse)
    return np.repeat(np.repeat(coarse, f, axis=0), f, axis=1)[: grid.ny, : grid.nx]


def fill_dtm(z: np.ndarray, footprint: np.ndarray, max_search_px: int = 50):
    """Fills empty cells inside the footprint; returns (z_filled, source, distance_to_measured_m_in_cells).

    The distance is in cells (multiply by the resolution), NaN on measured cells and outside the footprint.
    """
    from rasterio.fill import fillnodata

    measured = np.isfinite(z)
    filled = fillnodata(np.where(measured, z, 0.0).astype(np.float32), mask=measured.astype(np.uint8),
                        max_search_distance=max_search_px).astype(np.float64)
    filled[~footprint] = np.nan
    source = np.full(z.shape, SOURCE_NODATA, dtype=np.uint8)
    source[np.isfinite(filled)] = SOURCE_INTERPOLATED
    source[measured & footprint] = SOURCE_MEASURED
    filled[source == SOURCE_NODATA] = np.nan
    distance = ndimage.distance_transform_edt(~measured).astype(np.float64)
    distance[source != SOURCE_INTERPOLATED] = np.nan
    return filled, source, distance


def slope_deg(z: np.ndarray, res: float) -> np.ndarray:
    """Slope in degrees by central differences (NaN propagates)."""
    gy, gx = np.gradient(z, res)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def hillshade(z: np.ndarray, res: float, azimuth_deg: float = 315.0, altitude_deg: float = 45.0) -> np.ndarray:
    gy, gx = np.gradient(z, res)
    slope = np.pi / 2.0 - np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    az, alt = np.radians(360.0 - azimuth_deg + 90.0), np.radians(altitude_deg)
    return np.clip(np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * np.cos(az - aspect), 0, 1)


def write_raster(path: str | Path, array: np.ndarray, grid: Grid, crs_wkt: str, nodata=np.nan) -> Path:
    """Single-band Cloud Optimized GeoTIFF on ``grid``."""
    import rasterio

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    dtype = "uint8" if array.dtype == np.uint8 else "float32"
    profile = {"driver": "COG", "height": grid.ny, "width": grid.nx, "count": 1, "dtype": dtype, "crs": crs_wkt,
               "transform": grid.transform, "nodata": 0 if dtype == "uint8" else nodata, "compress": "deflate"}
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array.astype(dtype), 1)
    return path


def read_raster(path: str | Path) -> np.ndarray:
    import rasterio

    with rasterio.open(path) as src:
        band = src.read(1, masked=True)
        return band.filled(np.nan if band.dtype.kind == "f" else 0)
