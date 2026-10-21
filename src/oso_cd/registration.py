"""Step 1: align the post-event ground cloud to the pre-event one on stable terrain.

Both clouds are thinned to one point per voxel (the point closest to the voxel centroid), moved to a local frame
centred on a reduction point r, and registered with py4dgeo's stable-area ICP (Anders et al.): supervoxels whose
distance stays below a level of detection are taken as stable and used for a point-to-plane ICP (Besl and McKay,
1992; Yang and Schwieger, 2022). Working in the local frame keeps the supervoxel fit precise at projected
coordinates of ~10⁵–10⁶ m. The transformation is x' = R (x − r) + t + r.

References: besl1992icp, yang2022supervoxel, anders2026py4dgeo (see references.bib).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def voxel_reduce(xyz: np.ndarray, voxel_size: float) -> np.ndarray:
    """Indices of one point per voxel: the point closest to the centroid of its voxel's points."""
    xyz = np.asarray(xyz, dtype=np.float64)
    ijk = np.floor((xyz - xyz.min(axis=0)) / voxel_size).astype(np.int64)
    n = ijk.max(axis=0) + 1
    keys = ijk[:, 0] + n[0] * (ijk[:, 1] + n[1] * ijk[:, 2])   # one integer per voxel
    _, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
    inverse = inverse.ravel()
    centroids = np.column_stack([np.bincount(inverse, weights=xyz[:, i]) for i in range(3)]) / counts[:, None]
    d2 = ((xyz - centroids[inverse]) ** 2).sum(axis=1)
    order = np.lexsort((d2, inverse))
    first = np.r_[0, np.flatnonzero(np.diff(inverse[order])) + 1]
    return np.sort(order[first])


def apply_transform(xyz: np.ndarray, affine: np.ndarray, reduction_point: np.ndarray) -> np.ndarray:
    """x' = R (x − r) + t + r for a 4×4 (or 3×4) affine matrix [R | t]."""
    affine = np.asarray(affine, dtype=np.float64)
    r = np.asarray(reduction_point, dtype=np.float64).reshape(3)
    return (np.asarray(xyz, dtype=np.float64) - r) @ affine[:3, :3].T + affine[:3, 3] + r


def stable_area_icp(reference_xyz: np.ndarray, moving_xyz: np.ndarray, reduction_point: np.ndarray,
                    voxel_size: float = 1.0, initial_distance_threshold: float = 2.0,
                    level_of_detection: float = 0.2, supervoxel_resolution: float = 15.0,
                    min_svp_num: int = 5) -> np.ndarray:
    """4×4 affine matrix that moves ``moving_xyz`` onto ``reference_xyz`` (local frame around ``reduction_point``)."""
    import py4dgeo

    r = np.asarray(reduction_point, dtype=np.float64).reshape(3)

    def local_epoch(xyz: np.ndarray):
        sub = np.ascontiguousarray(xyz[voxel_reduce(xyz, voxel_size)] - r, dtype=np.float64)
        epoch = py4dgeo.Epoch(sub)
        epoch.build_kdtree()
        epoch.normals_attachment(epoch.calculate_normals())
        return epoch

    transformation = py4dgeo.icp_with_stable_areas(
        reference_epoch=local_epoch(reference_xyz), epoch=local_epoch(moving_xyz),
        initial_distance_threshold=initial_distance_threshold, level_of_detection=level_of_detection,
        reference_supervoxel_resolution=supervoxel_resolution, supervoxel_resolution=supervoxel_resolution,
        min_svp_num=min_svp_num, reduction_point=None,
    )
    affine = np.eye(4)
    affine[:3, :] = np.asarray(transformation.affine_transformation, dtype=np.float64)[:3, :]
    return affine


def read_xyz(path: str | Path) -> np.ndarray:
    import laspy

    las = laspy.read(str(path))
    return np.column_stack((las.x, las.y, las.z)).astype(np.float64)


def register_file(reference: str | Path, moving: str | Path, output: str | Path, params: dict) -> dict:
    """Registers ``moving`` onto ``reference`` and writes the moved cloud; returns the transformation record."""
    import laspy

    reference_xyz, moving_las = read_xyz(reference), laspy.read(str(moving))
    moving_xyz = np.column_stack((moving_las.x, moving_las.y, moving_las.z)).astype(np.float64)
    r = (reference_xyz.min(axis=0) + reference_xyz.max(axis=0)) / 2.0   # bounding-box centre of the reference
    affine = stable_area_icp(reference_xyz, moving_xyz, r, **params)
    moved = apply_transform(moving_xyz, affine, r)
    moving_las.x, moving_las.y, moving_las.z = moved[:, 0], moved[:, 1], moved[:, 2]
    moving_las.write(str(output))
    record = {"reference": Path(reference).name, "moving": Path(moving).name, "reduction_point_m": r.tolist(),
              "affine": affine.tolist(), "translation_m": affine[:3, 3].tolist(), "parameters": params}
    Path(output).with_suffix(".transform.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record
