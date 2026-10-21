"""Step 0: fix the CRS units of the delivered tiles and keep the ground class.

Both epochs are labelled EPSG:32148 (NAD83 / Washington North, metres) in their GeoKeys but hold US survey feet in
X, Y and Z. Read unchecked, every distance would be in feet labelled as metres: vertical change ×3.28, volumes
×35.3, and a CRS check on the EPSG code would still pass. The PDAL pipeline built here therefore *overrides* the
header CRS with the true feet-based one (EPSG:2285 + NAVD88 ftUS, EPSG:6360) and reprojects to metres
(EPSG:32148 + NAVD88 m, EPSG:5703): a unit conversion on the same datum, no datum shift.

The pipeline runs with the ``pdal`` command-line tool (Butler et al.), so python-pdal is not needed.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Sequence


def build_pipeline(inputs: Sequence[str | Path], output: str | Path, in_srs: str, out_srs: str,
                   classes: Sequence[int] | None = (2,)) -> dict:
    """PDAL pipeline: read (CRS overridden) -> merge -> reproject -> keep ``classes`` -> LAS 1.4 point format 6."""
    if not inputs:
        raise ValueError("No input file.")
    stages: list[dict] = [{"type": "readers.las", "filename": str(Path(p).as_posix()), "override_srs": in_srs}
                          for p in inputs]
    if len(inputs) > 1:
        stages.append({"type": "filters.merge"})
    stages.append({"type": "filters.reprojection", "in_srs": in_srs, "out_srs": out_srs})
    if classes:
        stages.append({"type": "filters.range",
                       "limits": ",".join(f"Classification[{c}:{c}]" for c in classes)})
    stages.append({"type": "writers.las", "filename": str(Path(output).as_posix()), "minor_version": 4,
                   "dataformat_id": 6, "scale_x": 0.01, "scale_y": 0.01, "scale_z": 0.01,
                   "offset_x": "auto", "offset_y": "auto", "offset_z": "auto", "a_srs": out_srs})
    return {"pipeline": stages}


def run_pipeline(pipeline: dict, record: str | Path | None = None) -> None:
    """Runs ``pdal pipeline --stdin``; writes the pipeline JSON next to the output as a record if asked."""
    exe = shutil.which("pdal")
    if exe is None:
        raise FileNotFoundError("The 'pdal' command was not found: install PDAL from conda-forge (environment.yml).")
    text = json.dumps(pipeline, indent=2)
    if record is not None:
        Path(record).write_text(text + "\n", encoding="utf-8")
    subprocess.run([exe, "pipeline", "--stdin"], input=text, text=True, check=True)
