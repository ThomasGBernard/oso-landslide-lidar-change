"""Which flight days does a LAS/LAZ file hold? Histogram of GPS time by weekday and hour.

LAS 1.2 files with global-encoding bit 0 unset store GPS *week* time: seconds since the start of the GPS week
(Sunday 00:00). Week time alone cannot tell two weeks apart, so the days are matched against the vendor's list of
flights. For the 2014 WSDOT survey (Quantum Spatial report): 24 March (Mon), 2-3 April (Wed-Thu), 6-7 April
(Sun-Mon) and 9 April 2014 (Wed).

    python scripts/check_flight_dates.py data/raw/2014/*.laz
"""

from __future__ import annotations

import sys
from collections import Counter

import laspy
import numpy as np

DAYS = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"]
FLIGHTS_2014 = {"Mon": "24 Mar or 7 Apr", "Wed": "2 Apr or 9 Apr", "Thu": "3 Apr", "Sun": "6 Apr"}
PDT_OFFSET_H = -7          # local time in spring 2014 (UTC-7); GPS-UTC = 16 s, ignored at this resolution


def main(paths: list[str]) -> None:
    for path in paths:
        counts: Counter = Counter()
        with laspy.open(path) as f:
            if f.header.global_encoding.gps_time_type != laspy.header.GpsTimeType.WEEK_TIME:
                print(f"{path}: adjusted standard GPS time, convert with t + 1e9 instead")
                continue
            for chunk in f.chunk_iterator(5_000_000):
                t = np.asarray(chunk.gps_time) + PDT_OFFSET_H * 3600.0
                day = np.floor(t / 86400.0).astype(int) % 7
                hour = np.floor((t % 86400.0) / 3600.0).astype(int)
                counts.update(zip(day.tolist(), hour.tolist(), strict=True))
        total = sum(counts.values())
        print(f"\n{path}: {total:,} points (local time, UTC-7)")
        for (d, h), n in sorted(counts.items()):
            print(f"  {DAYS[d]} {h:02d}h  {n:>12,}  {100 * n / total:5.1f} %  {FLIGHTS_2014.get(DAYS[d], '')}")


if __name__ == "__main__":
    main(sys.argv[1:])
