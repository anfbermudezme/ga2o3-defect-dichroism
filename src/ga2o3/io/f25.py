from __future__ import annotations

import re
from pathlib import Path
from typing import List, Sequence

import numpy as np


_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)[eEdD][-\+]?\d+")


def read_f25_maps(path: str | Path, *, max_maps: int = 2) -> list[np.ndarray]:
    """Read up to *max_maps* 2D maps from a CASTEP .f25 file.

    The CASTEP formatted grid file typically contains multiple MAPN blocks:

        MAPN  NX  NY
        <origin line>
        <lattice line>
        <NX*NY floating point values, possibly in D exponent>

    This function returns a list of arrays shaped (NY, NX).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    maps: list[np.ndarray] = []
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        # iterate line-by-line and detect MAPN
        for line in fh:
            m = re.search(r"MAPN\s+(\d+)\s+(\d+)", line)
            if not m:
                continue
            nx, ny = int(m.group(1)), int(m.group(2))

            # skip origin & lattice lines (if present)
            try:
                next(fh)
                next(fh)
            except StopIteration:
                break

            data: list[float] = []
            # read until we have nx*ny floats
            while len(data) < nx * ny:
                try:
                    ln = next(fh)
                except StopIteration:
                    break
                nums = _FLOAT_RE.findall(ln)
                if not nums:
                    continue
                for s in nums:
                    data.append(float(s.replace("D", "E").replace("d", "E")))

            if len(data) < nx * ny:
                # incomplete map; stop
                break

            arr = np.asarray(data, dtype=float).reshape(ny, nx)
            maps.append(arr)

            if len(maps) >= max_maps:
                break

    return maps
