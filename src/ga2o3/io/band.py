from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

from ga2o3.constants import HARTREE_TO_EV


@dataclass(frozen=True)
class BandData:
    """CRYSTAL band-structure data loaded from a .BAND file."""

    k: np.ndarray                     # (N,) k-coordinate along the path
    energies_hartree: np.ndarray      # (N, NBANDS) energies in Hartree
    xaxis_ticks: Tuple[float, ...] = ()  # high-symmetry x positions if present

    @property
    def energies_ev(self) -> np.ndarray:
        return self.energies_hartree * HARTREE_TO_EV


def _parse_xaxis_ticks(lines: Sequence[str]) -> Tuple[float, ...]:
    """Parse '@ XAXIS TICK i, value' lines if present."""
    coords_map: dict[int, float] = {}
    for line in lines:
        if line.startswith("@ XAXIS TICK") and "," in line and "SPEC" not in line:
            try:
                parts = line.split(",")
                idx = int(parts[0].replace("@ XAXIS TICK", "").strip())
                val = float(parts[1].strip())
                coords_map[idx] = val
            except Exception:
                continue
    if not coords_map:
        return ()
    return tuple(coords_map[i] for i in sorted(coords_map))


def read_band(path: str | Path) -> BandData:
    """Read a CRYSTAL .BAND file.

    Notes
    -----
    The file usually contains:
    - comment/header lines starting with '#' or '@'
    - numeric lines: k  E1  E2  ...  (energies in Hartree)
    - optional x-axis tick markers: '@ XAXIS TICK i, x'

    Returns
    -------
    BandData
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    ticks = _parse_xaxis_ticks(lines)

    data: list[list[float]] = []
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("@"):
            continue
        try:
            data.append([float(x) for x in s.split()])
        except ValueError:
            continue

    if not data:
        raise ValueError(f"No numeric band data found in {p}")

    arr = np.asarray(data, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(f"Unexpected band data shape in {p}: {arr.shape}")

    k = arr[:, 0]
    energies = arr[:, 1:]
    return BandData(k=k, energies_hartree=energies, xaxis_ticks=ticks)
