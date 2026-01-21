from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import re

from ga2o3.constants import HARTREE_TO_EV


_META_KEYS = ("NEPTS", "NPROJ", "NSPIN")


def _parse_header_metadata(header_lines: list[str]) -> Dict[str, int]:
    text = " ".join(line.lstrip("#").strip() for line in header_lines)
    meta: Dict[str, int] = {}
    for key in _META_KEYS:
        m = re.search(rf"\b{key}\s+(\d+)", text)
        if m:
            meta[key.lower()] = int(m.group(1))
    return meta


@dataclass(frozen=True)
class DossData:
    """CRYSTAL density-of-states data loaded from a .DOSS file."""

    energies_hartree: np.ndarray  # (N,)
    dos: np.ndarray              # (N, M) states/Hartree/cell
    metadata: Dict[str, int]

    @property
    def energies_ev(self) -> np.ndarray:
        return self.energies_hartree * HARTREE_TO_EV

    @property
    def nspin(self) -> Optional[int]:
        return self.metadata.get("nspin")

    @property
    def nproj(self) -> Optional[int]:
        return self.metadata.get("nproj")

    @property
    def nepts(self) -> Optional[int]:
        return self.metadata.get("nepts")


def read_doss(path: str | Path) -> DossData:
    """Read a CRYSTAL .DOSS file.

    The reader is intentionally permissive:
    - skips blank lines and lines starting with '#', '@'
    - parses any remaining whitespace-separated floats

    Returns
    -------
    DossData
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    header_lines: list[str] = []
    energies: list[float] = []
    cols: list[list[float]] = []

    with p.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("#"):
                header_lines.append(s)
            if not s or s.startswith("#") or s.startswith("@"):
                continue
            try:
                vals = [float(x) for x in s.split()]
            except ValueError:
                continue
            if len(vals) < 2:
                continue
            energies.append(vals[0])
            cols.append(vals[1:])

    if not energies:
        raise ValueError(f"No DOS data found in {p}")

    e = np.asarray(energies, dtype=float)
    dos = np.asarray(cols, dtype=float)

    # sanity: consistent column count
    if dos.ndim != 2:
        raise ValueError(f"Unexpected DOS array shape in {p}: {dos.shape}")

    meta = _parse_header_metadata(header_lines)

    # If NEPTS is present, enforce length check (but do not hard-fail if mismatch; warn via ValueError?)
    nepts = meta.get("nepts")
    if nepts is not None and nepts != e.size:
        # keep going but store actual length as a hint
        meta["nepts_parsed"] = int(e.size)

    return DossData(energies_hartree=e, dos=dos, metadata=meta)
