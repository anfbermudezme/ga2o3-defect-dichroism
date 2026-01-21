from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple


def _norm(s: str) -> str:
    return s.strip().lower().replace("-", "_").replace(" ", "")


_CASE_SYNONYMS = {
    # common aliases -> folder names as in your tree
    "pristine": "pristine",
    "bulk": "pristine",
    "octa": "octa",
    "octahedral": "octa",
    "tetra": "tetra",
    "tetrahedral": "tetra",
    "o1": "O_1",
    "o_1": "O_1",
    "o2": "O_2",
    "o_2": "O_2",
    "o3": "O_3",
    "o_3": "O_3",
    "octa1x2x2": "octa_1x2x2",
    "octa_1x2x2": "octa_1x2x2",
    "tetra1x2x2": "tetra_1x2x2",
    "tetra_1x2x2": "tetra_1x2x2",
}


def _resolve_case_dir(parent: Path, case: str) -> Path:
    """Find a case subdirectory under *parent*, with case-insensitive matching.

    Parameters
    ----------
    parent:
        Directory that contains case folders (e.g. <ROOT>/BANDS).
    case:
        Case name. Accepts synonyms like "o1" -> "O_1".

    Returns
    -------
    Path to the case directory.

    Raises
    ------
    FileNotFoundError if no suitable directory exists.
    """
    if not parent.is_dir():
        raise FileNotFoundError(f"Expected directory not found: {parent}")

    wanted = _CASE_SYNONYMS.get(_norm(case), case)

    # direct hit
    direct = parent / wanted
    if direct.is_dir():
        return direct

    # case-insensitive scan
    wanted_norm = _norm(wanted)
    for p in parent.iterdir():
        if p.is_dir() and _norm(p.name) == wanted_norm:
            return p

    # last-resort: accept case as given
    case_norm = _norm(case)
    for p in parent.iterdir():
        if p.is_dir() and _norm(p.name) == case_norm:
            return p

    raise FileNotFoundError(f"Case '{case}' not found under {parent}")


def _find_first_existing(dirpath: Path, candidates: Iterable[str]) -> Path:
    for name in candidates:
        p = dirpath / name
        if p.exists():
            return p
    # case-insensitive search
    lowered = {c.lower() for c in candidates}
    for p in dirpath.iterdir():
        if p.is_file() and p.name.lower() in lowered:
            return p
    raise FileNotFoundError(
        f"None of the following files exist in {dirpath}: {', '.join(candidates)}"
    )


@dataclass(frozen=True)
class Ga2O3Dataset:
    """Utility object that knows your on-disk folder layout.

    This is intentionally lightweight: it only resolves paths and does not parse data.
    """

    root: Path

    @classmethod
    def from_env(cls, env_var: str = "GA2O3_DATA", default: str = ".") -> "Ga2O3Dataset":
        root = Path(os.getenv(env_var, default)).expanduser().resolve()
        return cls(root=root)

    # --- top-level directories -------------------------------------------------
    @property
    def bands_root(self) -> Path:
        return self.root / "BANDS"

    @property
    def dos_root(self) -> Path:
        return self.root / "DOS"

    @property
    def cphf_root(self) -> Path:
        # your tree calls this CPKS_CPHF
        return self.root / "CPKS_CPHF"

    @property
    def anbd_root(self) -> Path:
        return self.root / "ANBD"

    @property
    def formation_energy_root(self) -> Path:
        return self.root / "formation_energy"

    # --- BANDS -----------------------------------------------------------------
    def bands_case_dir(self, case: str) -> Path:
        return _resolve_case_dir(self.bands_root, case)

    def band_files(self, case: str) -> Tuple[Path, Path]:
        """Return (alpha_band, beta_band) paths for a case."""
        d = self.bands_case_dir(case)
        alpha = _find_first_existing(d, ["ALPHA.BAND", "ALPHA.band"])
        beta = _find_first_existing(d, ["BETA.BAND", "BETA.band"])
        return alpha, beta

    # --- DOS -------------------------------------------------------------------
    def dos_spin_files(self, case: str) -> Tuple[Path, Path]:
        """Attempt to locate spin-resolved DOS files for a case.

        Heuristics (case-insensitive):
        - Prefer files that contain the case name *and* 'ALPHA' / 'BETA' in the filename.
        - If nothing matches, fall back to a single *.DOSS file that contains the case name.
        - If still nothing matches, fall back to the first *.DOSS file in the folder and return it for both.

        Returns
        -------
        (alpha_path, beta_path)
        """
        dos_dir = self.dos_root
        if not dos_dir.is_dir():
            raise FileNotFoundError(f"DOS directory not found: {dos_dir}")

        files = [p for p in dos_dir.iterdir() if p.is_file() and p.suffix.lower() == ".doss"]
        if not files:
            raise FileNotFoundError(f"No .DOSS files found under {dos_dir}")

        case_norm = _norm(case)

        def matches_case(p: Path) -> bool:
            return case_norm in _norm(p.name)

        alpha_candidates = [p for p in files if matches_case(p) and "alpha" in p.name.lower()]
        beta_candidates  = [p for p in files if matches_case(p) and "beta"  in p.name.lower()]

        if alpha_candidates and beta_candidates:
            alpha_candidates.sort()
            beta_candidates.sort()
            return alpha_candidates[0], beta_candidates[0]

        # fall back: any DOSS containing case name
        case_files = [p for p in files if matches_case(p)]
        if case_files:
            case_files.sort()
            return case_files[0], case_files[0]

        # last-resort: any .DOSS in folder
        files.sort()
        return files[0], files[0]

    # --- CPHF ------------------------------------------------------------------
    def cphf_case_dir(self, case: str) -> Path:
        return _resolve_case_dir(self.cphf_root, case)

    def cphf_outputs(self, case: str, pattern: str = "*.out") -> list[Path]:
        d = self.cphf_case_dir(case)
        files = sorted(d.glob(pattern))
        return files

