from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Literal, Optional

import numpy as np
import pandas as pd

from ga2o3.constants import HC_EVN_M, NM_TO_CM


_float = r"([-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?)"
_row_re = re.compile(
    rf"^\s*{_float}\s+([XYZ]{{2}})\s+{_float}\s+{_float}\s+{_float}\s+{_float}(?:\s+{_float})?"
)


def sqrt_physical(z: complex) -> complex:
    """Complex sqrt with a physical branch (enforce Im(sqrt)>=0)."""
    w = np.sqrt(z)
    return w if w.imag >= 0 else -w


def parse_eps_components(path: str | Path, *, round_wavelength: bool = True) -> pd.DataFrame:
    """Parse ALL dielectric-table blocks in a CRYSTAL CPHF output file.

    The parser searches for blocks whose header line contains:
        "FR.(eV)" and "COMP." and "EPSILON"

    It then expects rows that match the regex defined in `_row_re`, and collects
    XX/YY/ZZ complex epsilon values for each energy.

    Returns
    -------
    DataFrame with columns:
      wavelength_nm, energy_eV, eps_xx, eps_yy, eps_zz, source
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()

    records: list[dict] = []
    for i, line in enumerate(lines):
        if ("FR.(eV)" in line) and ("COMP." in line) and ("EPSILON" in line):
            comps: dict[str, complex] = {}
            energy_ev: Optional[float] = None

            for j in range(i + 1, min(i + 250, len(lines))):
                m = _row_re.match(lines[j])
                if not m:
                    if lines[j].strip() == "":
                        break
                    continue

                energy_ev = float(m.group(1))
                tag = m.group(2)
                eps_re = float(m.group(5))
                eps_im = float(m.group(6))
                comps[tag] = complex(eps_re, eps_im)

            if energy_ev is None:
                continue
            if not all(k in comps for k in ("XX", "YY", "ZZ")):
                continue

            wl_nm = HC_EVN_M / energy_ev
            if round_wavelength:
                wl_nm = float(int(round(wl_nm)))

            records.append(
                {
                    "wavelength_nm": wl_nm,
                    "energy_eV": float(energy_ev),
                    "eps_xx": comps["XX"],
                    "eps_yy": comps["YY"],
                    "eps_zz": comps["ZZ"],
                    "source": p.name,
                }
            )

    if not records:
        return pd.DataFrame(
            columns=["wavelength_nm", "energy_eV", "eps_xx", "eps_yy", "eps_zz", "source"]
        )

    return pd.DataFrame(records).sort_values("wavelength_nm").reset_index(drop=True)


def add_absorption_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Compute n, k and absorption coefficient alpha for XX/YY/ZZ and isotropic.

    Absorption coefficient:
        alpha [cm^-1] = 4*pi*k / lambda[cm]

    The wavelength is taken from df["wavelength_nm"]. If you later apply an energy correction,
    recompute alpha using :func:`apply_energy_correction`.
    """
    if df.empty:
        return df.copy()

    out = df.copy()

    lam_cm = out["wavelength_nm"].to_numpy(dtype=float) * NM_TO_CM

    for comp in ("xx", "yy", "zz"):
        ntilde = out[f"eps_{comp}"].apply(sqrt_physical)
        out[f"n_{comp}"] = ntilde.apply(np.real).astype(float)
        out[f"k_{comp}"] = ntilde.apply(np.imag).astype(float)
        out[f"alpha_{comp}_cm^-1"] = 4 * math.pi * out[f"k_{comp}"].to_numpy(dtype=float) / lam_cm

    out["alpha_iso_cm^-1"] = out[["alpha_xx_cm^-1", "alpha_yy_cm^-1", "alpha_zz_cm^-1"]].mean(axis=1)
    return out


CorrectionMode = Literal["none", "scissor", "ratio"]


def apply_energy_correction(
    df: pd.DataFrame,
    *,
    mode: CorrectionMode = "none",
    scissor_shift_ev: float = 0.0,
    ratio_factor: float = 1.0,
    extra_shift_ev: float = 0.0,
) -> pd.DataFrame:
    """Correct the energy/wavelength axis and recompute alpha using existing k columns.

    Parameters
    ----------
    mode:
        - "none": no correction
        - "scissor": additive correction: E' = E + scissor_shift_ev
        - "ratio": multiplicative correction: E' = E * ratio_factor + extra_shift_ev
    scissor_shift_ev:
        Only used when mode="scissor".
    ratio_factor, extra_shift_ev:
        Only used when mode="ratio".

    Notes
    -----
    This function assumes `add_absorption_columns(df)` has already been called (so that
    k_xx/k_yy/k_zz exist). If not, it will raise a KeyError.
    """
    if df.empty:
        return df.copy()

    out = df.copy()
    out["wavelength_orig_nm"] = out.get("wavelength_orig_nm", out["wavelength_nm"])
    out["energy_orig_eV"] = out.get("energy_orig_eV", out["energy_eV"])

    if mode == "none":
        return out

    E = out["energy_orig_eV"].to_numpy(dtype=float)

    if mode == "scissor":
        E_new = E + float(scissor_shift_ev)
    elif mode == "ratio":
        E_new = E * float(ratio_factor) + float(extra_shift_ev)
    else:
        raise ValueError(f"Unknown correction mode: {mode}")

    # guard against non-physical zeros
    E_new = np.clip(E_new, 1e-9, None)

    out["energy_eV"] = E_new
    out["wavelength_nm"] = HC_EVN_M / out["energy_eV"]

    lam_cm = out["wavelength_nm"].to_numpy(dtype=float) * NM_TO_CM
    for comp in ("xx", "yy", "zz"):
        out[f"alpha_{comp}_cm^-1"] = 4 * math.pi * out[f"k_{comp}"].to_numpy(dtype=float) / lam_cm
    out["alpha_iso_cm^-1"] = out[["alpha_xx_cm^-1", "alpha_yy_cm^-1", "alpha_zz_cm^-1"]].mean(axis=1)

    return out
