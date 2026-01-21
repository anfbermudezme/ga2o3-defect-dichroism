from pathlib import Path
import numpy as np

from ga2o3.io.cphf import parse_eps_components, add_absorption_columns, apply_energy_correction
from ga2o3.constants import HC_EVN_M


def test_parse_eps_components_finds_two_blocks():
    p = Path(__file__).parent / "data" / "sample_cphf.out"
    df = parse_eps_components(p, round_wavelength=False)

    assert len(df) == 2
    assert set(df.columns) >= {"wavelength_nm", "energy_eV", "eps_xx", "eps_yy", "eps_zz"}

    # wavelength from hc/E
    assert np.isclose(df.loc[0, "wavelength_nm"], HC_EVN_M / df.loc[0, "energy_eV"])


def test_absorption_and_correction_runs():
    p = Path(__file__).parent / "data" / "sample_cphf.out"
    df = parse_eps_components(p, round_wavelength=False)
    df = add_absorption_columns(df)

    assert "alpha_iso_cm^-1" in df.columns
    assert np.all(df["alpha_iso_cm^-1"].to_numpy() > 0)

    df2 = apply_energy_correction(df, mode="scissor", scissor_shift_ev=1.0)
    assert "wavelength_orig_nm" in df2.columns
    assert np.all(df2["energy_eV"].to_numpy() > df2["energy_orig_eV"].to_numpy())
