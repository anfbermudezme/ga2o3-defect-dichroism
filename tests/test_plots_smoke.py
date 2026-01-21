from pathlib import Path

from ga2o3.dataset import Ga2O3Dataset
from ga2o3.plot.bands import plot_bands_only
from ga2o3.plot.dos import plot_spin_resolved_dos
from ga2o3.plot.density import plot_charge_density_f25
from ga2o3.plot.cphf import plot_absorption_overlays


def test_plot_bands_only_smoke():
    data_dir = Path(__file__).parent / "data"
    fig, ax = plot_bands_only(data_dir / "sample.BAND", data_dir / "sample.BAND", band_legend=False)
    assert fig is not None
    assert ax is not None


def test_plot_spin_resolved_dos_smoke():
    data_dir = Path(__file__).parent / "data"
    fig, ax = plot_spin_resolved_dos(
        data_dir / "sample_alpha.DOSS",
        data_dir / "sample_beta.DOSS",
        drop_last_projection=False,
    )
    assert fig is not None
    assert ax is not None


def test_plot_charge_density_smoke(tmp_path: Path):
    data_dir = Path(__file__).parent / "data"
    out = tmp_path / "density.png"
    fig, ax = plot_charge_density_f25(data_dir / "sample.f25", out_path=out, show_colorbar=False)
    assert out.exists()
    assert fig is not None
    assert ax is not None


def test_plot_absorption_overlays_smoke(tmp_path: Path):
    # build minimal dataset with CPKS_CPHF/pristine/*.out
    root = tmp_path
    cdir = root / "CPKS_CPHF" / "pristine"
    cdir.mkdir(parents=True)
    sample = Path(__file__).parent / "data" / "sample_cphf.out"
    (cdir / "GaO_DISCRET_PRISTINE_CPHF_425.out").write_text(sample.read_text(encoding="utf-8"), encoding="utf-8")

    ds = Ga2O3Dataset(root=root)
    figs = plot_absorption_overlays(ds, cases=["Pristine"], correction="none", y_scale="linear", use_global_ylims=False, add_legend=False, dpi=80, figsize=3.0)
    assert set(figs.keys()) >= {"iso", "xx", "yy", "zz"}
