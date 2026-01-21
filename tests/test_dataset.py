from pathlib import Path

from ga2o3.dataset import Ga2O3Dataset


def test_dataset_resolves_band_files_case_insensitive(tmp_path: Path):
    root = tmp_path
    bands = root / "BANDS" / "pristine"
    bands.mkdir(parents=True)

    (bands / "ALPHA.band").write_text("0 0\n", encoding="utf-8")
    (bands / "BETA.band").write_text("0 0\n", encoding="utf-8")

    ds = Ga2O3Dataset(root=root)
    a, b = ds.band_files("Pristine")
    assert a.name.lower().startswith("alpha")
    assert b.name.lower().startswith("beta")
