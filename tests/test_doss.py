from pathlib import Path
import numpy as np

from ga2o3.io.doss import read_doss


def test_read_doss_parses_numeric_and_metadata():
    p = Path(__file__).parent / "data" / "sample_alpha.DOSS"
    d = read_doss(p)

    assert d.energies_hartree.shape == (3,)
    assert d.dos.shape == (3, 3)
    assert d.metadata.get("nepts") == 3
    assert d.metadata.get("nproj") == 3
    assert d.metadata.get("nspin") == 2

    # conversion property
    assert np.isclose(d.energies_ev[1], 0.0)
