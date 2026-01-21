from pathlib import Path
import numpy as np

from ga2o3.io.band import read_band
from ga2o3.constants import HARTREE_TO_EV


def test_read_band_parses_ticks_and_data():
    p = Path(__file__).parent / "data" / "sample.BAND"
    b = read_band(p)

    assert b.k.shape == (3,)
    assert b.energies_hartree.shape == (3, 2)
    assert b.xaxis_ticks == (0.0, 0.5, 1.0)

    # check conversion
    assert np.isclose(b.energies_ev[0, 0], -0.010 * HARTREE_TO_EV)
