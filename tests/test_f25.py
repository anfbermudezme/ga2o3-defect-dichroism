from pathlib import Path
import numpy as np

from ga2o3.io.f25 import read_f25_maps


def test_read_f25_maps_reads_two_maps():
    p = Path(__file__).parent / "data" / "sample.f25"
    maps = read_f25_maps(p, max_maps=2)

    assert len(maps) == 2
    assert maps[0].shape == (2, 2)
    assert maps[1].shape == (2, 2)

    assert np.isclose(maps[0][0, 0], 1.0)
    assert np.isclose(maps[1][0, 0], -1.0)
