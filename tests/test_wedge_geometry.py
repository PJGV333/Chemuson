import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from gui.wedge_geometry import compute_wedge_points


def _midpoint(p0, p1):
    return ((p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0)


def test_wedge_points_ring_no_trim():
    tip, base1, base2 = compute_wedge_points((0.0, 0.0), (10.0, 0.0), 4.0)
    assert tip == (0.0, 0.0)
    mid = _midpoint(base1, base2)
    assert mid[0] == pytest.approx(10.0, abs=1e-6)
    assert mid[1] == pytest.approx(0.0, abs=1e-6)


def test_wedge_points_non_ring_trim():
    tip, base1, base2 = compute_wedge_points(
        (0.0, 0.0),
        (10.0, 0.0),
        4.0,
        trim_start=2.0,
        trim_end=1.0,
    )
    assert tip[0] == pytest.approx(2.0, abs=1e-6)
    assert tip[1] == pytest.approx(0.0, abs=1e-6)
    mid = _midpoint(base1, base2)
    assert mid[0] == pytest.approx(9.0, abs=1e-6)
    assert mid[1] == pytest.approx(0.0, abs=1e-6)
