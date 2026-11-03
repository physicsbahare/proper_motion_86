import numpy as np
import pytest
from dataclasses import dataclass
from pm86.cross_filter_injection import normalized_psf, inject_at_sky, empirical_tail

class LinearWCS:
    def invert(self, ra, dec):
        return ra, dec

@dataclass
class Exposure:
    sci: object
    err: object
    dq: object
    gwcs: object
    x0: int
    y0: int
    ra_deg: float
    dec_deg: float
    x_full: float
    y_full: float

def fixture():
    return Exposure(np.zeros((40, 40)), np.ones((40, 40)),
                    np.zeros((40, 40), dtype=np.uint32),
                    LinearWCS(), 0, 0, 20, 20, 20, 20)

def test_injection_is_copy_and_flux_conserving():
    a = fixture()
    b, meta = inject_at_sky(a, 20, 21, np.ones((5, 5)), 50)
    assert a.sci.sum() == 0
    assert b.sci.sum() == pytest.approx(50)
    assert b.sci[21, 20] == pytest.approx(2)
    assert meta["valid_for_mas_astrometry"] is False

def test_dq_rejection_and_edge():
    a = fixture()
    a.dq[20, 20] = 1
    with pytest.raises(ValueError, match="DQ"):
        inject_at_sky(a, 20, 20, np.ones((5, 5)), 50)
    with pytest.raises(ValueError, match="outside"):
        inject_at_sky(a, 1, 1, np.ones((5, 5)), 50)

def test_psf_rejects_invalid():
    with pytest.raises(ValueError):
        normalized_psf(np.array([[1., -1.], [2., 3.]]))

def test_tail_resolution_not_gaussian_sigma():
    v = np.arange(500)
    result = empirical_tail(v, 999)
    assert result["n_exceed"] == 0
    assert result["p_empirical"] == pytest.approx(1/501)
    assert result["minimum_resolvable_p"] == pytest.approx(1/501)
