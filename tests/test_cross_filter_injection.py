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


def test_cross_filter_fluxes_and_templates_must_be_independent():
    """A red target cannot be simulated with the same flux in both bands."""
    a = fixture()
    psf_444 = np.array([[0., 1., 0.], [1., 8., 1.], [0., 1., 0.]])
    psf_356 = np.ones((5, 5))
    early, m1 = inject_at_sky(a, 20, 20, psf_444, 120.)
    late, m2 = inject_at_sky(a, 20, 20, psf_356, 8.)
    assert early.sci.sum() == pytest.approx(120.)
    assert late.sci.sum() == pytest.approx(8.)
    assert early.sci[20, 20] > late.sci[20, 20]
    assert m1["flux_injected"] == pytest.approx(120.)
    assert m2["flux_injected"] == pytest.approx(8.)
    assert a.sci.sum() == 0


def test_fractional_phase_injection_preserves_flux_and_position_metadata():
    a = fixture()
    psf = np.array([[0., 1., 0.], [1., 10., 1.], [0., 1., 0.]])
    b, meta = inject_at_sky(a, 20.37, 21.18, psf, 15.)
    assert b.sci.sum() == pytest.approx(15.)
    assert meta["subpixel_phase_x"] == pytest.approx(.37)
    assert meta["subpixel_phase_y"] == pytest.approx(.18)
    assert meta["valid_for_mas_astrometry"] is False


def test_noisy_injection_is_seed_reproducible_and_updates_err():
    a = fixture()
    p = np.ones((5, 5))
    b, _ = inject_at_sky(a, 20, 20, p, 40., rng=np.random.default_rng(42),
                         poisson_sigma=0.3)
    c, _ = inject_at_sky(a, 20, 20, p, 40., rng=np.random.default_rng(42),
                         poisson_sigma=0.3)
    assert np.array_equal(b.sci, c.sci)
    assert np.array_equal(b.err, c.err)
    assert b.err[20, 20] == pytest.approx(np.hypot(1., .3))
    assert a.err[20, 20] == 1.


def test_invalid_injection_inputs_fail_closed():
    a = fixture()
    p = np.ones((5, 5))
    with pytest.raises(ValueError, match="positive"):
        inject_at_sky(a, 20, 20, p, 0.)
    with pytest.raises(ValueError, match="seeded RNG"):
        inject_at_sky(a, 20, 20, p, 1., poisson_sigma=1.)
    a.sci[20, 20] = np.nan
    with pytest.raises(ValueError, match="nonfinite"):
        inject_at_sky(a, 20, 20, p, 1.)
