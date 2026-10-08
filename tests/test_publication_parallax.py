"""Pure-geometry and statistical tests; no archive/network requests."""
import numpy as np
import pytest
from scripts.publication_parallax_audit import (
    tangent_basis, parallax_factors_for_xyz,
    pm_after_assumed_parallax, radial_significance_to_2sided_gaussian,
)

def test_tangent_orthogonal():
    east, north = tangent_basis(150.3, 2.3)
    assert np.linalg.norm(east) == pytest.approx(1)
    assert np.linalg.norm(north) == pytest.approx(1)
    assert np.dot(east, north) == pytest.approx(0, abs=1e-12)

def test_projected_parallax_sign_and_units():
    factors = parallax_factors_for_xyz(0, 0, [0, 1, 0])
    assert factors == pytest.approx([-1, 0], abs=1e-12)
    factors_n = parallax_factors_for_xyz(0, 0, [0, 0, 1])
    assert factors_n == pytest.approx([0, -1], abs=1e-12)

def test_zero_differential_parallax_leaves_pm_unchanged():
    corrected = pm_after_assumed_parallax(-55, 16, 1.36, [0, 0], 11.2)
    assert corrected == pytest.approx([-55, 16])
    no_pi = pm_after_assumed_parallax(-55, 16, 1.36, [1, -1], None)
    assert no_pi == pytest.approx([-55, 16])

def test_nonzero_parallax_correction_is_distance_dependent():
    a = pm_after_assumed_parallax(50, 10, 2, [0.1, -0.2], 10)
    b = pm_after_assumed_parallax(50, 10, 2, [0.1, -0.2], 100)
    assert a == pytest.approx([45, 20])
    assert b == pytest.approx([49.5, 11])

def test_two_dimensional_radial_statistic_not_exact_1d_sigma():
    p, s = radial_significance_to_2sided_gaussian(5)
    assert p == pytest.approx(np.exp(-12.5))
    assert 4.0 < s < 5.0

def test_bad_significance_is_not_silently_converted():
    p, s = radial_significance_to_2sided_gaussian(-1)
    assert np.isnan(p) and np.isnan(s)
