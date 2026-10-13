"""Tests for Eq. 9 of Li & Fynbo (2026), without Gaia-to-JWST threshold transfer."""
import math

import pytest
from pm86.mahalanobis import mahalanobis_zero_pm


def test_no_motion_is_zero():
    out = mahalanobis_zero_pm(0, 0, 3, 5, 0.2)
    assert out["mahalanobis_d2"] == 0
    assert out["mahalanobis_p_two_sided"] == 1
    assert out["mahalanobis_gaussian_sigma"] == 0
    assert not out["mahalanobis_calibrated"]


def test_covariance_changes_significance():
    # x=1, y=2, rho=0.5 => D^2=(1+4-2)/(1-0.25) = 4.
    corr = mahalanobis_zero_pm(2, 8, 2, 4, 0.5)
    diag = mahalanobis_zero_pm(2, 8, 2, 4, 0.0)
    assert math.isclose(corr["mahalanobis_d2"], 4, abs_tol=1e-12)
    assert math.isclose(diag["mahalanobis_d2"], 5, abs_tol=1e-12)
    assert math.isclose(corr["mahalanobis_p_two_sided"], math.exp(-2))


@pytest.mark.parametrize("values", [
    (1, 1, 0, 2, 0), (1, 1, 2, 0, 0),
    (1, 1, 2, 2, 1), (1, 1, 2, 2, -1),
    (1, 1, 2, 2, float("nan")),
])
def test_invalid_covariance_rejected(values):
    with pytest.raises(ValueError):
        mahalanobis_zero_pm(*values)
