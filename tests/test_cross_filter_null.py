"""Fast synthetic tests for the Dave cross-filter stationary-control null.

No MAST data, scientific results, or network are required by these tests.
"""
import numpy as np
import pandas as pd

from pm86.cross_filter_null import (
    conservative_sensitivity, leave_one_out_null, summarize_null,
)


def _controls(xy, flux=100.0):
    xy = np.asarray(xy, dtype=float)
    return pd.DataFrame({
        "east_arcsec_raw_wcs": xy[:, 0],
        "north_arcsec_raw_wcs": xy[:, 1],
        "snr": np.full(len(xy), 20.0),
        "flux": np.linspace(flux, flux * 2.0, len(xy)),
    })


def test_leave_one_out_catches_a_single_moving_reference():
    rng = np.random.default_rng(31)
    coords = rng.uniform(-8, 8, size=(24, 2))
    # Keep all controls mutually separated: no ambiguous deblended neighbors.
    for i in range(1, len(coords)):
        while np.min(np.linalg.norm(coords[i] - coords[:i], axis=1)) < 0.8:
            coords[i] = rng.uniform(-8, 8, size=2)
    early = _controls(coords)
    late_coords = coords + np.array([0.025, -0.013])
    late_coords[5] += np.array([0.070, 0.0])  # 70 mas moving contaminant
    late = _controls(late_coords, flux=60.0)
    null, quality = leave_one_out_null(early, late)
    assert quality["status"] == "LOO_MEASURED"
    assert len(null) >= 20
    assert null.radius_mas.max() >= 45.0
    summary = summarize_null(null)
    assert summary["status"] == "PRELIMINARY_FIELD_CONTROL_NULL"
    assert summary["n"] >= 15
    assert not summary["stationarity_individually_verified"]


def test_insufficient_controls_is_explicit():
    xy = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [2, 0], [0, 2]])
    null, quality = leave_one_out_null(_controls(xy), _controls(xy))
    assert null.empty
    assert quality["status"] == "INSUFFICIENT_ISOLATED_CONTROLS"


def test_conservative_sensitivity_is_not_adopted_significance():
    null = {
        "status": "PRELIMINARY_FIELD_CONTROL_NULL",
        "bias_east_mas": 1.0, "bias_north_mas": -2.0,
        "covariance_mas2": [[100, 0], [0, 100]],
    }
    x = conservative_sensitivity(-55.2, 16.4, 11.5, 10, 497.64, null)
    assert x["status"] == "CONSERVATIVE_SENSITIVITY_ONLY"
    assert x["conservative_2d_mahalanobis_score"] < x["nominal_2d_quadrature_score"]
    assert x["not_a_gaussian_sigma"]
    assert x["registration_noise_may_be_double_counted"]
