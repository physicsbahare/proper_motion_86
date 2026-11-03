"""Fast synthetic tests for the Dave cross-filter stationary-control null.

No MAST data, scientific results, or network are required by these tests.
"""
import numpy as np
import pandas as pd

from pm86.cross_filter_null import (
    annotate_robust_retention, conservative_sensitivity, leave_one_out_null, summarize_null,
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


def test_color_correlation_uses_retained_controls_and_reports_outlier_leverage():
    rng = np.random.default_rng(42)
    n = 20
    color = np.linspace(-1, 1, n)
    east = rng.normal(0, 3, n)
    north = rng.normal(0, 3, n)
    # One control creates a spurious strong color correlation.
    color[-1] = 4
    east[-1] = 100
    north[-1] = 100
    residuals = pd.DataFrame({
        "east_mas": east,
        "north_mas": north,
        "instrumental_color_proxy": color,
    })
    summary = summarize_null(residuals)
    assert summary["status"] == "PRELIMINARY_FIELD_CONTROL_NULL"
    assert summary["n_clipped"] >= 1
    assert summary["color_corr_east_all_pearson_r"] > 0.8
    assert abs(summary["color_corr_east_retained_pearson_r"]) < 0.7
    assert summary["color_corr_east_all_n"] > summary["color_corr_east_retained_n"]
    assert len(summary["color_corr_east_all_leave_one_out_pearson_range"]) == 2


def test_reference_audit_coordinates_and_retention_match_covariance():
    rng = np.random.default_rng(52)
    coords = rng.uniform(-8, 8, size=(24, 2))
    for i in range(1, len(coords)):
        while np.min(np.linalg.norm(coords[i] - coords[:i], axis=1)) < 0.8:
            coords[i] = rng.uniform(-8, 8, size=2)
    early = _controls(coords)
    late = _controls(coords + np.array([0.02, -0.015]), flux=75.0)
    null, qc = leave_one_out_null(early, late)
    assert qc["status"] == "LOO_MEASURED"
    for col in ("early_east_arcsec", "early_north_arcsec",
                "late_east_arcsec", "late_north_arcsec",
                "early_flux_raw", "late_flux_raw"):
        assert col in null.columns
        assert np.isfinite(null[col]).all()
    audit = annotate_robust_retention(null)
    summary = summarize_null(null)
    assert audit["retained_for_covariance"].sum() == summary["n"]
    assert len(audit) == summary["n"] + summary["n_clipped"]
    assert not audit.empty
