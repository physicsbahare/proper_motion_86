"""Audit frozen versus independently refitted two-epoch proper motions.

The comparison is diagnostic only: two reductions of the same photons are
correlated, so a quadrature difference significance is NOT valid.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def reconcile_frozen_refit(frozen: dict, refit: dict,
                           exposure_measurements: pd.DataFrame,
                           registered: pd.DataFrame | None = None) -> dict:
    axes = (
        ("east", "mu_alpha_cosdec_masyr", "sigma_mu_alpha_cosdec_masyr"),
        ("north", "mu_delta_masyr", "sigma_mu_delta_masyr"),
    )
    result = {
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "frozen_refit_errors_correlated": True,
        "quadrature_difference_significance_invalid": True,
        "refit_classification": refit.get("classification", "UNKNOWN"),
        "refit_centroid_method_disagreement": refit.get("centroid_method_disagreement"),
        "refit_cross_filter_color_flag": refit.get("cross_filter_color_systematic_flag"),
        "refit_covariance_status": refit.get("pm_covariance_status"),
    }
    if refit.get("pm_status") != "PM_MEASURED":
        return {**result, "status": "REFIT_NOT_MEASURED"}
    for axis, value, error in axes:
        f, r = float(frozen[value]), float(refit[value])
        fe, re = float(frozen[error]), float(refit[error])
        result[f"frozen_mu_{axis}_masyr"] = f
        result[f"refit_mu_{axis}_masyr"] = r
        result[f"delta_mu_{axis}_masyr"] = r - f
        result[f"frozen_sigma_mu_{axis}_masyr"] = fe
        result[f"refit_sigma_mu_{axis}_masyr"] = re
        result[f"error_ratio_{axis}_refit_to_frozen"] = re / fe if fe > 0 else None
    result["frozen_baseline_days"] = float(frozen["baseline_days"])
    result["refit_baseline_days"] = float(refit["baseline_days"])
    result["baseline_difference_days"] = result["refit_baseline_days"] - result["frozen_baseline_days"]
    if exposure_measurements is not None and len(exposure_measurements):
        m = exposure_measurements.copy()
        m["astrometric_snr"] = pd.to_numeric(m.get("astrometric_snr", m.get("clean_snr_err")), errors="coerce")
        m["centroid_finite"] = np.isfinite(pd.to_numeric(m["east_2dg_arcsec_raw_wcs"], errors="coerce")) & np.isfinite(pd.to_numeric(m["north_2dg_arcsec_raw_wcs"], errors="coerce"))
        m["target_snr_ge3"] = m.astrometric_snr >= 3
        result["all_exposure_diagnostics"] = []
        for _, row in m.iterrows():
            result["all_exposure_diagnostics"].append({
                "filename": str(row["filename"]),
                "epoch": str(row["epoch_label"]),
                "filter": str(row["filter"]),
                "mjd": float(row["mjd"]),
                "snr": float(row["astrometric_snr"]),
                "target_snr_ge3": bool(row["target_snr_ge3"]),
                "centroid_finite": bool(row["centroid_finite"]),
                "astrometry_source": str(row.get("astrometry_source", "unknown")),
                "seed_offset_pix": float(row.get("target_seed_offset_pix", np.nan)),
                "forced_fit_accepted": bool(row.get("forced_fit_accepted", False)),
            })
        result["usable_by_epoch"] = {
            epoch: int(sum(row["epoch"] == epoch and row["target_snr_ge3"] and row["centroid_finite"]
                           for row in result["all_exposure_diagnostics"]))
            for epoch in ("early", "late")
        }
    if registered is not None and len(registered):
        result["registered_exposure_filenames"] = registered["filename"].astype(str).tolist()
        result["registered_by_epoch"] = {
            epoch: int((registered["epoch_label"] == epoch).sum())
            for epoch in ("early", "late")
        }
    result["unresolved"] = [
        "Identify original frozen exposure IDs and centroid method for a like-for-like comparison.",
        "Check forced-fit and segmentation recentering associations visually in each epoch.",
        "Compare registered positions and PSF residuals before choosing an estimator.",
        "Do not infer a significance for the frozen/refit difference without their covariance.",
    ]
    return result
