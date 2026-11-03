"""Empirical cross-filter astrometric null from local field controls.

This module does not assume that all field sources are stationary galaxies.
It measures leave-one-out residuals after *independent* local affine fits.
The resulting scatter is a diagnostic upper-bound sensitivity budget, NOT a
validated target-color PSF systematic or a Gaussian significance calibration.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from scipy.spatial import cKDTree

from .astrometry import _fit_affine, _match_controls
from .config import CONFIG

XY = ["east_arcsec_raw_wcs", "north_arcsec_raw_wcs"]


def _isolated(controls: pd.DataFrame, minimum_separation_arcsec: float = 0.4) -> pd.DataFrame:
    """Avoid close/blended controls without claiming a morphology classification."""
    if controls.empty or not set(XY + ["snr"]).issubset(controls.columns):
        return controls.iloc[0:0].copy()
    c = controls.copy()
    for name in XY + ["snr"]:
        c[name] = pd.to_numeric(c[name], errors="coerce")
    c = c[np.isfinite(c[XY + ["snr"]]).all(axis=1) & (c["snr"] >= 10)].copy()
    if len(c) < 2:
        return c.reset_index(drop=True)
    xy = c[XY].to_numpy(float)
    d, _ = cKDTree(xy).query(xy, k=2)
    return c[d[:, 1] >= minimum_separation_arcsec].reset_index(drop=True)


def leave_one_out_null(early_controls: pd.DataFrame, late_controls: pd.DataFrame,
                       min_fitting_controls: int = 6) -> tuple[pd.DataFrame, dict]:
    """Transform each late control to the early frame, excluding it from the fit.

    The same stars/galaxies are not assumed stationary merely because matched.
    Source proper motions, residual PSF mismatch, centroid noise and registration
    uncertainties may all contribute to the empirical residual distribution.
    """
    early = _isolated(early_controls)
    late = _isolated(late_controls)
    status = {"n_early_controls": len(early), "n_late_controls": len(late)}
    if len(early) < min_fitting_controls + 1 or len(late) < min_fitting_controls + 1:
        return pd.DataFrame(), {**status, "status": "INSUFFICIENT_ISOLATED_CONTROLS"}
    try:
        li, ei = _match_controls(late, early)
    except RuntimeError as exc:
        return pd.DataFrame(), {**status, "status": "MATCH_FAILURE", "reason": str(exc)}
    src = late.iloc[li][XY].to_numpy(float)
    dst = early.iloc[ei][XY].to_numpy(float)
    if len(src) < min_fitting_controls + 1:
        return pd.DataFrame(), {**status, "status": "INSUFFICIENT_MATCHES", "n_matches": len(src)}
    records = []
    for i in range(len(src)):
        mask = np.arange(len(src)) != i
        try:
            coef, good, _, _, _ = _fit_affine(src[mask], dst[mask],
                                              min_matches=min_fitting_controls)
            if good.sum() < min_fitting_controls:
                continue
            transformed = np.r_[1.0, src[i]] @ coef
            residual_mas = (transformed - dst[i]) * 1000.0
            ef = float(early.iloc[ei[i]].get("flux", np.nan))
            lf = float(late.iloc[li[i]].get("flux", np.nan))
            color = -2.5 * np.log10(lf / ef) if ef > 0 and lf > 0 else np.nan
            records.append({
                "early_index": int(ei[i]), "late_index": int(li[i]),
                "east_mas": float(residual_mas[0]),
                "north_mas": float(residual_mas[1]),
                "radius_mas": float(np.hypot(*residual_mas)),
                "instrumental_color_proxy": float(color),
                "snr_early": float(early.iloc[ei[i]]["snr"]),
                "snr_late": float(late.iloc[li[i]]["snr"]),
                "n_registration_controls": int(good.sum()),
            })
        except (RuntimeError, np.linalg.LinAlgError, ValueError):
            continue
    out = pd.DataFrame(records)
    return out, {**status, "status": "LOO_MEASURED" if len(out) else "LOO_FAILED",
                 "n_cross_filter_matches": int(len(src)), "n_loo_success": len(out)}


def summarize_null(residuals: pd.DataFrame) -> dict:
    """Robust central bias and covariance, with an explicit data-sufficiency flag."""
    if residuals.empty or len(residuals) < 7:
        return {"status": "INSUFFICIENT_RESIDUALS", "n": len(residuals)}
    finite_residuals = residuals.loc[np.isfinite(residuals[["east_mas", "north_mas"]]).all(axis=1)].copy()
    values = finite_residuals[["east_mas", "north_mas"]].to_numpy(float)
    if len(values) < 7:
        return {"status": "INSUFFICIENT_FINITE_RESIDUALS", "n": len(values)}
    median = np.median(values, axis=0)
    radial = np.linalg.norm(values - median, axis=1)
    med_r = np.median(radial)
    scale_r = 1.4826 * np.median(np.abs(radial - med_r))
    keep = radial <= med_r + 4.0 * max(scale_r, 1.0)
    clean = values[keep]
    retained_residuals = finite_residuals.iloc[np.flatnonzero(keep)].copy()
    if len(clean) < 7:
        return {"status": "INSUFFICIENT_AFTER_CLIPPING", "n": len(clean)}
    bias = np.median(clean, axis=0)
    cov = np.cov(clean, rowvar=False, ddof=1)
    if cov.shape != (2, 2) or not np.isfinite(cov).all():
        return {"status": "INVALID_COVARIANCE"}
    cdiag = np.sqrt(np.maximum(np.diag(cov), 0))
    result = {
        "status": "PRELIMINARY_FIELD_CONTROL_NULL",
        "n": int(len(clean)),
        "n_clipped": int(len(values) - len(clean)),
        "bias_east_mas": float(bias[0]),
        "bias_north_mas": float(bias[1]),
        "scatter_east_mas": float(cdiag[0]),
        "scatter_north_mas": float(cdiag[1]),
        "covariance_mas2": cov.tolist(),
        "min_empirical_p_resolution": float(1.0 / (len(clean) + 1)),
        "sufficient_for_diagnostic": bool(len(clean) >= 15),
        "stationarity_individually_verified": False,
        "control_morphology_verified": False,
    }
    # Color tests must use the SAME robustly retained controls as the
    # reported covariance. Also disclose pre-clipping influence explicitly.
    for axis, key in [("east_mas", "east"), ("north_mas", "north")]:
        for label, sample in [("retained", retained_residuals), ("all", finite_residuals)]:
            sel = sample[["instrumental_color_proxy", axis]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(sel) >= 8 and sel.instrumental_color_proxy.std() > 0 and sel[axis].std() > 0:
                r, p = pearsonr(sel.instrumental_color_proxy, sel[axis])
                rho, rank_p = spearmanr(sel.instrumental_color_proxy, sel[axis])
                result[f"color_corr_{key}_{label}_n"] = int(len(sel))
                result[f"color_corr_{key}_{label}_pearson_r"] = float(r)
                result[f"color_corr_{key}_{label}_pearson_p"] = float(p)
                result[f"color_corr_{key}_{label}_spearman_rho"] = float(rho)
                result[f"color_corr_{key}_{label}_spearman_p"] = float(rank_p)
                # A single extreme control can drive a small-N Pearson test.
                if len(sel) >= 9:
                    x = sel.instrumental_color_proxy.to_numpy(float)
                    y = sel[axis].to_numpy(float)
                    loo_r = [float(pearsonr(np.delete(x, i), np.delete(y, i)).statistic)
                             for i in range(len(x))
                             if np.std(np.delete(x, i)) > 0 and np.std(np.delete(y, i)) > 0]
                    if loo_r:
                        result[f"color_corr_{key}_{label}_leave_one_out_pearson_range"] = [
                            float(min(loo_r)), float(max(loo_r))
                        ]
    valid_color = retained_residuals["instrumental_color_proxy"].dropna()
    if len(valid_color):
        result["control_instrumental_color_min"] = float(valid_color.min())
        result["control_instrumental_color_max"] = float(valid_color.max())
    return result


def conservative_sensitivity(mu_e: float, mu_n: float, err_e: float, err_n: float,
                             baseline_days: float, null_summary: dict) -> dict:
    """Conservative scenario, NOT an adopted error budget or calibrated sigma.

    Adds total control null covariance to already reported PM errors. As the
    original errors contain registration scatter, this may double count
    registration/centroid noise. Use only to flag robustness for follow-up.
    """
    vals = np.array([mu_e, mu_n, err_e, err_n, baseline_days], dtype=float)
    if not np.isfinite(vals).all() or np.min(vals[2:]) <= 0:
        return {"status": "INVALID_INPUT"}
    if null_summary.get("status") != "PRELIMINARY_FIELD_CONTROL_NULL":
        return {"status": "NO_USABLE_NULL"}
    dt = baseline_days / 365.25
    measured_mas = vals[:2] * dt
    bias = np.array([null_summary["bias_east_mas"], null_summary["bias_north_mas"]])
    stat_cov = np.diag((vals[2:4] * dt)**2)
    sys_cov = np.asarray(null_summary["covariance_mas2"], dtype=float)
    total = stat_cov + sys_cov
    if np.min(np.linalg.eigvalsh(total)) <= 0:
        return {"status": "SINGULAR_TOTAL_COVARIANCE"}
    corrected = measured_mas - bias
    score = float(np.sqrt(max(0.0, corrected @ np.linalg.solve(total, corrected))))
    return {
        "status": "CONSERVATIVE_SENSITIVITY_ONLY",
        "corrected_mu_east_masyr": float(corrected[0] / dt),
        "corrected_mu_north_masyr": float(corrected[1] / dt),
        "nominal_2d_quadrature_score": float(np.hypot(mu_e / err_e, mu_n / err_n)),
        "conservative_2d_mahalanobis_score": score,
        "total_displacement_covariance_mas2": total.tolist(),
        "not_a_gaussian_sigma": True,
        "registration_noise_may_be_double_counted": True,
        "target_color_psf_bias_not_validated": True,
    }
