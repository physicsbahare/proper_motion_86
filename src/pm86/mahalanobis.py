"""Covariance-aware significance against the zero proper-motion hypothesis.

Li & Fynbo (2026), arXiv:2610.13742, Eq. 9. Their *Gaia* empirical
scale (1.1485) and thresholds are NOT JWST calibration constants.
No significance is returned when the RA/Dec correlation is unavailable.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.stats import chi2, norm


def mahalanobis_zero_pm(mu_ra_cosdec, mu_dec, sigma_ra, sigma_dec, rho,
                        null_scale=1.0):
    """Return D_M^2 and the Gaussian-equivalent two-sided significance.

    The returned p-value assumes a scaled chi-square with 2 d.o.f.;
    null_scale=1 is an uncalibrated *theoretical* JWST diagnostic, not a
    validated detection threshold. Passing Gaia's 1.1485 is inappropriate
    unless it has been separately calibrated on the selected JWST sample.

    Raises ValueError when the covariance is unavailable, invalid or singular.
    """
    vals = np.array([mu_ra_cosdec, mu_dec, sigma_ra, sigma_dec, rho, null_scale],
                    dtype=float)
    if not np.all(np.isfinite(vals)):
        raise ValueError("Finite proper motions, errors, correlation and scale required")
    if sigma_ra <= 0 or sigma_dec <= 0 or null_scale <= 0:
        raise ValueError("Proper-motion errors and scale must be positive")
    if abs(rho) >= 1:
        raise ValueError("The RA/Dec covariance matrix must be positive definite")
    x = mu_ra_cosdec / sigma_ra
    y = mu_dec / sigma_dec
    d2 = float((x*x + y*y - 2*rho*x*y) / (1-rho*rho))
    d2 = max(0.0, d2)  # floating-point rounding at the null
    p = float(chi2.sf(d2 / null_scale, df=2))
    # Avoid infinity from 1-p rounding at the null.
    sig = float(norm.isf(p / 2)) if p < 1 else 0.0
    return {"mahalanobis_d2": d2, "mahalanobis_p_two_sided": p,
            "mahalanobis_gaussian_sigma": sig, "mahalanobis_null_scale": null_scale,
            "mahalanobis_calibrated": False}


def estimate_local_registration_rho(registered, matched_controls,
                                    sigma_e, sigma_n, pm_result):
    """Estimate RA/Dec correlation from the *same* local field-registration residuals.

    Propagates cross-covariances of per-exposure control residuals into
    inverse-variance weighted early/late epoch means. This is a provisional
    local estimate: the full reference-frame and centroid correlation should
    be verified with stationary-source injection/bootstrap before publication.
    """
    if matched_controls is None or len(matched_controls) == 0:
        return math.nan, "NO_CONTROL_RESIDUALS"
    reg = registered.reset_index(drop=True).copy()
    ctrl = matched_controls.copy()
    required = {"filename", "resid_east_arcsec", "resid_north_arcsec", "used"}
    if not required.issubset(ctrl.columns):
        return math.nan, "MISSING_CONTROL_COVARIANCE_COLUMNS"

    cov_by_file = {}
    for fname, group in ctrl.groupby("filename"):
        good = group[group["used"].astype(bool)]
        x = np.asarray(good["resid_east_arcsec"], dtype=float)
        y = np.asarray(good["resid_north_arcsec"], dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 6 or np.std(x[ok]) <= 0 or np.std(y[ok]) <= 0:
            continue
        corr = float(np.corrcoef(x[ok], y[ok])[0, 1])
        if not np.isfinite(corr):
            continue
        cov_by_file[str(fname)] = np.clip(corr, -0.95, 0.95)

    if not cov_by_file:
        return math.nan, "INSUFFICIENT_CORRELATED_CONTROLS"

    rho_displacement_num = 0.0
    for epoch in ("early", "late"):
        subset = reg[reg["epoch_label"].eq(epoch)]
        idx = subset.index.to_numpy()
        if len(idx) == 0:
            return math.nan, "NO_TWO_EPOCH_ESTIMATE"
        # _weighted_epoch uses weights proportional to inverse variance.
        e_err = np.asarray(sigma_e, dtype=float)[idx]
        n_err = np.asarray(sigma_n, dtype=float)[idx]
        valid_e = np.isfinite(e_err) & (e_err > 0)
        valid_n = np.isfinite(n_err) & (n_err > 0)
        if not (np.all(valid_e) and np.all(valid_n)):
            return math.nan, "INVALID_PER_EXPOSURE_UNCERTAINTY"
        we = (1 / e_err**2)
        wn = (1 / n_err**2)
        we /= we.sum()
        wn /= wn.sum()
        for i, (_, row) in enumerate(subset.iterrows()):
            fname = str(row["filename"])
            # Reference exposure's affine transform is the identity.
            corr = cov_by_file.get(fname, 0.0)
            s_reg_e = float(row["registration_scatter_east_mas"])
            s_reg_n = float(row["registration_scatter_north_mas"])
            if np.isfinite(s_reg_e) and np.isfinite(s_reg_n):
                rho_displacement_num += (we[i] * wn[i] * corr
                                         * max(0.0, s_reg_e) * max(0.0, s_reg_n))

    denom = (float(pm_result["sigma_delta_east_mas"]) *
             float(pm_result["sigma_delta_north_mas"]))
    if not np.isfinite(denom) or denom <= 0:
        return math.nan, "INVALID_DISPLACEMENT_ERROR"
    rho = rho_displacement_num / denom
    if not np.isfinite(rho):
        return math.nan, "INVALID_COVARIANCE"
    # Empirical, provisional correlation. Never allow a singular matrix.
    return float(np.clip(rho, -0.95, 0.95)), "PROVISIONAL_LOCAL_CONTROL_COVARIANCE"
