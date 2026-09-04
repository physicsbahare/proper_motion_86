"""Uncertainty-aware forced Gaussian astrometry for faint first26 targets.

This module provides a conservative fallback when aperture S/N or segmentation-based
recentering is insufficient.  It fits a compact 2-D Gaussian plus constant background
in a local stamp, using the CAL ERR array as pixel weights.  A fit is considered
astrometrically usable only if it satisfies explicit significance, width, positional
uncertainty, reduced-chi2, displacement, and compactness bounds.

The intent is to recover real faint point sources, not to turn arbitrary noise peaks
or extended residuals into proper-motion detections.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


MIN_FORCED_SNR = 4.0
MAX_POSITION_SIGMA_PIX = 1.2
MIN_SIGMA_PIX = 0.55
# NIRCam point sources in the red filters used here have Gaussian-equivalent
# sigmas of order ~1 pixel.  A 3-pixel sigma (the previous bound) corresponds
# to FWHM ~7 pixels and can admit extended/noise residuals as astrometric
# targets.  1.8 pixels is deliberately generous while still enforcing a
# point-source-like rescue measurement.
MAX_SIGMA_PIX = 1.8
MAX_AXIS_RATIO = 2.0
MAX_REDUCED_CHI2 = 3.0
FIT_HALF_SIZE = 7
FIT_CENTER_FREEDOM_PIX = 4.0


def _model(params, xx, yy):
    amp, x0, y0, sx, sy, bkg = params
    return bkg + amp * np.exp(
        -0.5 * (((xx - x0) / sx) ** 2 + ((yy - y0) / sy) ** 2)
    )


def fit_forced_gaussian(data_sub, err, bad, x_seed, y_seed):
    """Fit a compact Gaussian around ``(x_seed, y_seed)``.

    Returns a dictionary containing fit diagnostics and ``accepted``.  The fit uses
    only unmasked finite pixels and never assumes the source is detected beforehand.
    """
    out = {
        "forced_fit_attempted": True,
        "forced_fit_accepted": False,
        "forced_x": np.nan,
        "forced_y": np.nan,
        "forced_amp": np.nan,
        "forced_amp_err": np.nan,
        "forced_snr": np.nan,
        "forced_sigma_x_pix": np.nan,
        "forced_sigma_y_pix": np.nan,
        "forced_axis_ratio": np.nan,
        "forced_x_err_pix": np.nan,
        "forced_y_err_pix": np.nan,
        "forced_reduced_chi2": np.nan,
        "forced_npix": 0,
    }

    ny, nx = data_sub.shape
    xc, yc = int(round(x_seed)), int(round(y_seed))
    xa, xb = max(0, xc - FIT_HALF_SIZE), min(nx, xc + FIT_HALF_SIZE + 1)
    ya, yb = max(0, yc - FIT_HALF_SIZE), min(ny, yc + FIT_HALF_SIZE + 1)
    if xb - xa < 9 or yb - ya < 9:
        return out

    stamp = np.asarray(data_sub[ya:yb, xa:xb], dtype=float)
    estamp = np.asarray(err[ya:yb, xa:xb], dtype=float)
    m = np.asarray(bad[ya:yb, xa:xb], dtype=bool)
    good = (~m) & np.isfinite(stamp) & np.isfinite(estamp) & (estamp > 0)
    if good.sum() < 35:
        return out

    yy, xx = np.indices(stamp.shape, dtype=float)
    xx += xa
    yy += ya
    xg, yg = xx[good], yy[good]
    zg, eg = stamp[good], estamp[good]
    out["forced_npix"] = int(good.sum())

    med = float(np.nanmedian(zg))
    positive_peak = float(max(np.nanpercentile(zg, 99) - med, np.nanmax(zg) - med, 0.0))
    if not np.isfinite(positive_peak) or positive_peak <= 0:
        return out

    p0 = np.array([positive_peak, x_seed, y_seed, 1.0, 1.0, med], dtype=float)
    lo = np.array([
        0.0,
        x_seed - FIT_CENTER_FREEDOM_PIX,
        y_seed - FIT_CENTER_FREEDOM_PIX,
        MIN_SIGMA_PIX,
        MIN_SIGMA_PIX,
        med - 10.0 * np.nanstd(zg),
    ])
    hi = np.array([
        max(positive_peak * 20.0, np.nanmax(np.abs(zg)) * 20.0, 1e-12),
        x_seed + FIT_CENTER_FREEDOM_PIX,
        y_seed + FIT_CENTER_FREEDOM_PIX,
        MAX_SIGMA_PIX,
        MAX_SIGMA_PIX,
        med + 10.0 * np.nanstd(zg),
    ])

    def resid(p):
        return (_model(p, xg, yg) - zg) / eg

    try:
        fit = least_squares(resid, p0, bounds=(lo, hi), method="trf", max_nfev=800)
    except Exception:
        return out
    if not fit.success or not np.all(np.isfinite(fit.x)):
        return out

    amp, x0, y0, sx, sy, _ = fit.x
    dof = max(1, len(zg) - len(fit.x))
    chi2 = float(np.sum(fit.fun ** 2))
    rchi2 = chi2 / dof

    perr = np.full(len(fit.x), np.nan)
    try:
        jtj = fit.jac.T @ fit.jac
        cov = np.linalg.pinv(jtj) * rchi2
        diag = np.diag(cov)
        perr = np.sqrt(np.where(diag >= 0, diag, np.nan))
    except Exception:
        pass

    amp_err, xerr, yerr = perr[0], perr[1], perr[2]
    snr = float(amp / amp_err) if np.isfinite(amp_err) and amp_err > 0 else np.nan
    axis_ratio = float(max(sx, sy) / min(sx, sy)) if min(sx, sy) > 0 else np.inf
    accepted = bool(
        np.isfinite(snr)
        and snr >= MIN_FORCED_SNR
        and np.isfinite(xerr) and xerr <= MAX_POSITION_SIGMA_PIX
        and np.isfinite(yerr) and yerr <= MAX_POSITION_SIGMA_PIX
        and MIN_SIGMA_PIX < sx < MAX_SIGMA_PIX
        and MIN_SIGMA_PIX < sy < MAX_SIGMA_PIX
        and np.isfinite(axis_ratio) and axis_ratio <= MAX_AXIS_RATIO
        and np.isfinite(rchi2) and rchi2 <= MAX_REDUCED_CHI2
    )

    out.update({
        "forced_fit_accepted": accepted,
        "forced_x": float(x0),
        "forced_y": float(y0),
        "forced_amp": float(amp),
        "forced_amp_err": float(amp_err) if np.isfinite(amp_err) else np.nan,
        "forced_snr": snr,
        "forced_sigma_x_pix": float(sx),
        "forced_sigma_y_pix": float(sy),
        "forced_axis_ratio": axis_ratio,
        "forced_x_err_pix": float(xerr) if np.isfinite(xerr) else np.nan,
        "forced_y_err_pix": float(yerr) if np.isfinite(yerr) else np.nan,
        "forced_reduced_chi2": float(rchi2),
    })
    return out
