"""Motion-aware target recentering for the first26 audit.

The base pipeline forced target photometry exactly at the catalog coordinate before
accepting a centroid. That can reject the very sources whose proper motion moves
them several NIRCam pixels between the catalog reference position and an archive
exposure. This patch searches a small local region for a compact source, recenters
photometry on the nearest plausible detection, and records the search offset. Field
controls and the downstream local-registration/PM inference remain unchanged.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources

import pm86.measurement as m
import pm86.pipeline as pipeline

SEARCH_RADIUS_PIX = 12.0
SEARCH_SIGMA = 2.5
SEARCH_MIN_PIXELS = 4


def _nearest_seed(data_sub, bad, tx, ty):
    _, _, std = sigma_clipped_stats(data_sub, mask=bad, sigma=3.0, maxiters=5)
    if not np.isfinite(std) or std <= 0:
        return tx, ty, 0.0, False
    seg = detect_sources(data_sub, SEARCH_SIGMA * std, npixels=SEARCH_MIN_PIXELS, mask=bad)
    if seg is None:
        return tx, ty, 0.0, False
    best = None
    for label in np.unique(seg.data):
        if label <= 0:
            continue
        sm = (seg.data == label) & (~bad) & np.isfinite(data_sub)
        if sm.sum() < SEARCH_MIN_PIXELS:
            continue
        yy, xx = np.nonzero(sm)
        vals = np.clip(data_sub[sm], 0.0, None)
        if vals.sum() <= 0:
            continue
        x = float(np.sum(xx * vals) / vals.sum())
        y = float(np.sum(yy * vals) / vals.sum())
        d = float(np.hypot(x - tx, y - ty))
        if d > SEARCH_RADIUS_PIX:
            continue
        # Prefer proximity first, then integrated positive signal.
        score = (d, -float(vals.sum()))
        if best is None or score < best[0]:
            best = (score, x, y, d)
    if best is None:
        return tx, ty, 0.0, False
    return best[1], best[2], best[3], True


def measure_exposure_recentered(exposure):
    row, controls = m.measure_exposure(exposure)
    # Keep an already secure forced-position measurement unchanged.
    if np.isfinite(row.get("clean_snr_err", np.nan)) and row["clean_snr_err"] >= 3:
        row["target_seed_recentered"] = False
        row["target_seed_offset_pix"] = 0.0
        return row, controls

    sci = np.asarray(exposure.sci, dtype=float)
    err = np.asarray(exposure.err, dtype=float)
    dq = np.asarray(exposure.dq, dtype=np.uint32)
    phot_bad, astrom_bad = m.masks(sci, err, dq)
    tx = exposure.x_full - exposure.x0
    ty = exposure.y_full - exposure.y0
    _, bkg, _ = sigma_clipped_stats(sci, mask=phot_bad, sigma=3.0, maxiters=5)
    data_sub = sci - float(bkg)
    sx, sy, offset, found = _nearest_seed(data_sub, astrom_bad, tx, ty)
    row["target_seed_recentered"] = bool(found)
    row["target_seed_offset_pix"] = float(offset)
    if not found:
        return row, controls

    raw_flux, raw_ferr, raw_snr, _ = m.aperture_measure(sci, err, phot_bad, sx, sy)
    clean_flux, clean_ferr, clean_snr, _ = m.aperture_measure(sci, err, astrom_bad, sx, sy)
    cen = m.centroid_stamp(
        data_sub, err, astrom_bad, sx, sy,
        n_mc=m.CONFIG.centroid_mc_draws,
        seed=m.stable_seed(f"recenter|{exposure.candidate_id}|{exposure.filename}"),
    )

    def sky(x, y):
        if not (np.isfinite(x) and np.isfinite(y)):
            return np.nan, np.nan
        try:
            ra, dec = m.local_pixel_to_sky(exposure, x, y)
            east, north = m.sky_to_tangent(ra, dec, exposure.ra_deg, exposure.dec_deg)
            return float(east), float(north)
        except Exception:
            return np.nan, np.nan

    e2, n2 = sky(cen["x_2dg"], cen["y_2dg"])
    ec, nc = sky(cen["x_com"], cen["y_com"])
    row.update({
        "raw_flux": raw_flux, "raw_fluxerr": raw_ferr, "raw_snr_err": raw_snr,
        "clean_flux": clean_flux, "clean_fluxerr": clean_ferr, "clean_snr_err": clean_snr,
        "x_2dg": cen["x_2dg"], "y_2dg": cen["y_2dg"],
        "x_com": cen["x_com"], "y_com": cen["y_com"],
        "sigma_x_pix": cen["sigma_x_pix"], "sigma_y_pix": cen["sigma_y_pix"],
        "n_mc_good": cen["n_mc_good"],
        "east_2dg_arcsec_raw_wcs": e2, "north_2dg_arcsec_raw_wcs": n2,
        "east_com_arcsec_raw_wcs": ec, "north_com_arcsec_raw_wcs": nc,
    })
    return row, controls


def install():
    pipeline.measure_exposure = measure_exposure_recentered
