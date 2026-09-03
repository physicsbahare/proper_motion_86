"""Motion-aware target recentering and faint-source forced astrometry.

The base pipeline measures at the catalog coordinate.  For first26 we first preserve
any secure measurement, then search locally for a plausible compact source, and
finally apply an uncertainty-aware forced Gaussian fit when ordinary aperture/
segmentation centroiding is inadequate.  The forced fit is accepted only under the
strict quality cuts in ``forced_astrometry.py``.
"""
from __future__ import annotations

import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources

import pm86.measurement as m
import pm86.pipeline as pipeline
from forced_astrometry import fit_forced_gaussian

# Association radius, not a generic motion search radius.  A much larger radius
# can latch onto an unrelated stationary neighbour in every epoch and produce a
# formally precise but scientifically invalid PM.  Four LW pixels (~0.25 arcsec)
# is deliberately conservative for autonomous rescue; larger motions require a
# trajectory/linking analysis rather than independent nearest-source recentering.
SEARCH_RADIUS_PIX = 4.0
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
        score = (d, -float(vals.sum()))
        if best is None or score < best[0]:
            best = (score, x, y, d)
    if best is None:
        return tx, ty, 0.0, False
    return best[1], best[2], best[3], True


def measure_exposure_recentered(exposure):
    row, controls = m.measure_exposure(exposure)
    base_snr = float(row.get("clean_snr_err", np.nan))
    row["aperture_clean_snr_err"] = base_snr
    row["astrometric_snr"] = base_snr
    row["forced_fit_attempted"] = False
    row["forced_fit_accepted"] = False
    row["forced_snr"] = np.nan

    # Keep an already secure forced-position measurement unchanged.
    if np.isfinite(base_snr) and base_snr >= 3 and np.isfinite(row.get("x_2dg", np.nan)):
        row["target_seed_recentered"] = False
        row["target_seed_offset_pix"] = 0.0
        row["astrometry_source"] = "catalog_position_2dg"
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

    # The local detection is only a seed.  The forced fit provides its own
    # uncertainty/significance test and can also operate directly at the catalog
    # prediction if segmentation finds nothing.
    seed_x, seed_y = (sx, sy) if found else (tx, ty)
    forced = fit_forced_gaussian(data_sub, err, astrom_bad, seed_x, seed_y)
    row.update(forced)

    if forced.get("forced_fit_accepted", False):
        fx = float(forced["forced_x"])
        fy = float(forced["forced_y"])
        total_offset = float(np.hypot(fx - tx, fy - ty))
        row["target_seed_offset_pix"] = total_offset
        # Never accept a fitted solution outside the conservative association
        # envelope.  Larger offsets need explicit cross-epoch trajectory linking.
        if total_offset <= SEARCH_RADIUS_PIX:
            raw_flux, raw_ferr, raw_snr, _ = m.aperture_measure(sci, err, phot_bad, fx, fy)
            clean_flux, clean_ferr, clean_snr, _ = m.aperture_measure(sci, err, astrom_bad, fx, fy)
            cen = m.centroid_stamp(
                data_sub, err, astrom_bad, fx, fy,
                n_mc=0,
                seed=m.stable_seed(f"forced-com|{exposure.candidate_id}|{exposure.filename}"),
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

            e2, n2 = sky(fx, fy)
            ec, nc = sky(cen["x_com"], cen["y_com"])
            forced_snr = float(forced["forced_snr"])
            row.update({
                "raw_flux": raw_flux, "raw_fluxerr": raw_ferr, "raw_snr_err": raw_snr,
                "clean_flux": clean_flux, "clean_fluxerr": clean_ferr,
                "aperture_clean_snr_err": clean_snr,
                # first26's existing registration gate consumes clean_snr_err.
                # For an accepted forced fit, expose its stricter model-fit S/N
                # there while preserving the literal aperture value separately.
                "clean_snr_err": forced_snr,
                "astrometric_snr": forced_snr,
                "astrometry_source": "forced_gaussian",
                "x_2dg": fx, "y_2dg": fy,
                "x_com": cen["x_com"], "y_com": cen["y_com"],
                "sigma_x_pix": float(forced["forced_x_err_pix"]),
                "sigma_y_pix": float(forced["forced_y_err_pix"]),
                "n_mc_good": 0,
                "east_2dg_arcsec_raw_wcs": e2, "north_2dg_arcsec_raw_wcs": n2,
                "east_com_arcsec_raw_wcs": ec, "north_com_arcsec_raw_wcs": nc,
            })
            return row, controls
        row["forced_fit_accepted"] = False

    # If forced fitting fails but segmentation found a real local source, retain the
    # original recentering fallback.  It must still satisfy the ordinary aperture
    # S/N >= 3 gate downstream.
    if not found:
        row["astrometry_source"] = "unusable_catalog_position"
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
        "clean_flux": clean_flux, "clean_fluxerr": clean_ferr,
        "aperture_clean_snr_err": clean_snr, "clean_snr_err": clean_snr,
        "astrometric_snr": clean_snr,
        "astrometry_source": "segmentation_recenter_2dg",
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
