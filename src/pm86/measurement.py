"""Per-exposure photometry, DQ diagnostics, centroids, and field controls."""

from __future__ import annotations

import zlib

import numpy as np
import pandas as pd
from scipy import ndimage
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clipped_stats
from photutils.centroids import centroid_2dg, centroid_com
from photutils.segmentation import detect_sources

from .archive import ExposureCutout
from .config import CONFIG


def robust_sigma(values) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return np.nan
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return float(1.4826 * mad)


def stable_seed(text: str) -> int:
    return int(zlib.crc32(text.encode("utf-8")) & 0xFFFFFFFF)


def masks(sci, err, dq):
    base = ~np.isfinite(sci) | ~np.isfinite(err) | (err <= 0)
    dq64 = np.asarray(dq, dtype=np.uint64)
    phot_bad = base | ((dq64 & np.uint64(CONFIG.phot_hard_bad_bits)) != 0)
    astrom_bad = base | ((dq64 & np.uint64(CONFIG.astrom_bad_bits)) != 0)
    return phot_bad, astrom_bad


def aperture_measure(sci, err, bad, x, y):
    rin = CONFIG.target_annulus_inner_pix
    rout = CONFIG.target_annulus_outer_pix
    rap = CONFIG.target_aperture_radius_pix
    ny, nx = sci.shape
    xmin, xmax = max(0, int(np.floor(x - rout - 1))), min(nx, int(np.ceil(x + rout + 2)))
    ymin, ymax = max(0, int(np.floor(y - rout - 1))), min(ny, int(np.ceil(y + rout + 2)))
    if xmin >= xmax or ymin >= ymax:
        return np.nan, np.nan, np.nan, np.nan

    ss, ee, bb = sci[ymin:ymax, xmin:xmax], err[ymin:ymax, xmin:xmax], bad[ymin:ymax, xmin:xmax]
    yy, xx = np.indices(ss.shape, dtype=float)
    xx += xmin
    yy += ymin
    rr = np.hypot(xx - x, yy - y)
    ap = (rr <= rap) & (~bb) & np.isfinite(ss) & np.isfinite(ee)
    ann = (rr >= rin) & (rr <= rout) & (~bb) & np.isfinite(ss)
    if ap.sum() < 5 or ann.sum() < 15:
        return np.nan, np.nan, np.nan, np.nan
    bkg = float(np.nanmedian(ss[ann]))
    flux = float(np.nansum(ss[ap] - bkg))
    ferr = float(np.sqrt(np.nansum(ee[ap] ** 2)))
    return flux, ferr, flux / ferr if ferr > 0 else np.nan, bkg


def dq_diagnostics(dq, x, y, radius=4.0):
    ny, nx = dq.shape
    xmin, xmax = max(0, int(x - radius - 1)), min(nx, int(x + radius + 2))
    ymin, ymax = max(0, int(y - radius - 1)), min(ny, int(y + radius + 2))
    sub = np.asarray(dq[ymin:ymax, xmin:xmax], dtype=np.uint64)
    yy, xx = np.indices(sub.shape, dtype=float)
    xx += xmin
    yy += ymin
    vals = sub[np.hypot(xx - x, yy - y) <= radius]
    dq_or = int(np.bitwise_or.reduce(vals)) if vals.size else 0
    xc = int(np.clip(round(x), 0, nx - 1))
    yc = int(np.clip(round(y), 0, ny - 1))
    center = int(np.asarray(dq[yc, xc], dtype=np.uint64))
    return {
        "dq_or_r4": dq_or,
        "dq_center": center,
        "dq_do_not_use": bool(dq_or & 1),
        "dq_saturated": bool(dq_or & 2),
        "dq_jump_det": bool(dq_or & 4),
        "dq_outlier": bool(dq_or & 16),
    }


def centroid_stamp(data_sub, err, bad, x_expected, y_expected, n_mc=0, seed=1):
    out = {
        "x_2dg": np.nan, "y_2dg": np.nan, "x_com": np.nan, "y_com": np.nan,
        "sigma_x_pix": np.nan, "sigma_y_pix": np.nan, "n_mc_good": 0,
    }
    h = CONFIG.centroid_half_size
    ny, nx = data_sub.shape
    xc, yc = int(round(x_expected)), int(round(y_expected))
    xa, xb = max(0, xc - h), min(nx, xc + h + 1)
    ya, yb = max(0, yc - h), min(ny, yc + h + 1)
    if xb - xa < 7 or yb - ya < 7:
        return out

    stamp = np.asarray(data_sub[ya:yb, xa:xb], dtype=float)
    estamp = np.asarray(err[ya:yb, xa:xb], dtype=float)
    m = np.asarray(bad[ya:yb, xa:xb], dtype=bool)
    m |= ~np.isfinite(stamp) | ~np.isfinite(estamp) | (estamp <= 0)

    try:
        x, y = centroid_2dg(stamp, error=estamp, mask=m)
        x, y = xa + float(x), ya + float(y)
        if abs(x - x_expected) <= h - 1 and abs(y - y_expected) <= h - 1:
            out["x_2dg"], out["y_2dg"] = x, y
    except Exception:
        pass

    try:
        p = np.array(stamp, copy=True)
        p[m] = 0.0
        p = np.clip(p, 0.0, None)
        if p.sum() > 0:
            x, y = centroid_com(p)
            out["x_com"], out["y_com"] = xa + float(x), ya + float(y)
    except Exception:
        pass

    if n_mc and np.isfinite(out["x_2dg"]) and np.isfinite(out["y_2dg"]):
        rng = np.random.default_rng(seed)
        xs, ys = [], []
        for _ in range(int(n_mc)):
            try:
                noisy = stamp + rng.normal(0.0, estamp)
                x, y = centroid_2dg(noisy, error=estamp, mask=m)
                x, y = xa + float(x), ya + float(y)
                if abs(x - out["x_2dg"]) < 3 and abs(y - out["y_2dg"]) < 3:
                    xs.append(x)
                    ys.append(y)
            except Exception:
                pass
        if len(xs) >= max(10, n_mc // 4):
            out["sigma_x_pix"] = float(np.std(xs, ddof=1))
            out["sigma_y_pix"] = float(np.std(ys, ddof=1))
            out["n_mc_good"] = len(xs)
    return out


def local_pixel_to_sky(exposure: ExposureCutout, x_local, y_local):
    result = exposure.gwcs(float(x_local + exposure.x0), float(y_local + exposure.y0))
    return float(np.asarray(result[0]).squeeze()), float(np.asarray(result[1]).squeeze())


def sky_to_tangent(ra, dec, origin_ra, origin_dec):
    origin = SkyCoord(float(origin_ra) * u.deg, float(origin_dec) * u.deg)
    point = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg)
    east, north = origin.spherical_offsets_to(point)
    return np.asarray(east.to_value(u.arcsec)), np.asarray(north.to_value(u.arcsec))


def empirical_blank_noise(sci, err, bad, source_mask, tx, ty, seed):
    rng = np.random.default_rng(seed)
    ny, nx = sci.shape
    margin = int(np.ceil(CONFIG.target_annulus_outer_pix + 2))
    values = []
    tries = 0
    while len(values) < CONFIG.blank_apertures and tries < CONFIG.blank_apertures * 60:
        tries += 1
        x = rng.uniform(margin, nx - 1 - margin)
        y = rng.uniform(margin, ny - 1 - margin)
        if np.hypot(x - tx, y - ty) < 20:
            continue
        xi, yi = int(round(x)), int(round(y))
        r = int(np.ceil(CONFIG.target_aperture_radius_pix + 2))
        if np.any(source_mask[max(0, yi-r):min(ny, yi+r+1), max(0, xi-r):min(nx, xi+r+1)]):
            continue
        flux, _, _, _ = aperture_measure(sci, err, bad, x, y)
        if np.isfinite(flux):
            values.append(flux)
    if len(values) < 20:
        return np.nan, len(values)
    sig = robust_sigma(values)
    if not np.isfinite(sig) or sig <= 0:
        sig = float(np.std(values, ddof=1))
    return float(sig), len(values)


def detect_controls(exposure, data_sub, err, bad, target_x, target_y):
    _, _, std = sigma_clipped_stats(data_sub, mask=bad, sigma=3.0, maxiters=5)
    if not np.isfinite(std) or std <= 0:
        return pd.DataFrame(), np.zeros_like(data_sub, dtype=bool)

    seg = detect_sources(data_sub, CONFIG.control_detect_sigma * std, npixels=CONFIG.control_min_pixels, mask=bad)
    if seg is None:
        return pd.DataFrame(), np.zeros_like(data_sub, dtype=bool)

    source_mask = ndimage.binary_dilation(seg.data > 0, iterations=4)
    rows = []
    ny, nx = data_sub.shape
    for label in np.unique(seg.data):
        if label <= 0:
            continue
        sm = (seg.data == label) & (~bad) & np.isfinite(data_sub) & np.isfinite(err)
        if sm.sum() < CONFIG.control_min_pixels:
            continue
        yy, xx = np.nonzero(sm)
        vals = data_sub[sm]
        pos = np.clip(vals, 0, None)
        if pos.sum() <= 0:
            continue
        x0 = float(np.sum(xx * pos) / pos.sum())
        y0 = float(np.sum(yy * pos) / pos.sum())
        if (
            x0 < CONFIG.control_edge_pix or y0 < CONFIG.control_edge_pix
            or x0 > nx - 1 - CONFIG.control_edge_pix or y0 > ny - 1 - CONFIG.control_edge_pix
            or np.hypot(x0 - target_x, y0 - target_y) < CONFIG.control_target_exclusion_pix
        ):
            continue
        flux = float(np.sum(vals))
        ferr = float(np.sqrt(np.sum(err[sm] ** 2)))
        snr = flux / ferr if ferr > 0 else np.nan
        if not np.isfinite(snr) or snr < CONFIG.control_min_snr:
            continue
        cen = centroid_stamp(data_sub, err, bad, x0, y0, n_mc=0)
        x = cen["x_2dg"] if np.isfinite(cen["x_2dg"]) else x0
        y = cen["y_2dg"] if np.isfinite(cen["y_2dg"]) else y0
        try:
            ra, dec = local_pixel_to_sky(exposure, x, y)
        except Exception:
            continue
        rows.append({"x_local": x, "y_local": y, "ra_deg": ra, "dec_deg": dec, "flux": flux, "fluxerr": ferr, "snr": snr})
    return pd.DataFrame(rows), source_mask


def measure_exposure(exposure: ExposureCutout) -> tuple[dict, pd.DataFrame]:
    sci = np.asarray(exposure.sci, dtype=float)
    err = np.asarray(exposure.err, dtype=float)
    dq = np.asarray(exposure.dq, dtype=np.uint32)
    phot_bad, astrom_bad = masks(sci, err, dq)

    tx = exposure.x_full - exposure.x0
    ty = exposure.y_full - exposure.y0
    _, bkg, _ = sigma_clipped_stats(sci, mask=phot_bad, sigma=3.0, maxiters=5)
    data_sub = sci - float(bkg)

    controls, source_mask = detect_controls(exposure, data_sub, err, astrom_bad, tx, ty)
    raw_flux, raw_ferr, raw_snr, _ = aperture_measure(sci, err, phot_bad, tx, ty)
    clean_flux, clean_ferr, clean_snr, _ = aperture_measure(sci, err, astrom_bad, tx, ty)
    emp_noise, n_blank = empirical_blank_noise(sci, err, phot_bad, source_mask, tx, ty, stable_seed(f"{exposure.candidate_id}|{exposure.filename}"))
    raw_snr_emp = raw_flux / emp_noise if np.isfinite(emp_noise) and emp_noise > 0 else np.nan

    clean_cen = centroid_stamp(
        data_sub, err, astrom_bad, tx, ty, n_mc=CONFIG.centroid_mc_draws,
        seed=stable_seed(f"centroid|{exposure.candidate_id}|{exposure.filename}"),
    )

    def centroid_sky(x, y):
        if not (np.isfinite(x) and np.isfinite(y)):
            return np.nan, np.nan, np.nan, np.nan
        try:
            ra, dec = local_pixel_to_sky(exposure, x, y)
            east, north = sky_to_tangent(ra, dec, exposure.ra_deg, exposure.dec_deg)
            return ra, dec, float(east), float(north)
        except Exception:
            return np.nan, np.nan, np.nan, np.nan

    _, _, e2, n2 = centroid_sky(clean_cen["x_2dg"], clean_cen["y_2dg"])
    _, _, ec, nc = centroid_sky(clean_cen["x_com"], clean_cen["y_com"])

    row = {
        "candidate_id": exposure.candidate_id,
        "epoch_label": exposure.epoch_label,
        "filter": exposure.filter_name,
        "filename": exposure.filename,
        "dataURI": exposure.data_uri,
        "mjd": exposure.mjd,
        "detector": exposure.detector,
        "x_full": exposure.x_full,
        "y_full": exposure.y_full,
        "raw_flux": raw_flux,
        "raw_fluxerr": raw_ferr,
        "raw_snr_err": raw_snr,
        "empirical_ap_sigma": emp_noise,
        "raw_snr_emp": raw_snr_emp,
        "n_blank_apertures": n_blank,
        "clean_flux": clean_flux,
        "clean_fluxerr": clean_ferr,
        "clean_snr_err": clean_snr,
        "x_2dg": clean_cen["x_2dg"], "y_2dg": clean_cen["y_2dg"],
        "x_com": clean_cen["x_com"], "y_com": clean_cen["y_com"],
        "sigma_x_pix": clean_cen["sigma_x_pix"], "sigma_y_pix": clean_cen["sigma_y_pix"],
        "n_mc_good": clean_cen["n_mc_good"],
        "east_2dg_arcsec_raw_wcs": e2, "north_2dg_arcsec_raw_wcs": n2,
        "east_com_arcsec_raw_wcs": ec, "north_com_arcsec_raw_wcs": nc,
        **dq_diagnostics(dq, tx, ty),
    }

    if len(controls):
        controls = controls.copy()
        controls["candidate_id"] = exposure.candidate_id
        controls["epoch_label"] = exposure.epoch_label
        controls["filter"] = exposure.filter_name
        controls["filename"] = exposure.filename
        east, north = sky_to_tangent(controls["ra_deg"].values, controls["dec_deg"].values, exposure.ra_deg, exposure.dec_deg)
        controls["east_arcsec_raw_wcs"] = east
        controls["north_arcsec_raw_wcs"] = north

    return row, controls
