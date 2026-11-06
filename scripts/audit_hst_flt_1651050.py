#!/usr/bin/env python3
"""Original HST FLT/FLC pixel detectability check for COSMOS source 1651050.

This is a *screening* test. Even a strong forced signal or Gaussian centroid
is not an independent astrometric epoch until source identification, local
registration, and cross-instrument systematics are independently validated.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.mast import Observations
from scipy.optimize import least_squares
from pm86.archive import mast_download_url

RA, DEC = 150.307303724545, 2.27642841317716


def forced_aperture(sci, err, dq, x, y, radius=2.5):
    """Background-subtracted *forced* S/N, not a blind source detection."""
    sci = np.asarray(sci, dtype=float)
    err = np.asarray(err, dtype=float)
    dq = np.asarray(dq)
    yy, xx = np.indices(sci.shape)
    rr = np.hypot(xx-x, yy-y)
    valid = np.isfinite(sci) & np.isfinite(err) & (err > 0) & (dq == 0)
    ap = valid & (rr <= radius)
    ring = valid & (rr >= 7) & (rr <= 12)
    if ap.sum() < 8 or ring.sum() < 30:
        return {"status": "INSUFFICIENT_VALID_PIXELS"}
    bg = float(np.median(sci[ring]))
    scatter = float(1.4826 * np.median(np.abs(sci[ring] - bg)))
    signal = float(np.sum(sci[ap] - bg))
    pipe_unc = float(np.sqrt(np.sum(err[ap]**2)))
    empirical_unc = float(scatter * np.sqrt(ap.sum() + ap.sum()**2 / ring.sum()))
    uncertainty = max(pipe_unc, empirical_unc)
    return {
        "status": "FORCED_PHOTOMETRY_ONLY",
        "aperture_pixels": int(ap.sum()),
        "background_pixels": int(ring.sum()),
        "net_flux_native": signal,
        "pipeline_unc_native": pipe_unc,
        "empirical_unc_native": empirical_unc,
        "conservative_snr": float(signal / uncertainty) if uncertainty > 0 else np.nan,
        "background_native": bg,
        "bad_aperture_pixels": int(np.count_nonzero((rr <= radius) & ~valid)),
    }


def gaussian_screen(sci, err, dq, x, y):
    """Exploratory Gaussian fit; NEVER an accepted astrometric centroid."""
    yy, xx = np.indices(sci.shape)
    rr = np.hypot(xx-x, yy-y)
    good = (rr <= 6) & (dq == 0) & np.isfinite(sci) & np.isfinite(err) & (err > 0)
    if good.sum() < 30:
        return {"fit_status": "INSUFFICIENT_FIT_PIXELS"}
    xg, yg, z, e = xx[good], yy[good], sci[good], err[good]
    background = float(np.median(z))
    amp = max(float(np.max(z) - background), 0.01)
    def residual(p):
        a, xc, yc, sigma, bg = p
        model = bg + a * np.exp(-0.5*((xg-xc)**2+(yg-yc)**2)/sigma**2)
        return (model-z)/e
    try:
        fit = least_squares(residual, [amp, x, y, 1.2, background],
            bounds=([0, x-2, y-2, 0.7, -np.inf],
                    [np.inf, x+2, y+2, 3.0, np.inf]),
            max_nfev=150)
        a, xc, yc, sigma, bg = fit.x
        chi2 = float(np.sum(fit.fun**2))
        dof = max(1, len(z)-len(fit.x))
        return {"fit_status": "EXPLORATORY_GAUSSIAN_NOT_ASTROMETRY",
                "gauss_x_stamp": float(xc), "gauss_y_stamp": float(yc),
                "gauss_offset_from_wcs_pix": float(np.hypot(xc-x, yc-y)),
                "gauss_sigma_pix": float(sigma), "gauss_amplitude_native": float(a),
                "gauss_chi2_reduced": chi2/dof}
    except Exception as exc:
        return {"fit_status": "GAUSSIAN_FIT_ERROR", "fit_error": str(exc)[:200]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="publication_check/hst_original_flt_1651050")
    args = parser.parse_args()
    root = Path(args.output)
    (root/"panels").mkdir(parents=True, exist_ok=True)
    coord = SkyCoord(RA*u.deg, DEC*u.deg)
    observations = Observations.query_region(coord, radius=3*u.arcsec).to_pandas()
    observations = observations[
        observations.obs_collection.astype(str).str.upper().eq("HST")
    ].copy()
    if "dataRights" in observations:
        observations = observations[
            observations.dataRights.astype(str).str.upper().isin(["PUBLIC", "", "NAN", "NONE"])
        ]
    if observations.empty:
        (root/"SUMMARY.md").write_text("No public HST observations returned. No HST centroid established.\n")
        return
    ids = observations.obsid.dropna().astype("int64").astype(str).drop_duplicates().tolist()
    products = []
    for start in range(0, len(ids), 15):
        tab = Observations.get_product_list(ids[start:start+15])
        if len(tab):
            products.append(tab.to_pandas())
    if not products:
        raise RuntimeError("HST observation inventory exists but product query returned none")
    p = pd.concat(products, ignore_index=True)
    p = p[p.productFilename.astype(str).str.lower().str.endswith(("_flt.fits", "_flc.fits"))].copy()
    if "productType" in p:
        p = p[p.productType.astype(str).str.upper().eq("SCIENCE")]
    p = p.drop_duplicates("dataURI")
    p.to_csv(root/"HST_ORIGINAL_FLT_FLC_PRODUCTS.csv", index=False)
    print("Original HST FLT/FLC products:", len(p), flush=True)
    receipts = []
    for _, row in p.iterrows():
        filename, uri = str(row.productFilename), str(row.dataURI)
        record = {"filename": filename, "dataURI": uri}
        print("Checking HST original pixels", filename, flush=True)
        try:
            with fits.open(mast_download_url(uri), lazy_load_hdus=True, memmap=False,
                    use_fsspec=True, fsspec_kwargs={"block_size": 1024*1024,
                                                     "cache_type": "readahead"}) as hdul:
                record["expstart_mjd"] = float(hdul[0].header.get("EXPSTART", np.nan))
                record["exptime_s"] = float(hdul[0].header.get("EXPTIME", np.nan))
                record["filter"] = str(hdul[0].header.get("FILTER", hdul[0].header.get("FILTER1", "")))
                record["rootname"] = str(hdul[0].header.get("ROOTNAME", filename.split("_")[0]))
                record["detector"] = str(hdul[0].header.get("DETECTOR", ""))
                record["coverage_status"] = "OUTSIDE_DETECTOR"
                for i, hdu in enumerate(hdul):
                    if hdu.name != "SCI":
                        continue
                    nx = int(hdu.header.get("NAXIS1", 0))
                    ny = int(hdu.header.get("NAXIS2", 0))
                    if not nx or not ny:
                        continue
                    try:
                        x, y = WCS(hdu.header).all_world2pix(RA, DEC, 0)
                        x, y = float(x), float(y)
                    except Exception:
                        continue
                    if not (np.isfinite([x,y]).all() and 18 <= x < nx-18 and 18 <= y < ny-18):
                        continue
                    record.update(coverage_status="COVERS_WITH_STAMP", sci_ext=i,
                                  wcs_x=x, wcs_y=y)
                    xc, yc = int(round(x)), int(round(y))
                    x0, y0 = xc-18, yc-18
                    sl = np.s_[y0:yc+19, x0:xc+19]
                    sci = np.asarray(hdu.data[sl], dtype=float)
                    extver = int(hdu.header.get("EXTVER", 1))
                    try:
                        err = np.asarray(hdul["ERR", extver].data[sl], dtype=float)
                    except (KeyError, IndexError):
                        record["measurement_status"] = "NO_MATCHED_ERR_EXTENSION"
                        break
                    try:
                        dq = np.asarray(hdul["DQ", extver].data[sl])
                    except (KeyError, IndexError):
                        dq = np.zeros_like(sci, dtype=np.uint16)
                        record["dq_note"] = "DQ extension missing; no formal detection"
                    tx, ty = x-x0, y-y0
                    measurement = forced_aperture(sci, err, dq, tx, ty)
                    record.update(measurement)
                    if measurement["status"] == "FORCED_PHOTOMETRY_ONLY":
                        record.update(gaussian_screen(sci, err, dq, tx, ty))
                        fig, axes = plt.subplots(1, 2, figsize=(9, 4))
                        bg = measurement["background_native"]
                        axes[0].imshow(sci-bg, origin="lower", cmap="gray")
                        sn = np.divide(sci-bg, err, out=np.full_like(sci, np.nan),
                                       where=(err > 0))
                        axes[1].imshow(sn, origin="lower", cmap="RdBu_r", vmin=-5, vmax=5)
                        for ax in axes:
                            ax.scatter([tx], [ty], marker="x", c="lime", s=90)
                            ax.set_xlim(tx-12, tx+12)
                            ax.set_ylim(ty-12, ty+12)
                        axes[0].set_title("SCI - local sky")
                        axes[1].set_title("Per-pixel S/N (not detection)")
                        fig.suptitle(f"{filename} | forced S/N={measurement['conservative_snr']:.2f}\n"
                                     "HST WCS catalog location; no registered centroid")
                        fig.tight_layout()
                        fig.savefig(root/"panels"/f"{filename.replace('.fits','')}.png", dpi=130)
                        plt.close(fig)
                    break
        except Exception as exc:
            record.update(coverage_status="ANALYSIS_ERROR",
                          error=f"{type(exc).__name__}: {exc}"[:300])
        receipts.append(record)
        pd.DataFrame(receipts).to_csv(root/"HST_ORIGINAL_PIXEL_SCREEN.csv", index=False)
    df = pd.DataFrame(receipts)
    if df.empty:
        ncover = 0
    else:
        ncover = int(df.coverage_status.eq("COVERS_WITH_STAMP").sum())
    text = [
        "# Original HST FLT/FLC screening for 1651050", "",
        f"Public original HST FLT/FLC files examined: {len(df)}.",
        f"Detector-WCS covering images with full 37x37 pixel stamps: {ncover}.",
        "",
        "This is an original-pixel *detectability screen*, not a third-epoch astrometric solution.",
        "Forced-aperture S/N and an exploratory Gaussian fit are NOT source identifications.",
        "HST cross-instrument reference-frame registration, PSF/ePSF comparison,",
        "color terms and a reliable independently detected centroid remain mandatory.",
        "Duplicate filenames/rootnames/timestamps do not count as independent visits.",
        "No 5-parameter motion-plus-parallax fit is attempted.",
    ]
    (root/"SUMMARY.md").write_text("\n".join(text)+"\n")
    print("\n".join(text), flush=True)


if __name__ == "__main__":
    main()
