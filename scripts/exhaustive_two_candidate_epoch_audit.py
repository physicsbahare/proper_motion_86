#!/usr/bin/env python3
"""Exhaustively audit public JWST/NIRCam detections for candidates 210975/575676.

Unlike the ranked-pair PM search, this script does not stop after a fixed number
of epoch pairs.  It visits every public NIRCam observation in the MAST inventory,
resolves CAL products, verifies detector coverage with gWCS, and runs the same
quality-controlled forced photometry/centroid measurement used by the PM
pipeline.  The purpose is to answer a narrower sanity-check question: is the
source actually detected/centroidable in more than one independent epoch, and
if so in which filters/exposures?
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from first26.recenter_measurement import install as install_recenter_patch
from first26.remote_io_patch import install as install_remote_io_patch
from pm86.archive import (
    _get_products_with_retry,
    load_covering_cutouts,
    query_candidate_inventory,
)
from pm86.measurement import measure_exposure
from pm86.pipeline import _group_epochs

TARGET_IDS = (210975, 575676)


def _finite(v):
    try:
        return bool(np.isfinite(float(v)))
    except Exception:
        return False


def _snr(row):
    for key in ("clean_snr_err", "raw_snr_emp", "raw_snr_err"):
        v = row.get(key, np.nan)
        if _finite(v):
            return float(v)
    return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default="data/candidates_first26_090630.csv")
    ap.add_argument("--output", default="results_two_candidate_exhaustive")
    ap.add_argument("--ids", nargs="+", type=int, default=list(TARGET_IDS))
    args = ap.parse_args()

    install_remote_io_patch()
    install_recenter_patch()

    cat = pd.read_csv(args.catalog)
    wanted = cat[cat.candidate_id.astype(int).isin(args.ids)].copy()
    missing = sorted(set(args.ids) - set(wanted.candidate_id.astype(int)))
    if missing:
        raise SystemExit(f"IDs absent from catalog: {missing}")

    outroot = Path(args.output)
    outroot.mkdir(parents=True, exist_ok=True)
    global_summary = []

    for _, crow in wanted.sort_values("source_row").iterrows():
        cid = int(crow.candidate_id)
        ra = float(crow.ra_deg)
        dec = float(crow.dec_deg)
        cdir = outroot / f"candidate_{cid}"
        cdir.mkdir(parents=True, exist_ok=True)

        inv = query_candidate_inventory(cid, ra, dec)
        inv.to_csv(cdir / "archive_inventory_all.csv", index=False)

        coll = inv.get("obs_collection", pd.Series(index=inv.index, dtype=object)).astype(str).str.upper()
        jw = inv[coll.eq("JWST")].copy()
        if "instrument_name" in jw:
            jw = jw[jw.instrument_name.astype(str).str.upper().str.contains("NIRCAM")]
        jw = jw[jw.filter_norm.notna()].copy()
        if "t_min_numeric" not in jw:
            jw["t_min_numeric"] = pd.to_numeric(jw.get("t_min"), errors="coerce")

        epoch_rows, measurement_rows, coverage_rows = [], [], []
        epoch_index = 0
        for filt, frows in jw.groupby("filter_norm"):
            for epoch in _group_epochs(frows):
                epoch_index += 1
                label = f"epoch_{epoch_index:03d}"
                products = _get_products_with_retry(epoch["rows"])
                products.to_csv(cdir / f"{label}_{filt}_products.csv", index=False)
                cutouts, audit = load_covering_cutouts(
                    cid, ra, dec, label, str(filt), products
                )
                for a in audit:
                    a["epoch_mjd_inventory"] = float(epoch["mjd"])
                    a["filter_inventory"] = str(filt)
                    coverage_rows.append(a)

                local = []
                for exp in cutouts:
                    row, _controls = measure_exposure(exp)
                    row["epoch_mjd_inventory"] = float(epoch["mjd"])
                    row["filter_inventory"] = str(filt)
                    row["candidate_id"] = cid
                    local.append(row)
                    measurement_rows.append(row)

                mdf = pd.DataFrame(local)
                if len(mdf):
                    snrs = mdf.apply(_snr, axis=1).to_numpy(float)
                    centroid_ok = (
                        (pd.to_numeric(mdf.get("clean_snr_err"), errors="coerce") >= 3)
                        & np.isfinite(pd.to_numeric(mdf.get("east_2dg_arcsec_raw_wcs"), errors="coerce"))
                    )
                    n_centroid = int(centroid_ok.sum())
                    max_snr = float(np.nanmax(snrs)) if np.isfinite(snrs).any() else np.nan
                else:
                    n_centroid, max_snr = 0, np.nan

                epoch_rows.append({
                    "candidate_id": cid,
                    "epoch_label": label,
                    "filter": str(filt),
                    "epoch_mjd_inventory": float(epoch["mjd"]),
                    "n_observation_rows": int(len(epoch["rows"])),
                    "n_cal_products": int(len(products)),
                    "n_covering_cal_exposures": int(len(cutouts)),
                    "n_secure_centroids_snr3": n_centroid,
                    "max_measurement_snr": max_snr,
                    "secure_detection_epoch": bool(n_centroid > 0),
                })

        epochs = pd.DataFrame(epoch_rows).sort_values("epoch_mjd_inventory") if epoch_rows else pd.DataFrame()
        meas = pd.DataFrame(measurement_rows)
        cov = pd.DataFrame(coverage_rows)
        epochs.to_csv(cdir / "epoch_detection_matrix.csv", index=False)
        meas.to_csv(cdir / "all_exposure_measurements.csv", index=False)
        cov.to_csv(cdir / "all_product_coverage_audit.csv", index=False)

        secure = epochs[epochs.secure_detection_epoch.astype(bool)] if len(epochs) else pd.DataFrame()
        if len(secure) >= 2:
            baseline = float(secure.epoch_mjd_inventory.max() - secure.epoch_mjd_inventory.min())
            status = "MULTI_EPOCH_SECURE_DETECTION"
        elif len(secure) == 1:
            baseline = 0.0
            status = "ONE_SECURE_DETECTION_EPOCH"
        else:
            baseline = np.nan
            status = "NO_SECURE_CENTROID_EPOCH"

        summary = {
            "candidate_id": cid,
            "ra_deg": ra,
            "dec_deg": dec,
            "n_jwst_nircam_inventory_rows": int(len(jw)),
            "n_distinct_filter_epochs": int(len(epochs)),
            "n_secure_detection_epochs": int(len(secure)),
            "secure_detection_baseline_days": None if not np.isfinite(baseline) else baseline,
            "status": status,
            "secure_epochs": secure[["epoch_label", "filter", "epoch_mjd_inventory", "n_secure_centroids_snr3", "max_measurement_snr"]].to_dict("records") if len(secure) else [],
        }
        (cdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        global_summary.append(summary)
        print(json.dumps(summary, indent=2), flush=True)

    pd.DataFrame(global_summary).to_csv(outroot / "two_candidate_summary.csv", index=False)


if __name__ == "__main__":
    main()
