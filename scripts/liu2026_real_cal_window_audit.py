#!/usr/bin/env python3
"""Small real-CAL centroid-window sensitivity check for two Dave UCD candidates.

Uses the frozen epoch pairs and the existing pm86 exposure loader/recentering
method. This is a *detector-level model-sensitivity diagnostic only*: changing
fit windows is NOT an independent proper-motion estimator or a chromatic ePSF.
No catalog/paper result is written or changed.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.stats import sigma_clipped_stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "first26"))
import pm86.measurement as measure
from pm86.archive import query_candidate_inventory, get_pair_products, load_covering_cutouts
from pm86.config import CONFIG
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install
from run_cross_filter_null import choose_frozen_pair, EXPECTED_PAIRS, target_is_usable

WINDOW_HALF_SIZES = (4, 5, 7, 9)


def numerical_row(exposure, baseline, half_size):
    """Compare raw-WCS detector centroid at fixed measured source seed."""
    orig_x = float(baseline.get("x_2dg", np.nan))
    orig_y = float(baseline.get("y_2dg", np.nan))
    if not np.isfinite([orig_x, orig_y]).all():
        return {"half_size": half_size, "status": "NO_BASE_CENTROID"}
    sci = np.asarray(exposure.sci, dtype=float)
    err = np.asarray(exposure.err, dtype=float)
    dq = np.asarray(exposure.dq, dtype=np.uint32)
    phot_bad, astrom_bad = measure.masks(sci, err, dq)
    _, bkg, _ = sigma_clipped_stats(sci, mask=phot_bad, sigma=3., maxiters=5)
    original_config = measure.CONFIG
    try:
        measure.CONFIG = replace(original_config, centroid_half_size=int(half_size))
        fitted = measure.centroid_stamp(sci - float(bkg), err, astrom_bad,
                                       orig_x, orig_y, n_mc=0)
    finally:
        measure.CONFIG = original_config
    out = {"half_size": int(half_size), "status": "FIT_FAILED",
           "x_2dg": float(fitted["x_2dg"]), "y_2dg": float(fitted["y_2dg"]),
           "x_com": float(fitted["x_com"]), "y_com": float(fitted["y_com"])}
    for method in ("2dg", "com"):
        x, y = out["x_" + method], out["y_" + method]
        if not np.isfinite([x, y]).all():
            continue
        ra, dec = measure.local_pixel_to_sky(exposure, x, y)
        east, north = measure.sky_to_tangent(
            ra, dec, exposure.ra_deg, exposure.dec_deg)
        out["east_" + method + "_mas"] = float(east) * 1000
        out["north_" + method + "_mas"] = float(north) * 1000
        out["delta_east_" + method + "_mas"] = (
            float(east) * 1000 - float(baseline["east_2dg_arcsec_raw_wcs"]) * 1000)
        out["delta_north_" + method + "_mas"] = (
            float(north) * 1000 - float(baseline["north_2dg_arcsec_raw_wcs"]) * 1000)
        if method == "2dg":
            out["status"] = "CENTROID_OK"
    return out


def summarize_candidate(cid, output):
    out = output / str(cid)
    out.mkdir(parents=True, exist_ok=True)
    cat = pd.read_csv(ROOT / "data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    frozen = pd.read_csv(ROOT / "data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    target = cat.loc[cid]
    original = frozen.loc[cid]
    ra, dec = float(target.ra_deg), float(target.dec_deg)
    filters = EXPECTED_PAIRS[cid]
    pair = choose_frozen_pair(query_candidate_inventory(cid, ra, dec),
                              filters, float(original.baseline_days))
    ep, lp = get_pair_products(pair)
    early, a = load_covering_cutouts(cid, ra, dec, "early", filters[0], ep)
    late, b = load_covering_cutouts(cid, ra, dec, "late", filters[1], lp)
    pd.DataFrame(a + b).to_csv(out / "coverage.csv", index=False)
    observations, windows = [], []
    for e in early + late:
        obs, controls = measure_exposure_recentered(e)
        usable = bool(target_is_usable(obs))
        r = {k: v for k, v in obs.items()
             if k in ("filename", "filter", "filter_name", "mjd", "epoch_label",
                      "astrometric_snr", "clean_snr_err", "astrometry_source",
                      "forced_fit_attempted", "forced_fit_accepted",
                      "x_2dg", "y_2dg", "x_com", "y_com",
                      "dq_center", "dq_or_r4", "dq_jump_det",
                      "east_2dg_arcsec_raw_wcs", "north_2dg_arcsec_raw_wcs")}
        r.update({"filename": e.filename, "filter": e.filter_name,
                  "epoch": e.epoch_label, "usable_frozen_gate": usable,
                  "n_field_controls": len(controls)})
        observations.append(r)
        if not usable:
            continue
        for h in WINDOW_HALF_SIZES:
            try:
                w = numerical_row(e, obs, h)
                w.update({"filename": e.filename, "filter": e.filter_name,
                          "epoch": e.epoch_label, "target_id": cid})
                windows.append(w)
            except Exception as error:
                windows.append({"filename": e.filename, "half_size": h,
                                "status": "ERROR:" + repr(error)})
    pd.DataFrame(observations).to_csv(out / "original_target_quality.csv", index=False)
    pd.DataFrame(windows).to_csv(out / "window_centroids.csv", index=False)

    summary = {"target_id": cid, "status": "DIAGNOSTIC_ONLY",
               "frozen_filters": list(filters),
               "n_covering": {"early": len(early), "late": len(late)},
               "n_target_usable": {
                   "early": sum(x["usable_frozen_gate"] for x in observations if x["epoch"] == "early"),
                   "late": sum(x["usable_frozen_gate"] for x in observations if x["epoch"] == "late")},
               "per_exposure": [],
               "restrictions": [
                   "Raw-WCS within-exposure differences only; NOT a new registered PM.",
                   "Window sizes change the Gaussian/COM estimator; NO color-matched PSF validation.",
                   "Cannot resolve two-epoch parallax/intrinsic-motion degeneracy.",
                   "No significance or publication classification inferred from this sensitivity test."
               ]}
    for row in observations:
        if not row["usable_frozen_gate"]:
            continue
        fits = [w for w in windows if w.get("filename") == row["filename"]
                and w.get("status") == "CENTROID_OK"
                and np.isfinite(w.get("delta_east_2dg_mas", np.nan))
                and np.isfinite(w.get("delta_north_2dg_mas", np.nan))]
        spans = [float(np.ptp([w["east_2dg_mas"] for w in fits])) if fits else None,
                 float(np.ptp([w["north_2dg_mas"] for w in fits])) if fits else None]
        baseline_com = np.nan
        if np.isfinite(row.get("x_com", np.nan)) and np.isfinite(row.get("y_com", np.nan)):
            baseline_com = float(np.hypot(row["x_com"] - row["x_2dg"],
                                          row["y_com"] - row["y_2dg"])) * 63.0
        summary["per_exposure"].append({
            "filename": row["filename"], "filter": row["filter"],
            "epoch": row["epoch"], "snr": row.get("astrometric_snr"),
            "n_fit_windows_success": len(fits),
            "gaussian_window_span_east_mas": spans[0],
            "gaussian_window_span_north_mas": spans[1],
            "com_vs_gaussian_approx_separation_mas": baseline_com,
        })
    (out / "summary.json").write_text(json.dumps(summary, indent=2,
                                               default=lambda x: x.item() if isinstance(x, np.generic) else str(x)))
    print("AUDIT_SUMMARY " + json.dumps(summary, default=str), flush=True)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", type=int, required=True, choices=(387132,1651050))
    ap.add_argument("--out", default="results_liu2026_window_audit")
    a = ap.parse_args()
    install()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    try:
        result = summarize_candidate(a.candidate, out)
        n = result["n_target_usable"]
        if not n["early"] or not n["late"]:
            raise RuntimeError("NO_USABLE_BOTH_EPOCHS")
    except Exception:
        dest = out / str(a.candidate)
        dest.mkdir(parents=True, exist_ok=True)
        (dest / "failure.txt").write_text(traceback.format_exc())
        raise


if __name__ == "__main__":
    main()
