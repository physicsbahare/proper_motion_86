#!/usr/bin/env python3
"""Reproduce the two Dave epoch pairs and measure cross-filter control-source nulls.

Usage: python scripts/run_cross_filter_null.py --out results_cross_filter_null
Reads only public MAST Stage-2 CAL cutouts. Original paper numbers are never
overwritten. This run is a preliminary method-sensitivity audit, NOT a certified
PSF/chromatic correction or a new proper-motion measurement.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "first26"))

from pm86.archive import get_pair_products, load_covering_cutouts, query_candidate_inventory
from pm86.astrometry import infer_proper_motion, register_exposures
from pm86.cross_filter_null import (
    conservative_sensitivity, leave_one_out_null, summarize_null,
)
from pm86.pipeline import _rank_epoch_pairs
from pm86.refit_reconciliation import reconcile_frozen_refit
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install as install_bounded_remote

FROZEN_CANDIDATE_IDS = (387132, 1651050)
EXPECTED_PAIRS = {
    387132: ("F444W", "F410M"),
    1651050: ("F444W", "F356W"),
}


def plain(obj):
    if isinstance(obj, dict):
        return {str(k): plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [plain(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return plain(float(obj))
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def dump(path: Path, obj):
    path.write_text(json.dumps(plain(obj), indent=2, sort_keys=True), encoding="utf-8")


def choose_frozen_pair(inventory: pd.DataFrame, expected: tuple[str, str],
                       baseline_days: float) -> dict:
    pairs = _rank_epoch_pairs(inventory, max_pairs=100000)
    relevant = [p for p in pairs if
                (p["filter_early"], p["filter_late"]) == expected]
    if not relevant:
        raise RuntimeError(f"Cannot locate public {expected} epoch pair")
    best = min(relevant, key=lambda x: abs(x["baseline_days_inventory"] - baseline_days))
    if abs(best["baseline_days_inventory"] - baseline_days) > 5:
        raise RuntimeError("No matching epoch pair within 5 days of frozen measured baseline")
    return best


def target_is_usable(m: dict) -> bool:
    snr = float(m.get("astrometric_snr", m.get("clean_snr_err", np.nan)))
    return bool(np.isfinite(snr) and snr >= 3 and
                np.isfinite(float(m.get("east_2dg_arcsec_raw_wcs", np.nan))) and
                np.isfinite(float(m.get("north_2dg_arcsec_raw_wcs", np.nan))))


def run_candidate(cid: int, catalog_row: pd.Series, frozen: pd.Series,
                  destination: Path) -> dict:
    destination.mkdir(parents=True, exist_ok=True)
    ra, dec = float(catalog_row.ra_deg), float(catalog_row.dec_deg)
    filters = EXPECTED_PAIRS[cid]
    outcome = {
        "candidate_id": cid,
        "catalog_ra_deg": ra, "catalog_dec_deg": dec,
        "expected_filters": filters,
        "published_baseline_days": float(frozen.baseline_days),
        "published_nominal_2d_score": float(frozen.significance_2d),
        "status": "STARTED",
        "provenance": "data/two_pm_candidates_125.csv and public MAST original CAL cutouts",
        "calibration_type": "preliminary_exposure_pair_leave_one_out_field_control_null",
    }
    try:
        inv = query_candidate_inventory(cid, ra, dec)
        pair = choose_frozen_pair(inv, filters, float(frozen.baseline_days))
        outcome["selected_inventory_pair"] = {
            "filter_early": pair["filter_early"],
            "filter_late": pair["filter_late"],
            "inventory_early_mjd": float(pair["early"]["mjd"]),
            "inventory_late_mjd": float(pair["late"]["mjd"]),
            "inventory_baseline_days": float(pair["baseline_days_inventory"]),
        }
        early_products, late_products = get_pair_products(pair)
        for epoch, products in [("early", early_products), ("late", late_products)]:
            products.to_csv(destination / f"{epoch}_products.csv", index=False)
        early, early_audit = load_covering_cutouts(cid, ra, dec, "early", filters[0], early_products)
        late, late_audit = load_covering_cutouts(cid, ra, dec, "late", filters[1], late_products)
        pd.DataFrame(early_audit + late_audit).to_csv(destination / "product_coverage_audit.csv", index=False)
        outcome["n_covering_early"] = len(early)
        outcome["n_covering_late"] = len(late)
        if not early or not late:
            raise RuntimeError("No CAL cutout covers candidate in one of the selected epochs")

        records = []
        controls = {}
        for exposure in early + late:
            measurement, ctrl = measure_exposure_recentered(exposure)
            records.append(measurement)
            controls[exposure.filename] = ctrl
        measurements = pd.DataFrame(records)
        measurements.to_csv(destination / "exposure_measurements.csv", index=False)
        outcome["n_target_usable_early"] = int(sum(target_is_usable(r) for r in records if r["epoch_label"] == "early"))
        outcome["n_target_usable_late"] = int(sum(target_is_usable(r) for r in records if r["epoch_label"] == "late"))

        # Refit the target using ALL usable exposures as a reproducibility check.
        # It remains separate from the frozen manuscript values.
        all_ctrl = pd.concat([c for c in controls.values() if len(c)], ignore_index=True) if any(len(c) for c in controls.values()) else pd.DataFrame()
        if not all_ctrl.empty:
            try:
                reg, qual, matched = register_exposures(measurements, all_ctrl)
                reg.to_csv(destination / "registered_target_positions.csv", index=False)
                qual.to_csv(destination / "registration_quality.csv", index=False)
                matched.to_csv(destination / "registration_control_residuals.csv", index=False)
                if not reg.empty:
                    pm, pairwise = infer_proper_motion(reg, qual, "CROSS_FILTER_JWST", matched)
                    dump(destination / "independent_target_refit.json", pm)
                    reconciliation = reconcile_frozen_refit(frozen.to_dict(), pm, measurements, reg)
                    dump(destination / "frozen_refit_reconciliation.json", reconciliation)
                    outcome["refit_reconciliation"] = {k: v for k, v in reconciliation.items() if k != "all_exposure_diagnostics"}
                    pairwise.to_csv(destination / "target_pairwise_diagnostics.csv", index=False)
                    outcome["target_refit"] = pm
            except Exception as exc:
                outcome["target_refit_status"] = f"UNAVAILABLE:{type(exc).__name__}:{exc}"

        # Choose one target-centroidable early/late exposure pair. Do not treat
        # dither combinations as independent null measurements.
        early_rows = [r for r in records if r["epoch_label"] == "early"]
        late_rows = [r for r in records if r["epoch_label"] == "late"]
        candidates = []
        for a in early_rows:
            for b in late_rows:
                if not (target_is_usable(a) and target_is_usable(b)):
                    continue
                na, nb = len(controls[a["filename"]]), len(controls[b["filename"]])
                snr_a = float(a.get("astrometric_snr", a.get("clean_snr_err", np.nan)))
                snr_b = float(b.get("astrometric_snr", b.get("clean_snr_err", np.nan)))
                candidates.append((min(na, nb), min(snr_a, snr_b), a, b))
        if not candidates:
            raise RuntimeError("No target-centroidable exposure pair; cannot calibrate same target solution")
        _, _, a, b = max(candidates, key=lambda z: z[:2])
        outcome["representative_exposure_pair"] = {
            "early_filename": a["filename"], "late_filename": b["filename"],
            "early_snr": a.get("astrometric_snr", a.get("clean_snr_err")),
            "late_snr": b.get("astrometric_snr", b.get("clean_snr_err")),
            "not_independent_exposure_ensemble": True,
        }
        residuals, qc = leave_one_out_null(controls[a["filename"]], controls[b["filename"]])
        residuals.to_csv(destination / "leave_one_out_residuals.csv", index=False)
        outcome["null_qc"] = qc
        summary = summarize_null(residuals)
        outcome["null_summary"] = summary
        if summary.get("status") == "PRELIMINARY_FIELD_CONTROL_NULL":
            outcome["sensitivity"] = conservative_sensitivity(
                float(frozen.mu_alpha_cosdec_masyr),
                float(frozen.mu_delta_masyr),
                float(frozen.sigma_mu_alpha_cosdec_masyr),
                float(frozen.sigma_mu_delta_masyr),
                float(frozen.baseline_days),
                summary,
            )
            # A raw segmented flux ratio is an instrumental proxy, NOT AB color.
            ea, eb = float(a.get("clean_flux", np.nan)), float(b.get("clean_flux", np.nan))
            if ea > 0 and eb > 0 and "control_instrumental_color_min" in summary:
                proxy = -2.5 * math.log10(eb / ea)
                outcome["target_instrumental_color_proxy"] = proxy
                outcome["target_proxy_within_controls"] = bool(
                    summary["control_instrumental_color_min"] <= proxy
                    <= summary["control_instrumental_color_max"]
                )
            else:
                outcome["target_proxy_within_controls"] = None
        outcome["status"] = "PRELIMINARY_NULL_AVAILABLE" if summary.get("status") == "PRELIMINARY_FIELD_CONTROL_NULL" else "INSUFFICIENT_FIELD_CONTROL_NULL"
    except Exception as exc:
        outcome["status"] = "UNABLE_TO_CALIBRATE"
        outcome["error"] = f"{type(exc).__name__}: {exc}"
        (destination / "traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
    dump(destination / "summary.json", outcome)
    return outcome


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results_cross_filter_null")
    args = ap.parse_args()
    install_bounded_remote()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    catalog = pd.read_csv(ROOT / "data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    frozen = pd.read_csv(ROOT / "data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    summaries = []
    for cid in FROZEN_CANDIDATE_IDS:
        if cid not in catalog.index or cid not in frozen.index:
            raise RuntimeError(f"Candidate {cid} absent from frozen input catalogs")
        print(f"Running cross-filter null for candidate {cid}", flush=True)
        result = run_candidate(cid, catalog.loc[cid], frozen.loc[cid], outdir / str(cid))
        summaries.append(result)
        print(f"{cid}: {result['status']}", flush=True)
    dump(outdir / "summary.json", {
        "method": "leave_one_out_local_registration_field_control_null",
        "results": summaries,
        "scientific_status": "PRELIMINARY_ONLY_NO_SECURE_PM_RECLASSIFICATION",
        "limitations": [
            "Only one representative exposure pair enters the empirical null.",
            "Field control galaxies are not morphologically vetted; stars may move.",
            "Control centroid noise and registration scatter already occur in nominal PM errors; quadrature sensitivity may double-count.",
            "Instrumental color proxy is not calibrated AB color and may not represent ultracool dwarf spectra.",
            "No ePSF injection-recovery, source-color chromatic correction, or stationary-population null calibration.",
            "Do not update nominal paper significance or claim final Gaussian sigma from this output."
        ],
    })
    pd.DataFrame([{
        "candidate_id": r["candidate_id"],
        "status": r["status"],
        "n_null": r.get("null_summary", {}).get("n"),
        "bias_east_mas": r.get("null_summary", {}).get("bias_east_mas"),
        "bias_north_mas": r.get("null_summary", {}).get("bias_north_mas"),
        "scatter_east_mas": r.get("null_summary", {}).get("scatter_east_mas"),
        "scatter_north_mas": r.get("null_summary", {}).get("scatter_north_mas"),
        "conservative_2d_score": r.get("sensitivity", {}).get("conservative_2d_mahalanobis_score"),
    } for r in summaries]).to_csv(outdir / "summary_table.csv", index=False)
    if not any(r["status"] == "PRELIMINARY_NULL_AVAILABLE" for r in summaries):
        raise SystemExit("No empirical null obtained: examine candidate artifact and MAST accessibility")


if __name__ == "__main__":
    main()
