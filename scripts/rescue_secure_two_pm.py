#!/usr/bin/env python3
"""Exhaustive proper-motion rescue for secure candidates 835615 and 1140578.

These two objects were independently labelled ``confirmed_secure`` by the
single-exposure audit but the standard 45-candidate PM run returned
INSUFFICIENT_DATA.  This script therefore removes the ranked-pair cap and tests
*every* scientifically admissible independent JWST/NIRCam epoch pair.

For every pair it:
  * resolves Stage-2 CAL products and verifies true detector coverage with gWCS;
  * runs the quality-controlled target measurement plus the faint-source forced
    Gaussian/recentring rescue used by the first26 audit;
  * records DQ/artifact diagnostics and target S/N per exposure;
  * requires a usable centroid in both independent epochs;
  * performs local-field affine registration using control sources;
  * infers PM with propagated centroid + registration uncertainties;
  * stores pairwise 2DG/COM PM measurements and registration quality;
  * compares all successful pair solutions for cross-pair consistency.

No external team's PM labels/measurements are used as inputs.  The physical
inputs are only the catalog coordinates and public MAST imaging.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "first26"))
from recenter_measurement import install as install_recenter_patch
from remote_io_patch import install as install_remote_io_patch

from pm86.archive import get_pair_products, load_covering_cutouts, query_candidate_inventory
from pm86.astrometry import infer_proper_motion, register_exposures
from pm86.measurement import measure_exposure
from pm86.pipeline import _pair_metadata, _rank_epoch_pairs

TARGET_IDS = (835615, 1140578)


def finite(v):
    try:
        return bool(np.isfinite(float(v)))
    except Exception:
        return False


def target_snr(row: pd.Series) -> float:
    """Use the astrometric forced-fit S/N when the rescue patch provides it."""
    for key in ("astrometric_snr", "clean_snr_err", "raw_snr_emp", "raw_snr_err"):
        v = row.get(key, np.nan)
        if finite(v):
            return float(v)
    return np.nan


def write_json(path: Path, obj):
    def clean(x):
        if isinstance(x, dict):
            return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        if isinstance(x, (np.integer,)):
            return int(x)
        if isinstance(x, (np.floating,)):
            x = float(x)
        if isinstance(x, float) and not np.isfinite(x):
            return None
        if isinstance(x, (np.bool_,)):
            return bool(x)
        return x
    path.write_text(json.dumps(clean(obj), indent=2, sort_keys=True), encoding="utf-8")


def pair_success_score(row: dict) -> tuple:
    """Prefer same-filter, long-baseline, precise, well-registered solutions."""
    same = 1 if row.get("pair_type") == "SAME_FILTER_JWST" else 0
    sig = float(row.get("significance_2d", -np.inf)) if finite(row.get("significance_2d")) else -np.inf
    baseline = float(row.get("baseline_days", 0.0)) if finite(row.get("baseline_days")) else 0.0
    se = float(row.get("sigma_mu_alpha_cosdec_masyr", np.inf)) if finite(row.get("sigma_mu_alpha_cosdec_masyr")) else np.inf
    sn = float(row.get("sigma_mu_delta_masyr", np.inf)) if finite(row.get("sigma_mu_delta_masyr")) else np.inf
    precision = -(math.hypot(se, sn) if np.isfinite(se) and np.isfinite(sn) else np.inf)
    return (same, baseline, precision, sig)


def cross_pair_consistency(success: pd.DataFrame) -> dict:
    if len(success) == 0:
        return {"n_successful_pairs": 0, "status": "NO_PM_SOLUTION"}
    if len(success) == 1:
        return {"n_successful_pairs": 1, "status": "ONE_PM_SOLUTION"}

    mua = pd.to_numeric(success["mu_alpha_cosdec_masyr"], errors="coerce").to_numpy(float)
    mud = pd.to_numeric(success["mu_delta_masyr"], errors="coerce").to_numpy(float)
    sea = pd.to_numeric(success["sigma_mu_alpha_cosdec_masyr"], errors="coerce").to_numpy(float)
    sed = pd.to_numeric(success["sigma_mu_delta_masyr"], errors="coerce").to_numpy(float)
    good = np.isfinite(mua) & np.isfinite(mud) & np.isfinite(sea) & np.isfinite(sed) & (sea > 0) & (sed > 0)
    if good.sum() < 2:
        return {"n_successful_pairs": int(len(success)), "status": "MULTIPLE_SOLUTIONS_UNQUANTIFIED"}

    wa = 1.0 / sea[good] ** 2
    wd = 1.0 / sed[good] ** 2
    ma = float(np.sum(wa * mua[good]) / np.sum(wa))
    md = float(np.sum(wd * mud[good]) / np.sum(wd))
    chi2 = float(np.sum(((mua[good] - ma) / sea[good]) ** 2 + ((mud[good] - md) / sed[good]) ** 2))
    dof = max(1, 2 * int(good.sum()) - 2)
    red = chi2 / dof
    return {
        "n_successful_pairs": int(len(success)),
        "n_quantified_pairs": int(good.sum()),
        "weighted_mu_alpha_cosdec_masyr": ma,
        "weighted_mu_delta_masyr": md,
        "cross_pair_chi2": chi2,
        "cross_pair_dof": dof,
        "cross_pair_reduced_chi2": red,
        "status": "CROSS_PAIR_CONSISTENT" if red <= 3.0 else "CROSS_PAIR_TENSION",
    }


def copy_best(adir: Path, cdir: Path):
    for name in (
        "pair.json", "products_early.csv", "products_late.csv", "product_coverage_audit.csv",
        "exposure_measurements.csv", "field_controls.csv", "registered_positions.csv",
        "registration_quality.csv", "registration_control_residuals.csv", "pairwise_pm.csv",
        "pm.json",
    ):
        src = adir / name
        if src.exists():
            shutil.copy2(src, cdir / f"best_{name}")


def audit_candidate(crow: pd.Series, outroot: Path) -> dict:
    cid = int(crow.candidate_id)
    ra, dec = float(crow.ra_deg), float(crow.dec_deg)
    cdir = outroot / f"candidate_{cid}"
    cdir.mkdir(parents=True, exist_ok=True)

    inv = query_candidate_inventory(cid, ra, dec)
    inv.to_csv(cdir / "archive_inventory.csv", index=False)

    # 100000 is effectively uncapped for a single COSMOS coordinate.
    pairs = _rank_epoch_pairs(inv, max_pairs=100000)
    write_json(cdir / "inventory_summary.json", {
        "candidate_id": cid,
        "ra_deg": ra,
        "dec_deg": dec,
        "catalog_f444w_mag_ab": float(crow.f444w_mag_ab),
        "n_inventory_rows": int(len(inv)),
        "n_all_independent_nircam_pairs": int(len(pairs)),
    })

    attempts = []
    successful = []
    all_measurements = []
    best = None

    for rank, pair in enumerate(pairs, 1):
        adir = cdir / "all_pair_attempts" / f"pair_{rank:04d}"
        adir.mkdir(parents=True, exist_ok=True)
        meta = {**_pair_metadata(pair), "attempt_rank": rank}
        write_json(adir / "pair.json", meta)
        base = dict(meta)
        try:
            ep, lp = get_pair_products(pair)
            ep.to_csv(adir / "products_early.csv", index=False)
            lp.to_csv(adir / "products_late.csv", index=False)
            ec, ea = load_covering_cutouts(cid, ra, dec, "early", pair["filter_early"], ep)
            lc, la = load_covering_cutouts(cid, ra, dec, "late", pair["filter_late"], lp)
            pd.DataFrame(ea + la).to_csv(adir / "product_coverage_audit.csv", index=False)

            base.update(n_covering_early=len(ec), n_covering_late=len(lc))
            if not ec or not lc:
                attempts.append({**base, "stage": "COVERAGE_FAIL", "reason": "no covering CAL exposure in one epoch"})
                continue

            mrows, cframes = [], []
            for exp in ec + lc:
                m, controls = measure_exposure(exp)
                m["pair_attempt_rank"] = rank
                mrows.append(m)
                all_measurements.append(m.copy())
                if len(controls):
                    cframes.append(controls)
            mdf = pd.DataFrame(mrows).sort_values("mjd").reset_index(drop=True)
            cdf = pd.concat(cframes, ignore_index=True) if cframes else pd.DataFrame()
            mdf.to_csv(adir / "exposure_measurements.csv", index=False)
            cdf.to_csv(adir / "field_controls.csv", index=False)

            snr = mdf.apply(target_snr, axis=1)
            east = pd.to_numeric(mdf.get("east_2dg_arcsec_raw_wcs"), errors="coerce")
            north = pd.to_numeric(mdf.get("north_2dg_arcsec_raw_wcs"), errors="coerce")
            usable = (snr >= 3.0) & np.isfinite(east) & np.isfinite(north)
            ue = usable & mdf.epoch_label.eq("early")
            ul = usable & mdf.epoch_label.eq("late")
            base.update(
                n_astrometric_early=int(ue.sum()),
                n_astrometric_late=int(ul.sum()),
                max_snr_early=float(np.nanmax(snr[mdf.epoch_label.eq("early")])) if mdf.epoch_label.eq("early").any() else np.nan,
                max_snr_late=float(np.nanmax(snr[mdf.epoch_label.eq("late")])) if mdf.epoch_label.eq("late").any() else np.nan,
            )
            if not ue.any() or not ul.any():
                attempts.append({**base, "stage": "CENTROID_FAIL", "reason": "no secure centroid in both independent epochs"})
                continue

            reg, qual, matched = register_exposures(mdf, cdf)
            reg.to_csv(adir / "registered_positions.csv", index=False)
            qual.to_csv(adir / "registration_quality.csv", index=False)
            matched.to_csv(adir / "registration_control_residuals.csv", index=False)
            if reg.empty:
                attempts.append({**base, "stage": "REGISTRATION_FAIL", "reason": "no registered target positions"})
                continue

            pm, pairwise = infer_proper_motion(reg, qual, pair["pair_type"])
            pairwise.to_csv(adir / "pairwise_pm.csv", index=False)
            write_json(adir / "pm.json", pm)
            classification = pm.get("classification") or "INSUFFICIENT_DATA"
            if classification == "INSUFFICIENT_DATA":
                attempts.append({**base, "stage": "PM_INFERENCE_FAIL", "reason": pm.get("reason", "PM inference insufficient")})
                continue

            result = {**base, **pm, "stage": "PM_MEASURED"}
            successful.append(result)
            attempts.append(result)
            if best is None or pair_success_score(result) > pair_success_score(best[0]):
                best = (result, adir)

        except Exception as exc:
            attempts.append({**base, "stage": "ERROR", "reason": f"{type(exc).__name__}: {exc}"})

    adf = pd.DataFrame(attempts)
    adf.to_csv(cdir / "all_pair_results.csv", index=False)
    if all_measurements:
        pd.DataFrame(all_measurements).to_csv(cdir / "all_exposure_measurements_across_pairs.csv", index=False)
    sdf = pd.DataFrame(successful)
    sdf.to_csv(cdir / "successful_pm_solutions.csv", index=False)
    consistency = cross_pair_consistency(sdf)
    write_json(cdir / "cross_pair_consistency.json", consistency)

    stage_counts = adf.stage.value_counts().to_dict() if len(adf) else {}
    if best is None:
        summary = {
            "candidate_id": cid,
            "ra_deg": ra,
            "dec_deg": dec,
            "catalog_f444w_mag_ab": float(crow.f444w_mag_ab),
            "n_independent_pairs_tested": int(len(pairs)),
            "n_successful_pm_solutions": 0,
            "classification": "INSUFFICIENT_DATA",
            "pm_status": "INSUFFICIENT_DATA",
            "reason": "exhaustive all-pair rescue found no defensible two-epoch PM solution",
            "attempt_stage_counts": stage_counts,
            "cross_pair_consistency": consistency,
        }
    else:
        b, bdir = best
        copy_best(bdir, cdir)
        summary = {
            "candidate_id": cid,
            "ra_deg": ra,
            "dec_deg": dec,
            "catalog_f444w_mag_ab": float(crow.f444w_mag_ab),
            "n_independent_pairs_tested": int(len(pairs)),
            "n_successful_pm_solutions": int(len(sdf)),
            "best_pair_rank": int(b["attempt_rank"]),
            "best_pair_type": b.get("pair_type"),
            "best_filter_early": b.get("filter_early"),
            "best_filter_late": b.get("filter_late"),
            "classification": b.get("classification"),
            "pm_status": b.get("pm_status"),
            "mu_alpha_cosdec_masyr": b.get("mu_alpha_cosdec_masyr"),
            "mu_delta_masyr": b.get("mu_delta_masyr"),
            "sigma_mu_alpha_cosdec_masyr": b.get("sigma_mu_alpha_cosdec_masyr"),
            "sigma_mu_delta_masyr": b.get("sigma_mu_delta_masyr"),
            "significance_2d": b.get("significance_2d"),
            "baseline_days": b.get("baseline_days"),
            "attempt_stage_counts": stage_counts,
            "cross_pair_consistency": consistency,
        }
    write_json(cdir / "summary.json", summary)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default="data/candidates_45.csv")
    ap.add_argument("--ids", nargs="+", type=int, default=list(TARGET_IDS))
    ap.add_argument("--output", default="results_secure_two_rescue")
    args = ap.parse_args()

    install_remote_io_patch()
    install_recenter_patch()

    cat = pd.read_csv(args.catalog)
    sub = cat[cat.candidate_id.astype(int).isin(args.ids)].copy()
    missing = sorted(set(args.ids) - set(sub.candidate_id.astype(int)))
    if missing:
        raise SystemExit(f"candidate IDs missing from catalog: {missing}")

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    summaries = []
    for _, row in sub.sort_values("candidate_id").iterrows():
        s = audit_candidate(row, out)
        summaries.append(s)
        print(json.dumps(s, indent=2, default=str), flush=True)
    pd.DataFrame(summaries).to_csv(out / "secure_two_rescue_summary.csv", index=False)


if __name__ == "__main__":
    main()
