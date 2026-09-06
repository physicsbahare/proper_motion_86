"""Epoch-cached exhaustive PM audit for candidates with many legal epoch pairs.

Scientifically this is the same search as :mod:`defensible.exhaustive_pm`: every legal
independent JWST/NIRCam epoch pair is still tested.  The difference is purely
computational.  Archive product discovery, CAL cutout loading, DQ-aware target
measurement, and field-control measurement are performed once per distinct epoch and
then reused across all pairs containing that epoch.  This removes the O(N_pairs)
repetition that caused dense candidates to hit GitHub Actions time limits.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import pm86.archive as archive
import pm86.astrometry as astrometry

from defensible.exhaustive_pm import (
    enumerate_all_pairs,
    measure_exposure_dq_aware,
    _load_all_covering,
)


def _epoch_key(epoch: dict) -> str:
    return str(epoch["epoch_id"])


def _prepare_epoch(candidate_id: int, ra_deg: float, dec_deg: float, epoch: dict, outdir: Path) -> dict:
    """Discover/load/measure one physical epoch exactly once."""
    key = _epoch_key(epoch)
    edir = outdir / "epochs" / key
    edir.mkdir(parents=True, exist_ok=True)

    products = archive._get_products_with_retry(epoch["rows"])
    products.to_csv(edir / "products.csv", index=False)

    exposures, coverage = _load_all_covering(
        candidate_id, ra_deg, dec_deg, key, str(epoch["filter"]), products
    )
    pd.DataFrame(coverage).to_csv(edir / "product_coverage_audit.csv", index=False)

    rows, ctrs = [], []
    for exp in exposures:
        row, controls = measure_exposure_dq_aware(exp)
        rows.append(row)
        if len(controls):
            ctrs.append(controls)

    meas = pd.DataFrame(rows)
    controls = pd.concat(ctrs, ignore_index=True) if ctrs else pd.DataFrame()
    meas.to_csv(edir / "exposure_measurements.csv", index=False)
    controls.to_csv(edir / "field_controls.csv", index=False)

    n_astrom = 0
    if not meas.empty and "astrometrically_usable" in meas:
        n_astrom = int(meas["astrometrically_usable"].astype(bool).sum())

    return {
        "epoch_id": key,
        "filter": str(epoch["filter"]),
        "mjd": float(epoch["mjd"]),
        "n_covering": len(exposures),
        "n_astrometric": n_astrom,
        "measurements": meas,
        "controls": controls,
    }


def _relabel_epoch(df: pd.DataFrame, label: str) -> pd.DataFrame:
    out = df.copy()
    if not out.empty and "epoch_label" in out.columns:
        out["epoch_label"] = label
    return out


def _audit_pair_from_cache(candidate_id: int, pair: dict, cache: dict[str, dict], outdir: Path) -> dict:
    pdir = outdir / "pairs" / pair["pair_id"]
    pdir.mkdir(parents=True, exist_ok=True)

    early = cache[_epoch_key(pair["early"])]
    late = cache[_epoch_key(pair["late"])]
    result = {
        "candidate_id": candidate_id,
        "pair_id": pair["pair_id"],
        "pair_type": pair["pair_type"],
        "filter_early": pair["filter_early"],
        "filter_late": pair["filter_late"],
        "inventory_baseline_days": pair["baseline_days_inventory"],
        "n_covering_early": int(early["n_covering"]),
        "n_covering_late": int(late["n_covering"]),
        "n_astrometric_early": int(early["n_astrometric"]),
        "n_astrometric_late": int(late["n_astrometric"]),
    }

    if early["n_covering"] == 0 or late["n_covering"] == 0:
        result.update(pm_status="INSUFFICIENT_DATA", stage="COVERAGE_FAIL", reason="epochs do not both have covering CAL exposures")
        return result
    if early["n_astrometric"] == 0 or late["n_astrometric"] == 0:
        result.update(pm_status="INSUFFICIENT_DATA", stage="CENTROID_FAIL", reason="target is not securely centroidable in both independent epochs")
        return result

    meas = pd.concat([
        _relabel_epoch(early["measurements"], "early"),
        _relabel_epoch(late["measurements"], "late"),
    ], ignore_index=True)
    ctr_frames = []
    for block, label in ((early["controls"], "early"), (late["controls"], "late")):
        if not block.empty:
            ctr_frames.append(_relabel_epoch(block, label))
    controls = pd.concat(ctr_frames, ignore_index=True) if ctr_frames else pd.DataFrame()

    # Save pair-specific joined tables only when the pair reaches registration; the
    # per-epoch files above retain the complete measurement audit for every epoch.
    try:
        reg, qual, matched = astrometry.register_exposures(meas, controls)
    except Exception as exc:
        result.update(pm_status="INSUFFICIENT_DATA", stage="REGISTRATION_FAIL", reason=f"{type(exc).__name__}: {exc}")
        return result

    reg.to_csv(pdir / "registered_target_positions.csv", index=False)
    qual.to_csv(pdir / "registration_quality.csv", index=False)
    matched.to_csv(pdir / "matched_controls.csv", index=False)

    summary, pairwise = astrometry.infer_proper_motion(reg, qual, pair["pair_type"])
    pairwise.to_csv(pdir / "pairwise_pm_diagnostics.csv", index=False)
    result.update(stage="PM_FIT", **summary)
    return result


def run_candidate_cached(candidate_id: int, ra_deg: float, dec_deg: float, output_root: str | Path) -> dict:
    outdir = Path(output_root) / f"candidate_{candidate_id}"
    outdir.mkdir(parents=True, exist_ok=True)

    inventory = archive.query_candidate_inventory(candidate_id, ra_deg, dec_deg)
    inventory.to_csv(outdir / "archive_inventory.csv", index=False)
    pairs = enumerate_all_pairs(inventory)
    pd.DataFrame([{k: v for k, v in p.items() if k not in {"early", "late"}} for p in pairs]).to_csv(
        outdir / "all_independent_epoch_pairs.csv", index=False
    )

    if not pairs:
        final = {
            "candidate_id": candidate_id, "ra_deg": ra_deg, "dec_deg": dec_deg,
            "pm_status": "INSUFFICIENT_DATA", "classification": "INSUFFICIENT_DATA",
            "reason": "NO_INDEPENDENT_JWST_EPOCH_PAIR", "n_independent_pairs": 0,
        }
        (outdir / "summary.json").write_text(json.dumps(final, indent=2))
        return final

    # Prepare each distinct physical epoch once.
    unique_epochs: dict[str, dict] = {}
    for pair in pairs:
        for side in ("early", "late"):
            epoch = pair[side]
            unique_epochs.setdefault(_epoch_key(epoch), epoch)

    cache: dict[str, dict] = {}
    for i, (key, epoch) in enumerate(sorted(unique_epochs.items(), key=lambda kv: kv[1]["mjd"]), 1):
        print(f"candidate={candidate_id} prepare epoch {i}/{len(unique_epochs)} {key}", flush=True)
        try:
            cache[key] = _prepare_epoch(candidate_id, ra_deg, dec_deg, epoch, outdir)
        except Exception as exc:
            # Preserve exhaustive pair accounting even when one epoch cannot be loaded.
            cache[key] = {
                "epoch_id": key, "filter": str(epoch["filter"]), "mjd": float(epoch["mjd"]),
                "n_covering": 0, "n_astrometric": 0,
                "measurements": pd.DataFrame(), "controls": pd.DataFrame(),
                "prepare_error": f"{type(exc).__name__}: {exc}",
            }

    epoch_qc = pd.DataFrame([{k: v for k, v in c.items() if k not in {"measurements", "controls"}} for c in cache.values()])
    epoch_qc.to_csv(outdir / "EPOCH_CACHE_QC.csv", index=False)

    pair_results = []
    for rank, pair in enumerate(pairs, 1):
        print(f"candidate={candidate_id} cached pair {rank}/{len(pairs)} {pair['pair_id']}", flush=True)
        try:
            r = _audit_pair_from_cache(candidate_id, pair, cache, outdir)
        except Exception as exc:
            r = {
                "candidate_id": candidate_id, "pair_id": pair["pair_id"],
                "pair_type": pair["pair_type"], "filter_early": pair["filter_early"],
                "filter_late": pair["filter_late"],
                "inventory_baseline_days": pair["baseline_days_inventory"],
                "pm_status": "INSUFFICIENT_DATA", "stage": "UNEXPECTED_ERROR",
                "reason": f"{type(exc).__name__}: {exc}",
            }
        r["attempt_rank"] = rank
        pair_results.append(r)

    pdf = pd.DataFrame(pair_results)
    pdf.to_csv(outdir / "ALL_PAIR_RESULTS.csv", index=False)
    measured = pdf[pdf["pm_status"].eq("PM_MEASURED")].copy()

    if measured.empty:
        final = {
            "candidate_id": candidate_id, "ra_deg": ra_deg, "dec_deg": dec_deg,
            "pm_status": "INSUFFICIENT_DATA", "classification": "INSUFFICIENT_DATA",
            "reason": "ALL_INDEPENDENT_JWST_PAIRS_EXHAUSTED",
            "n_independent_pairs": len(pairs), "n_measured_pairs": 0,
            "n_distinct_epochs_cached": len(unique_epochs),
        }
    else:
        measured["evidence_rank"] = measured["evidence_grade"].map({"A_SAME_FILTER_JWST": 0, "B_CROSS_FILTER_JWST": 1}).fillna(9)
        measured["max_pm_sigma"] = np.maximum(
            pd.to_numeric(measured["sigma_mu_alpha_cosdec_masyr"], errors="coerce"),
            pd.to_numeric(measured["sigma_mu_delta_masyr"], errors="coerce"),
        )
        best = measured.sort_values(["evidence_rank", "max_pm_sigma", "attempt_rank"]).iloc[0]
        final = {
            "candidate_id": candidate_id, "ra_deg": ra_deg, "dec_deg": dec_deg,
            "pm_status": "PM_MEASURED", "classification": best["classification"],
            "evidence_grade": best["evidence_grade"], "pair_id": best["pair_id"],
            "filter_early": best["filter_early"], "filter_late": best["filter_late"],
            "mu_alpha_cosdec_masyr": best["mu_alpha_cosdec_masyr"],
            "mu_delta_masyr": best["mu_delta_masyr"],
            "sigma_mu_alpha_cosdec_masyr": best["sigma_mu_alpha_cosdec_masyr"],
            "sigma_mu_delta_masyr": best["sigma_mu_delta_masyr"],
            "significance_2d": best.get("significance_2d", np.nan),
            "n_independent_pairs": len(pairs), "n_measured_pairs": len(measured),
            "n_distinct_epochs_cached": len(unique_epochs),
        }

    (outdir / "summary.json").write_text(json.dumps(final, indent=2, default=str))
    return final
