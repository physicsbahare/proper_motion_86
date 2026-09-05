"""Exhaustive, conservative JWST/NIRCam proper-motion audit.

This module reuses the validated pm86 measurement/registration machinery but changes
candidate handling in two important ways:

1. all independent JWST/NIRCam epoch pairs are enumerated and attempted;
2. a DQ-aware hard gate rejects target centroids landing on DO_NOT_USE, SATURATED,
   JUMP_DET or OUTLIER pixels before a PM fit is allowed.

The goal is scientific defensibility rather than speed.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd

import pm86.archive as archive
import pm86.measurement as measurement
import pm86.astrometry as astrometry
from pm86.config import CONFIG, FILTER_PRIORITY

from first26.forced_astrometry import fit_forced_gaussian
from first26.recenter_measurement import _nearest_seed

ASTROM_BAD_BITS = int(CONFIG.astrom_bad_bits)


def enumerate_all_pairs(inventory: pd.DataFrame) -> list[dict]:
    """Return every legal independent JWST/NIRCam epoch pair.

    The existing pipeline ranks useful pairs; here ranking controls execution order
    only. It never truncates the search, so INSUFFICIENT_DATA means all legal pairs
    were actually exhausted.
    """
    if inventory.empty:
        return []
    collection = inventory.get("obs_collection", pd.Series(index=inventory.index, dtype=object))
    jwst = inventory[collection.astype(str).str.upper().eq("JWST")].copy()
    if "instrument_name" in jwst:
        jwst = jwst[jwst["instrument_name"].astype(str).str.upper().str.contains("NIRCAM", na=False)]
    jwst = jwst[jwst["filter_norm"].notna()].copy()
    if jwst.empty:
        return []

    epochs = []
    for filt, grp in jwst.groupby("filter_norm"):
        for i, epoch in enumerate(archive._group_epochs(grp)):
            epochs.append({
                "epoch_id": f"{filt}_epoch_{i:03d}",
                "filter": str(filt),
                "mjd": float(epoch["mjd"]),
                "rows": epoch["rows"].copy(),
            })
    epochs.sort(key=lambda e: e["mjd"])

    pairs = []
    for i in range(len(epochs)):
        for j in range(i + 1, len(epochs)):
            a, b = epochs[i], epochs[j]
            dt = float(b["mjd"] - a["mjd"])
            if dt < CONFIG.min_pm_baseline_days:
                continue
            f1, f2 = a["filter"], b["filter"]
            cost = archive._pair_filter_cost(f1, f2)
            p1 = FILTER_PRIORITY.index(f1) if f1 in FILTER_PRIORITY else 999
            p2 = FILTER_PRIORITY.index(f2) if f2 in FILTER_PRIORITY else 999
            pairs.append({
                "pair_id": f"{a['epoch_id']}__{b['epoch_id']}",
                "status": "PAIR_FOUND",
                "pair_type": "SAME_FILTER_JWST" if f1 == f2 else "CROSS_FILTER_JWST",
                "filter_early": f1,
                "filter_late": f2,
                "early": a,
                "late": b,
                "baseline_days_inventory": dt,
                "pair_filter_cost": float(cost),
                "_sort": (0 if f1 == f2 else 1, cost, -dt, min(p1, p2), max(p1, p2)),
            })
    pairs.sort(key=lambda p: p["_sort"])
    for p in pairs:
        p.pop("_sort", None)
    return pairs


def _centroid_pixel_clean(exposure, x, y) -> bool:
    if not (np.isfinite(x) and np.isfinite(y)):
        return False
    dq = np.asarray(exposure.dq, dtype=np.uint64)
    xi, yi = int(round(x)), int(round(y))
    if not (0 <= yi < dq.shape[0] and 0 <= xi < dq.shape[1]):
        return False
    return (int(dq[yi, xi]) & ASTROM_BAD_BITS) == 0


def measure_exposure_dq_aware(exposure):
    """Validated recenter/forced-centroid logic plus a hard DQ gate at the adopted centroid."""
    sci = np.asarray(exposure.sci, dtype=float)
    err = np.asarray(exposure.err, dtype=float)
    dq = np.asarray(exposure.dq, dtype=np.uint32)
    phot_bad, astrom_bad = measurement.masks(sci, err, dq)

    tx = exposure.x_full - exposure.x0
    ty = exposure.y_full - exposure.y0
    _, bkg, _ = measurement.sigma_clipped_stats(sci, mask=phot_bad, sigma=3.0, maxiters=5)
    data_sub = sci - float(bkg)

    # First run the base measurement to preserve all diagnostics and field controls.
    row, controls = measurement.measure_exposure(exposure)
    base_snr = float(row.get("clean_snr_err", np.nan))
    row["aperture_clean_snr_err"] = base_snr
    row["astrometric_snr"] = np.nan
    row["astrometry_source"] = "UNUSABLE"
    row["target_seed_recentered"] = False
    row["target_seed_offset_pix"] = 0.0
    row["forced_fit_attempted"] = False
    row["forced_fit_accepted"] = False

    # Direct catalog-position solution: S/N + finite centroid + clean adopted pixel.
    if (
        np.isfinite(base_snr) and base_snr >= CONFIG.target_min_snr_astrometry
        and np.isfinite(row.get("x_2dg", np.nan))
        and np.isfinite(row.get("y_2dg", np.nan))
        and _centroid_pixel_clean(exposure, float(row["x_2dg"]), float(row["y_2dg"]))
    ):
        row["astrometric_snr"] = base_snr
        row["astrometry_source"] = "catalog_position_2dg"
        row["astrometrically_usable"] = True
        return row, controls

    # Conservative local association seed, bounded to 4 pixels by the validated helper.
    sx, sy, offset, found = _nearest_seed(data_sub, astrom_bad, tx, ty)
    seed_x, seed_y = (sx, sy) if found else (tx, ty)
    row["target_seed_recentered"] = bool(found)
    row["target_seed_offset_pix"] = float(offset)

    forced = fit_forced_gaussian(data_sub, err, astrom_bad, seed_x, seed_y)
    row.update(forced)
    if forced.get("forced_fit_accepted", False):
        fx, fy = float(forced["forced_x"]), float(forced["forced_y"])
        total_offset = float(np.hypot(fx - tx, fy - ty))
        if total_offset <= 4.0 and _centroid_pixel_clean(exposure, fx, fy):
            cen = measurement.centroid_stamp(data_sub, err, astrom_bad, fx, fy, n_mc=0)
            e2, n2 = _to_tangent(exposure, fx, fy)
            ec, nc = _to_tangent(exposure, cen["x_com"], cen["y_com"])
            row.update({
                "astrometric_snr": float(forced["forced_snr"]),
                "astrometry_source": "forced_gaussian",
                "target_seed_offset_pix": total_offset,
                "x_2dg": fx, "y_2dg": fy,
                "x_com": cen["x_com"], "y_com": cen["y_com"],
                "sigma_x_pix": float(forced["forced_x_err_pix"]),
                "sigma_y_pix": float(forced["forced_y_err_pix"]),
                "east_2dg_arcsec_raw_wcs": e2, "north_2dg_arcsec_raw_wcs": n2,
                "east_com_arcsec_raw_wcs": ec, "north_com_arcsec_raw_wcs": nc,
                "clean_snr_err": float(forced["forced_snr"]),
                "astrometrically_usable": True,
            })
            return row, controls

    # Final segmentation recenter fallback.
    if found:
        _, _, csnr, _ = measurement.aperture_measure(sci, err, astrom_bad, sx, sy)
        cen = measurement.centroid_stamp(
            data_sub, err, astrom_bad, sx, sy,
            n_mc=CONFIG.centroid_mc_draws,
            seed=measurement.stable_seed(f"dqrecenter|{exposure.candidate_id}|{exposure.filename}"),
        )
        if (
            np.isfinite(csnr) and csnr >= CONFIG.target_min_snr_astrometry
            and np.isfinite(cen["x_2dg"]) and np.isfinite(cen["y_2dg"])
            and float(np.hypot(cen["x_2dg"] - tx, cen["y_2dg"] - ty)) <= 4.0
            and _centroid_pixel_clean(exposure, float(cen["x_2dg"]), float(cen["y_2dg"]))
        ):
            e2, n2 = _to_tangent(exposure, cen["x_2dg"], cen["y_2dg"])
            ec, nc = _to_tangent(exposure, cen["x_com"], cen["y_com"])
            row.update({
                "astrometric_snr": float(csnr),
                "astrometry_source": "segmentation_recenter_2dg",
                "target_seed_offset_pix": float(np.hypot(cen["x_2dg"] - tx, cen["y_2dg"] - ty)),
                **cen,
                "east_2dg_arcsec_raw_wcs": e2, "north_2dg_arcsec_raw_wcs": n2,
                "east_com_arcsec_raw_wcs": ec, "north_com_arcsec_raw_wcs": nc,
                "clean_snr_err": float(csnr),
                "astrometrically_usable": True,
            })
            return row, controls

    row["astrometrically_usable"] = False
    return row, controls


def _to_tangent(exposure, x, y):
    if not (np.isfinite(x) and np.isfinite(y)):
        return np.nan, np.nan
    try:
        ra, dec = measurement.local_pixel_to_sky(exposure, x, y)
        e, n = measurement.sky_to_tangent(ra, dec, exposure.ra_deg, exposure.dec_deg)
        return float(e), float(n)
    except Exception:
        return np.nan, np.nan


def _load_all_covering(candidate_id, ra_deg, dec_deg, epoch_label, filt, products):
    """Use the existing CAL loader but remove the previous four-exposure cap."""
    old = CONFIG.max_covering_exposures_per_epoch
    # Config is frozen, so temporarily replace the module-level CONFIG object.
    cfg = copy.copy(CONFIG)
    object.__setattr__(cfg, "max_covering_exposures_per_epoch", 10_000)
    old_archive_cfg = archive.CONFIG
    archive.CONFIG = cfg
    try:
        return archive.load_covering_cutouts(candidate_id, ra_deg, dec_deg, epoch_label, filt, products)
    finally:
        archive.CONFIG = old_archive_cfg


def audit_pair(candidate_id: int, ra_deg: float, dec_deg: float, pair: dict, outdir: Path) -> dict:
    pdir = outdir / "pairs" / pair["pair_id"]
    pdir.mkdir(parents=True, exist_ok=True)

    pe = archive._get_products_with_retry(pair["early"]["rows"])
    pl = archive._get_products_with_retry(pair["late"]["rows"])
    pe.to_csv(pdir / "products_early.csv", index=False)
    pl.to_csv(pdir / "products_late.csv", index=False)

    ce, ae = _load_all_covering(candidate_id, ra_deg, dec_deg, "early", pair["filter_early"], pe)
    cl, al = _load_all_covering(candidate_id, ra_deg, dec_deg, "late", pair["filter_late"], pl)
    pd.DataFrame(ae + al).to_csv(pdir / "product_coverage_audit.csv", index=False)

    result = {
        "candidate_id": candidate_id,
        "pair_id": pair["pair_id"],
        "pair_type": pair["pair_type"],
        "filter_early": pair["filter_early"],
        "filter_late": pair["filter_late"],
        "inventory_baseline_days": pair["baseline_days_inventory"],
        "n_covering_early": len(ce),
        "n_covering_late": len(cl),
    }
    if not ce or not cl:
        result.update(pm_status="INSUFFICIENT_DATA", stage="COVERAGE_FAIL", reason="epochs do not both have covering CAL exposures")
        return result

    rows, ctrs = [], []
    for exp in [*ce, *cl]:
        r, c = measure_exposure_dq_aware(exp)
        rows.append(r)
        if len(c):
            ctrs.append(c)
    meas = pd.DataFrame(rows)
    controls = pd.concat(ctrs, ignore_index=True) if ctrs else pd.DataFrame()
    meas.to_csv(pdir / "exposure_measurements.csv", index=False)
    controls.to_csv(pdir / "field_controls.csv", index=False)

    early_ok = int(((meas["epoch_label"] == "early") & meas["astrometrically_usable"].astype(bool)).sum())
    late_ok = int(((meas["epoch_label"] == "late") & meas["astrometrically_usable"].astype(bool)).sum())
    result.update(n_astrometric_early=early_ok, n_astrometric_late=late_ok)
    if early_ok == 0 or late_ok == 0:
        result.update(pm_status="INSUFFICIENT_DATA", stage="CENTROID_FAIL", reason="target is not securely centroidable in both independent epochs")
        return result

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


def run_candidate(candidate_id: int, ra_deg: float, dec_deg: float, output_root: str | Path) -> dict:
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
            "candidate_id": candidate_id,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "pm_status": "INSUFFICIENT_DATA",
            "classification": "INSUFFICIENT_DATA",
            "reason": "NO_INDEPENDENT_JWST_EPOCH_PAIR",
            "n_independent_pairs": 0,
        }
        (outdir / "summary.json").write_text(json.dumps(final, indent=2))
        return final

    pair_results = []
    for rank, pair in enumerate(pairs, 1):
        print(f"candidate={candidate_id} pair {rank}/{len(pairs)} {pair['pair_id']}", flush=True)
        try:
            r = audit_pair(candidate_id, ra_deg, dec_deg, pair, outdir)
        except Exception as exc:
            r = {
                "candidate_id": candidate_id,
                "pair_id": pair["pair_id"],
                "pair_type": pair["pair_type"],
                "filter_early": pair["filter_early"],
                "filter_late": pair["filter_late"],
                "inventory_baseline_days": pair["baseline_days_inventory"],
                "pm_status": "INSUFFICIENT_DATA",
                "stage": "UNEXPECTED_ERROR",
                "reason": f"{type(exc).__name__}: {exc}",
            }
        r["attempt_rank"] = rank
        pair_results.append(r)

    pdf = pd.DataFrame(pair_results)
    pdf.to_csv(outdir / "ALL_PAIR_RESULTS.csv", index=False)
    measured = pdf[pdf["pm_status"].eq("PM_MEASURED")].copy()

    if measured.empty:
        final = {
            "candidate_id": candidate_id,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "pm_status": "INSUFFICIENT_DATA",
            "classification": "INSUFFICIENT_DATA",
            "reason": "ALL_INDEPENDENT_JWST_PAIRS_EXHAUSTED",
            "n_independent_pairs": len(pairs),
            "n_measured_pairs": 0,
        }
    else:
        measured["evidence_rank"] = measured["evidence_grade"].map({"A_SAME_FILTER_JWST": 0, "B_CROSS_FILTER_JWST": 1}).fillna(9)
        measured["max_pm_sigma"] = np.maximum(
            pd.to_numeric(measured["sigma_mu_alpha_cosdec_masyr"], errors="coerce"),
            pd.to_numeric(measured["sigma_mu_delta_masyr"], errors="coerce"),
        )
        best = measured.sort_values(["evidence_rank", "max_pm_sigma", "attempt_rank"]).iloc[0]
        final = {
            "candidate_id": candidate_id,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "pm_status": "PM_MEASURED",
            "classification": best["classification"],
            "evidence_grade": best["evidence_grade"],
            "pair_id": best["pair_id"],
            "filter_early": best["filter_early"],
            "filter_late": best["filter_late"],
            "mu_alpha_cosdec_masyr": best["mu_alpha_cosdec_masyr"],
            "mu_delta_masyr": best["mu_delta_masyr"],
            "sigma_mu_alpha_cosdec_masyr": best["sigma_mu_alpha_cosdec_masyr"],
            "sigma_mu_delta_masyr": best["sigma_mu_delta_masyr"],
            "significance_2d": best["significance_2d"],
            "n_independent_pairs": len(pairs),
            "n_measured_pairs": len(measured),
        }
    (outdir / "summary.json").write_text(json.dumps(final, indent=2, default=str))
    return final
