"""End-to-end independent candidate audit."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .archive import (
    choose_epoch_pair,
    get_pair_products,
    load_covering_cutouts,
    query_candidate_inventory,
)
from .astrometry import infer_proper_motion, register_exposures
from .config import CONFIG, FILTER_PRIORITY, FILTER_WAVELENGTH_NM
from .measurement import measure_exposure


INVENTORY_COLUMNS = [
    "candidate_id", "candidate_ra_deg", "candidate_dec_deg",
    "obs_collection", "obsid", "obs_id", "instrument_name", "filters",
    "filter_norm", "t_min", "t_max", "t_exptime", "proposal_id",
    "target_name", "calib_level", "dataRights", "distance",
]

MAX_PAIR_ATTEMPTS = 8


def _write_json(path: Path, obj: dict):
    clean = {}
    for key, value in obj.items():
        if isinstance(value, (np.integer,)):
            value = int(value)
        elif isinstance(value, (np.floating,)):
            value = float(value)
        if isinstance(value, float) and not np.isfinite(value):
            value = None
        clean[key] = value
    path.write_text(json.dumps(clean, indent=2, sort_keys=True), encoding="utf-8")


def _pair_metadata(pair: dict) -> dict:
    if pair.get("status") != "PAIR_FOUND":
        return {"pair_status": pair.get("status", "UNKNOWN")}
    return {
        "pair_status": "PAIR_FOUND",
        "pair_type": pair["pair_type"],
        "filter_early": pair["filter_early"],
        "filter_late": pair["filter_late"],
        "inventory_early_mjd": pair["early"]["mjd"],
        "inventory_late_mjd": pair["late"]["mjd"],
        "inventory_baseline_days": pair["baseline_days_inventory"],
        "pair_filter_cost": pair.get("pair_filter_cost"),
    }


def _group_epochs(rows: pd.DataFrame) -> list[dict]:
    rows = rows[np.isfinite(rows["t_min_numeric"])].sort_values("t_min_numeric")
    groups = []
    current = []
    previous = None
    for idx, row in rows.iterrows():
        mjd = float(row["t_min_numeric"])
        if previous is None or (mjd - previous) <= CONFIG.epoch_gap_days:
            current.append(idx)
        else:
            grp = rows.loc[current]
            groups.append({"mjd": float(np.nanmedian(grp["t_min_numeric"])), "rows": grp.copy()})
            current = [idx]
        previous = mjd
    if current:
        grp = rows.loc[current]
        groups.append({"mjd": float(np.nanmedian(grp["t_min_numeric"])), "rows": grp.copy()})
    return groups


def _pair_filter_cost(f1: str, f2: str) -> float:
    w1 = FILTER_WAVELENGTH_NM.get(f1, 0)
    w2 = FILTER_WAVELENGTH_NM.get(f2, 0)
    cost = abs(w1 - 4440) + abs(w2 - 4440)
    if f1 == f2:
        cost -= 1000.0
    return float(cost)


def _rank_epoch_pairs(inventory: pd.DataFrame, max_pairs: int = MAX_PAIR_ATTEMPTS) -> list[dict]:
    """Return several scientifically plausible NIRCam pairs in ranked order.

    The archive inventory is only an observation-level proximity query; it does not
    guarantee that the candidate falls on a detector in a particular visit.  A
    single preselected pair can therefore fail detector coverage or target S/N even
    when another legitimate pair exists.  We rank alternatives here and let the
    full gWCS/measurement chain decide which one is actually usable.
    """
    if inventory.empty:
        return []
    collection = inventory.get("obs_collection", pd.Series(index=inventory.index, dtype=object))
    jwst = inventory[collection.astype(str).str.upper().eq("JWST")].copy()
    if "instrument_name" in jwst:
        jwst = jwst[jwst["instrument_name"].astype(str).str.upper().str.contains("NIRCAM")]
    jwst = jwst[jwst["filter_norm"].notna()].copy()
    if jwst.empty:
        return []
    if "t_min_numeric" not in jwst:
        jwst["t_min_numeric"] = pd.to_numeric(jwst.get("t_min"), errors="coerce")

    all_epochs = []
    for filt, grp in jwst.groupby("filter_norm"):
        for epoch in _group_epochs(grp):
            all_epochs.append((str(filt), epoch))

    ranked = []
    for i, (f1, e1) in enumerate(all_epochs):
        for f2, e2 in all_epochs[i + 1:]:
            dt = abs(e2["mjd"] - e1["mjd"])
            if dt < CONFIG.min_pm_baseline_days:
                continue
            cost = _pair_filter_cost(f1, f2)
            same_penalty = 0 if f1 == f2 else 1
            p1 = FILTER_PRIORITY.index(f1) if f1 in FILTER_PRIORITY else 999
            p2 = FILTER_PRIORITY.index(f2) if f2 in FILTER_PRIORITY else 999
            ranked.append((cost, same_penalty, -dt, min(p1, p2), max(p1, p2), f1, e1, f2, e2))

    ranked.sort(key=lambda x: x[:5])
    out = []
    seen = set()
    for cost, _, _, _, _, f1, e1, f2, e2 in ranked:
        if e1["mjd"] <= e2["mjd"]:
            early_f, early, late_f, late = f1, e1, f2, e2
        else:
            early_f, early, late_f, late = f2, e2, f1, e1
        key = (early_f, round(early["mjd"], 5), late_f, round(late["mjd"], 5))
        if key in seen:
            continue
        seen.add(key)
        out.append({
            "status": "PAIR_FOUND",
            "pair_type": "SAME_FILTER_JWST" if early_f == late_f else "CROSS_FILTER_JWST",
            "filter_early": early_f,
            "filter_late": late_f,
            "early": early,
            "late": late,
            "baseline_days_inventory": float(late["mjd"] - early["mjd"]),
            "pair_filter_cost": float(cost),
        })
        if len(out) >= max_pairs:
            break
    return out


def _diagnostic_plot(candidate, cutouts, measurements, registered, outpath: Path):
    n = len(cutouts)
    if n == 0:
        return
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols)) + 1
    fig = plt.figure(figsize=(4 * ncols, 3.7 * nrows))
    gs = fig.add_gridspec(nrows, ncols)

    for i, exposure in enumerate(cutouts):
        ax = fig.add_subplot(gs[i // ncols, i % ncols])
        tx = exposure.x_full - exposure.x0
        ty = exposure.y_full - exposure.y0
        half = 15
        xc, yc = int(round(tx)), int(round(ty))
        xa, xb = max(0, xc-half), min(exposure.sci.shape[1], xc+half+1)
        ya, yb = max(0, yc-half), min(exposure.sci.shape[0], yc+half+1)
        stamp = np.asarray(exposure.sci[ya:yb, xa:xb], dtype=float)
        finite = stamp[np.isfinite(stamp)]
        if len(finite):
            vmin, vmax = np.nanpercentile(finite, [5, 99.5])
        else:
            vmin, vmax = 0, 1
        ax.imshow(stamp, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
        ax.axvline(tx-xa, lw=0.8)
        ax.axhline(ty-ya, lw=0.8)
        m = measurements[measurements["filename"].eq(exposure.filename)]
        if len(m):
            r = m.iloc[0]
            snr = r["raw_snr_emp"] if np.isfinite(r["raw_snr_emp"]) else r["raw_snr_err"]
            title = f"{exposure.epoch_label} {exposure.filter_name}\nS/N={snr:.1f}, clean={r['clean_snr_err']:.1f}"
        else:
            title = f"{exposure.epoch_label} {exposure.filter_name}"
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

    ax1 = fig.add_subplot(gs[-1, :max(1, ncols // 2)])
    if len(registered):
        for epoch, grp in registered.groupby("epoch_label"):
            ax1.scatter(grp["east_2dg_mas"], grp["north_2dg_mas"], label=epoch)
        ax1.legend(fontsize=8)
        ax1.set_xlabel("registered east [mas]")
        ax1.set_ylabel("registered north [mas]")
        ax1.grid(alpha=0.25)
    ax1.set_title("2DG local-registration positions")

    ax2 = fig.add_subplot(gs[-1, max(1, ncols // 2):])
    if len(measurements):
        snr = measurements["raw_snr_emp"].astype(float).copy()
        snr[~np.isfinite(snr)] = measurements.loc[~np.isfinite(snr), "raw_snr_err"]
        x = np.arange(len(measurements))
        ax2.plot(x, snr, marker="o")
        ax2.axhline(3.0, ls=":", lw=1)
        ax2.axhline(5.0, ls="--", lw=1)
        ax2.set_xticks(x)
        ax2.set_xticklabels(measurements["epoch_label"].astype(str), rotation=30)
        ax2.set_ylabel("forced-photometry S/N")
        ax2.grid(alpha=0.25)
    ax2.set_title("Target detections")

    fig.suptitle(f"Candidate {int(candidate['candidate_id'])}")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _copy_attempt_to_root(adir: Path, cdir: Path):
    for name in [
        "products_early.csv", "products_late.csv", "product_coverage_audit.csv",
        "exposure_measurements.csv", "field_controls.csv", "registered_positions.csv",
        "registration_quality.csv", "registration_control_residuals.csv", "pairwise_pm.csv",
        "diagnostic.png",
    ]:
        src = adir / name
        if src.exists():
            shutil.copy2(src, cdir / name)


def run_candidate(candidate: dict, output_root: Path) -> dict:
    cid = int(candidate["candidate_id"])
    ra = float(candidate["ra_deg"])
    dec = float(candidate["dec_deg"])
    cdir = output_root / f"candidate_{cid}"
    cdir.mkdir(parents=True, exist_ok=True)

    base = {
        "candidate_id": cid,
        "ra_deg": ra,
        "dec_deg": dec,
        "input_f444w_mag_ab": candidate.get("f444w_mag_ab"),
        "input_f444w_flux_njy": candidate.get("f444w_flux_njy"),
    }
    for key, value in candidate.items():
        if key not in {"candidate_id", "ra_deg", "dec_deg", "f444w_mag_ab", "f444w_flux_njy"}:
            base[f"catalog_{key}"] = value

    try:
        inventory = query_candidate_inventory(cid, ra, dec)
        slim_cols = [c for c in INVENTORY_COLUMNS if c in inventory.columns]
        inventory[slim_cols].to_csv(cdir / "archive_inventory.csv", index=False)

        n_hst = int((inventory.get("obs_collection", pd.Series(dtype=str)).astype(str).str.upper() == "HST").sum())
        n_jwst = int((inventory.get("obs_collection", pd.Series(dtype=str)).astype(str).str.upper() == "JWST").sum())

        ranked_pairs = _rank_epoch_pairs(inventory)
        if not ranked_pairs:
            pair = choose_epoch_pair(inventory)
            pair_meta = _pair_metadata(pair)
            _write_json(cdir / "selected_epoch_pair.json", pair_meta)
            summary = {
                **base, **pair_meta,
                "n_jwst_inventory_rows": n_jwst,
                "n_hst_inventory_rows": n_hst,
                "pm_status": "INSUFFICIENT_DATA",
                "classification": "INSUFFICIENT_DATA",
                "reason": pair.get("status", "no usable epoch pair"),
                "n_pair_attempts": 0,
            }
            _write_json(cdir / "summary.json", summary)
            return summary

        attempts_dir = cdir / "pair_attempts"
        attempts_dir.mkdir(exist_ok=True)
        attempt_rows = []
        best_failure = None

        for rank, pair in enumerate(ranked_pairs, start=1):
            adir = attempts_dir / f"attempt_{rank:02d}"
            adir.mkdir(exist_ok=True)
            pair_meta = _pair_metadata(pair)
            _write_json(adir / "pair.json", pair_meta)

            try:
                early_products, late_products = get_pair_products(pair)
                early_products.to_csv(adir / "products_early.csv", index=False)
                late_products.to_csv(adir / "products_late.csv", index=False)

                early_cutouts, audit_early = load_covering_cutouts(
                    cid, ra, dec, "early", pair["filter_early"], early_products
                )
                late_cutouts, audit_late = load_covering_cutouts(
                    cid, ra, dec, "late", pair["filter_late"], late_products
                )
                product_audit = pd.DataFrame(audit_early + audit_late)
                product_audit.to_csv(adir / "product_coverage_audit.csv", index=False)
                cutouts = early_cutouts + late_cutouts

                if not early_cutouts or not late_cutouts:
                    reason = "selected archive epochs do not both have covering CAL exposures"
                    stage = 1
                    attempt_rows.append({**pair_meta, "attempt_rank": rank, "stage": stage, "reason": reason,
                                         "n_covering_early": len(early_cutouts), "n_covering_late": len(late_cutouts)})
                    candidate_failure = (stage, -rank, adir, pair_meta, reason, len(early_cutouts), len(late_cutouts), 0, 0)
                    if best_failure is None or candidate_failure[:2] > best_failure[:2]:
                        best_failure = candidate_failure
                    continue

                measurement_rows = []
                control_frames = []
                for exposure in cutouts:
                    row, controls = measure_exposure(exposure)
                    measurement_rows.append(row)
                    if len(controls):
                        control_frames.append(controls)

                measurements = pd.DataFrame(measurement_rows).sort_values("mjd").reset_index(drop=True)
                controls = pd.concat(control_frames, ignore_index=True) if control_frames else pd.DataFrame()
                measurements.to_csv(adir / "exposure_measurements.csv", index=False)
                controls.to_csv(adir / "field_controls.csv", index=False)

                usable_by_epoch = measurements[
                    (measurements["clean_snr_err"] >= 3)
                    & np.isfinite(measurements["east_2dg_arcsec_raw_wcs"])
                ].groupby("epoch_label").size()
                n_early = int(usable_by_epoch.get("early", 0))
                n_late = int(usable_by_epoch.get("late", 0))

                if n_early == 0 or n_late == 0:
                    reason = "target is not securely centroidable in both independent epochs"
                    stage = 2
                    _diagnostic_plot(candidate, cutouts, measurements, pd.DataFrame(), adir / "diagnostic.png")
                    attempt_rows.append({**pair_meta, "attempt_rank": rank, "stage": stage, "reason": reason,
                                         "n_covering_early": len(early_cutouts), "n_covering_late": len(late_cutouts),
                                         "n_astrometric_early": n_early, "n_astrometric_late": n_late})
                    candidate_failure = (stage, -rank, adir, pair_meta, reason, len(early_cutouts), len(late_cutouts), n_early, n_late)
                    if best_failure is None or candidate_failure[:2] > best_failure[:2]:
                        best_failure = candidate_failure
                    continue

                registered, quality, matched = register_exposures(measurements, controls)
                registered.to_csv(adir / "registered_positions.csv", index=False)
                quality.to_csv(adir / "registration_quality.csv", index=False)
                matched.to_csv(adir / "registration_control_residuals.csv", index=False)

                pm, pairwise = infer_proper_motion(registered, quality, pair["pair_type"])
                pairwise.to_csv(adir / "pairwise_pm.csv", index=False)
                _diagnostic_plot(candidate, cutouts, measurements, registered, adir / "diagnostic.png")

                classification = pm.get("classification") or "INSUFFICIENT_DATA"
                if classification == "INSUFFICIENT_DATA":
                    reason = pm.get("reason", "local registration/proper-motion inference insufficient")
                    stage = 3
                    attempt_rows.append({**pair_meta, "attempt_rank": rank, "stage": stage, "reason": reason,
                                         "n_covering_early": len(early_cutouts), "n_covering_late": len(late_cutouts),
                                         "n_astrometric_early": n_early, "n_astrometric_late": n_late})
                    candidate_failure = (stage, -rank, adir, pair_meta, reason, len(early_cutouts), len(late_cutouts), n_early, n_late)
                    if best_failure is None or candidate_failure[:2] > best_failure[:2]:
                        best_failure = candidate_failure
                    continue

                attempt_rows.append({**pair_meta, "attempt_rank": rank, "stage": 4, "reason": "PM_MEASURED",
                                     "n_covering_early": len(early_cutouts), "n_covering_late": len(late_cutouts),
                                     "n_astrometric_early": n_early, "n_astrometric_late": n_late})
                pd.DataFrame(attempt_rows).to_csv(cdir / "pair_attempts.csv", index=False)
                _write_json(cdir / "selected_epoch_pair.json", {**pair_meta, "attempt_rank": rank})
                _copy_attempt_to_root(adir, cdir)
                summary = {
                    **base, **pair_meta,
                    "n_jwst_inventory_rows": n_jwst,
                    "n_hst_inventory_rows": n_hst,
                    "n_covering_early": len(early_cutouts),
                    "n_covering_late": len(late_cutouts),
                    "n_astrometric_early": n_early,
                    "n_astrometric_late": n_late,
                    "n_pair_attempts": rank,
                    "selected_pair_rank": rank,
                    **pm,
                }
                _write_json(cdir / "summary.json", summary)
                return summary

            except Exception as exc:
                reason = f"PAIR_ATTEMPT_ERROR:{type(exc).__name__}:{exc}"
                attempt_rows.append({**pair_meta, "attempt_rank": rank, "stage": -1, "reason": reason})
                continue

        pd.DataFrame(attempt_rows).to_csv(cdir / "pair_attempts.csv", index=False)
        if best_failure is not None:
            stage, neg_rank, adir, pair_meta, reason, ncov_e, ncov_l, nast_e, nast_l = best_failure
            rank = -neg_rank
            _write_json(cdir / "selected_epoch_pair.json", {**pair_meta, "attempt_rank": rank, "selection_note": "best failed pair after alternate-pair search"})
            _copy_attempt_to_root(adir, cdir)
            summary = {
                **base, **pair_meta,
                "n_jwst_inventory_rows": n_jwst,
                "n_hst_inventory_rows": n_hst,
                "n_covering_early": ncov_e,
                "n_covering_late": ncov_l,
                "n_astrometric_early": nast_e,
                "n_astrometric_late": nast_l,
                "n_pair_attempts": len(ranked_pairs),
                "selected_pair_rank": rank,
                "pm_status": "INSUFFICIENT_DATA",
                "classification": "INSUFFICIENT_DATA",
                "reason": f"all {len(ranked_pairs)} ranked NIRCam pairs exhausted; best failure: {reason}",
            }
        else:
            summary = {
                **base,
                "n_jwst_inventory_rows": n_jwst,
                "n_hst_inventory_rows": n_hst,
                "n_pair_attempts": len(ranked_pairs),
                "pm_status": "ERROR",
                "classification": "INSUFFICIENT_DATA",
                "reason": "all ranked pair attempts failed with execution errors; inspect pair_attempts.csv",
            }
        _write_json(cdir / "summary.json", summary)
        return summary

    except Exception as exc:
        summary = {
            **base,
            "pm_status": "ERROR",
            "classification": "INSUFFICIENT_DATA",
            "reason": f"{type(exc).__name__}: {exc}",
        }
        _write_json(cdir / "summary.json", summary)
        return summary
