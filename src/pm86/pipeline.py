"""End-to-end independent candidate audit."""

from __future__ import annotations

import json
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
from .measurement import measure_exposure


INVENTORY_COLUMNS = [
    "candidate_id", "candidate_ra_deg", "candidate_dec_deg",
    "obs_collection", "obsid", "obs_id", "instrument_name", "filters",
    "filter_norm", "t_min", "t_max", "t_exptime", "proposal_id",
    "target_name", "calib_level", "dataRights", "distance",
]


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
    }


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

    try:
        inventory = query_candidate_inventory(cid, ra, dec)
        slim_cols = [c for c in INVENTORY_COLUMNS if c in inventory.columns]
        inventory[slim_cols].to_csv(cdir / "archive_inventory.csv", index=False)

        n_hst = int((inventory.get("obs_collection", pd.Series(dtype=str)).astype(str).str.upper() == "HST").sum())
        n_jwst = int((inventory.get("obs_collection", pd.Series(dtype=str)).astype(str).str.upper() == "JWST").sum())

        pair = choose_epoch_pair(inventory)
        pair_meta = _pair_metadata(pair)
        _write_json(cdir / "selected_epoch_pair.json", pair_meta)

        if pair.get("status") != "PAIR_FOUND":
            summary = {
                **base,
                **pair_meta,
                "n_jwst_inventory_rows": n_jwst,
                "n_hst_inventory_rows": n_hst,
                "pm_status": "INSUFFICIENT_DATA",
                "classification": "INSUFFICIENT_DATA",
                "reason": pair.get("status", "no usable epoch pair"),
            }
            _write_json(cdir / "summary.json", summary)
            return summary

        early_products, late_products = get_pair_products(pair)
        early_products.to_csv(cdir / "products_early.csv", index=False)
        late_products.to_csv(cdir / "products_late.csv", index=False)

        early_cutouts, audit_early = load_covering_cutouts(
            cid, ra, dec, "early", pair["filter_early"], early_products
        )
        late_cutouts, audit_late = load_covering_cutouts(
            cid, ra, dec, "late", pair["filter_late"], late_products
        )
        product_audit = pd.DataFrame(audit_early + audit_late)
        product_audit.to_csv(cdir / "product_coverage_audit.csv", index=False)
        cutouts = early_cutouts + late_cutouts

        if not early_cutouts or not late_cutouts:
            summary = {
                **base,
                **pair_meta,
                "n_jwst_inventory_rows": n_jwst,
                "n_hst_inventory_rows": n_hst,
                "n_covering_early": len(early_cutouts),
                "n_covering_late": len(late_cutouts),
                "pm_status": "INSUFFICIENT_DATA",
                "classification": "INSUFFICIENT_DATA",
                "reason": "selected archive epochs do not both have covering CAL exposures",
            }
            _write_json(cdir / "summary.json", summary)
            return summary

        measurement_rows = []
        control_frames = []
        for exposure in cutouts:
            row, controls = measure_exposure(exposure)
            measurement_rows.append(row)
            if len(controls):
                control_frames.append(controls)

        measurements = pd.DataFrame(measurement_rows).sort_values("mjd").reset_index(drop=True)
        controls = pd.concat(control_frames, ignore_index=True) if control_frames else pd.DataFrame()
        measurements.to_csv(cdir / "exposure_measurements.csv", index=False)
        controls.to_csv(cdir / "field_controls.csv", index=False)

        # A cross-filter non-detection is not evidence of motion. Require a
        # clean target centroid in both independent epochs before registration.
        usable_by_epoch = measurements[
            (measurements["clean_snr_err"] >= 3)
            & np.isfinite(measurements["east_2dg_arcsec_raw_wcs"])
        ].groupby("epoch_label").size()

        if usable_by_epoch.get("early", 0) == 0 or usable_by_epoch.get("late", 0) == 0:
            summary = {
                **base,
                **pair_meta,
                "n_jwst_inventory_rows": n_jwst,
                "n_hst_inventory_rows": n_hst,
                "n_covering_early": len(early_cutouts),
                "n_covering_late": len(late_cutouts),
                "n_astrometric_early": int(usable_by_epoch.get("early", 0)),
                "n_astrometric_late": int(usable_by_epoch.get("late", 0)),
                "pm_status": "INSUFFICIENT_DATA",
                "classification": "INSUFFICIENT_DATA",
                "reason": "target is not securely centroidable in both independent epochs",
            }
            _diagnostic_plot(candidate, cutouts, measurements, pd.DataFrame(), cdir / "diagnostic.png")
            _write_json(cdir / "summary.json", summary)
            return summary

        registered, quality, matched = register_exposures(measurements, controls)
        registered.to_csv(cdir / "registered_positions.csv", index=False)
        quality.to_csv(cdir / "registration_quality.csv", index=False)
        matched.to_csv(cdir / "registration_control_residuals.csv", index=False)

        pm, pairwise = infer_proper_motion(registered, quality, pair["pair_type"])
        pairwise.to_csv(cdir / "pairwise_pm.csv", index=False)
        _diagnostic_plot(candidate, cutouts, measurements, registered, cdir / "diagnostic.png")

        summary = {
            **base,
            **pair_meta,
            "n_jwst_inventory_rows": n_jwst,
            "n_hst_inventory_rows": n_hst,
            "n_covering_early": len(early_cutouts),
            "n_covering_late": len(late_cutouts),
            "n_astrometric_early": int(usable_by_epoch.get("early", 0)),
            "n_astrometric_late": int(usable_by_epoch.get("late", 0)),
            **pm,
        }
        if summary.get("classification") is None:
            summary["classification"] = "INSUFFICIENT_DATA"
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
