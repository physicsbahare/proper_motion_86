#!/usr/bin/env python3
"""Final threshold/centroid/pair-selection robustness audit for first26.

This is a conservative *screening* audit over the saved all-pair evidence.  It asks
whether any candidate currently classified INSUFFICIENT_DATA could acquire a usable
measurement in both epochs under reasonable perturbations of the astrometric gates.
It never allows >4 pixel source association (the 643980 neighbour-latching failure
mode), and never treats DO_NOT_USE or SATURATED pixels as acceptable.

If a candidate reaches both epochs in any sensitivity scenario it is flagged for a
fresh targeted PM rerun; otherwise its current INSUFFICIENT_DATA classification is
robust to the tested gate choices.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

SCENARIOS = {
    # Current-like gate reconstructed from saved diagnostics.
    "current_like": dict(snr=4.0, max_poserr=1.2, max_sigma=1.8, max_axis=2.0, max_rchi2=3.0),
    # One-at-a-time perturbations.
    "snr_3": dict(snr=3.0, max_poserr=1.2, max_sigma=1.8, max_axis=2.0, max_rchi2=3.0),
    "poserr_1p5": dict(snr=4.0, max_poserr=1.5, max_sigma=1.8, max_axis=2.0, max_rchi2=3.0),
    "sigma_2p2": dict(snr=4.0, max_poserr=1.2, max_sigma=2.2, max_axis=2.0, max_rchi2=3.0),
    "axis_2p5": dict(snr=4.0, max_poserr=1.2, max_sigma=1.8, max_axis=2.5, max_rchi2=3.0),
    "rchi2_5": dict(snr=4.0, max_poserr=1.2, max_sigma=1.8, max_axis=2.0, max_rchi2=5.0),
    # Deliberately generous combined sensitivity test, still constrained to a local
    # <=4-pixel association and clean/non-saturated pixels.
    "combined_relaxed": dict(snr=3.0, max_poserr=1.5, max_sigma=2.2, max_axis=2.5, max_rchi2=5.0),
}

BAD_HARD = ("dq_do_not_use", "dq_saturated")


def truthy(v):
    if isinstance(v, bool):
        return v
    if pd.isna(v):
        return False
    return str(v).strip().lower() in {"true", "1", "yes"}


def f(row, key):
    try:
        x = float(row.get(key, np.nan))
        return x
    except Exception:
        return np.nan


def local_association_ok(row):
    off = f(row, "target_seed_offset_pix")
    # Catalog-position rows generally carry 0/NaN; forced/recentered solutions must
    # remain inside the scientifically defended 4-pixel envelope.
    return (not np.isfinite(off)) or off <= 4.0 + 1e-9


def hard_dq_ok(row):
    return not any(truthy(row.get(k, False)) for k in BAD_HARD)


def ordinary_centroid_ok(row, snr_threshold=3.0, method="2dg"):
    if not local_association_ok(row) or not hard_dq_ok(row):
        return False
    snr = f(row, "aperture_clean_snr_err")
    if not np.isfinite(snr):
        snr = f(row, "clean_snr_err")
    if not np.isfinite(snr) or snr < snr_threshold:
        return False
    if method == "2dg":
        return np.isfinite(f(row, "x_2dg")) and np.isfinite(f(row, "y_2dg"))
    return np.isfinite(f(row, "x_com")) and np.isfinite(f(row, "y_com"))


def forced_ok(row, cfg):
    if not local_association_ok(row) or not hard_dq_ok(row):
        return False
    if not truthy(row.get("forced_fit_attempted", False)):
        return False
    snr = f(row, "forced_snr")
    sx, sy = f(row, "forced_sigma_x_pix"), f(row, "forced_sigma_y_pix")
    ar = f(row, "forced_axis_ratio")
    xe, ye = f(row, "forced_x_err_pix"), f(row, "forced_y_err_pix")
    rc = f(row, "forced_reduced_chi2")
    x, y = f(row, "forced_x"), f(row, "forced_y")
    vals = [snr, sx, sy, ar, xe, ye, rc, x, y]
    if not all(np.isfinite(v) for v in vals):
        return False
    return bool(
        snr >= cfg["snr"]
        and xe <= cfg["max_poserr"] and ye <= cfg["max_poserr"]
        and 0.55 < sx < cfg["max_sigma"] and 0.55 < sy < cfg["max_sigma"]
        and ar <= cfg["max_axis"] and rc <= cfg["max_rchi2"]
    )


def exposure_ok(row, scenario):
    cfg = SCENARIOS[scenario]
    # Sensitivity to centroid method: either ordinary 2DG or COM is enough for the
    # screen; a forced fit is evaluated with the scenario-specific quality gates.
    return (
        ordinary_centroid_ok(row, snr_threshold=max(3.0, cfg["snr"]), method="2dg")
        or ordinary_centroid_ok(row, snr_threshold=max(3.0, cfg["snr"]), method="com")
        or forced_ok(row, cfg)
    )


def audit_attempt(csv_path: Path):
    df = pd.read_csv(csv_path)
    out = {}
    for scenario in SCENARIOS:
        usable = df.apply(lambda r: exposure_ok(r, scenario), axis=1)
        tmp = df.loc[usable].copy()
        ne = int((tmp.get("epoch_label", pd.Series(dtype=str)).astype(str) == "early").sum())
        nl = int((tmp.get("epoch_label", pd.Series(dtype=str)).astype(str) == "late").sum())
        out[scenario] = {"n_early": ne, "n_late": nl, "two_epoch": ne > 0 and nl > 0}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_first26")
    ap.add_argument("--output", default="results_first26/final_robustness_audit")
    args = ap.parse_args()
    root = Path(args.results)
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    attempt_rows = []
    candidate_dirs = sorted((root / "candidates").glob("candidate_*"))
    for cdir in candidate_dirs:
        cid = int(cdir.name.split("_")[-1])
        sfile = cdir / "summary.json"
        summary = json.loads(sfile.read_text()) if sfile.exists() else {}
        attempts_root = cdir / "pair_attempts"
        any_two = {s: False for s in SCENARIOS}
        best_counts = {s: (0, 0) for s in SCENARIOS}
        n_attempts_with_measurements = 0

        if attempts_root.exists():
            for adir in sorted(attempts_root.glob("attempt_*")):
                mf = adir / "exposure_measurements.csv"
                if not mf.exists():
                    continue
                n_attempts_with_measurements += 1
                res = audit_attempt(mf)
                for scenario, z in res.items():
                    any_two[scenario] |= bool(z["two_epoch"])
                    prev = best_counts[scenario]
                    if min(z["n_early"], z["n_late"]) > min(prev):
                        best_counts[scenario] = (z["n_early"], z["n_late"])
                    attempt_rows.append({
                        "candidate_id": cid,
                        "attempt": adir.name,
                        "scenario": scenario,
                        **z,
                    })

        triggered = [s for s, yes in any_two.items() if yes]
        current_class = summary.get("classification", "")
        if current_class != "INSUFFICIENT_DATA":
            verdict = "CURRENT_NON_INSUFFICIENT_REVIEW"
        elif triggered:
            verdict = "SENSITIVITY_FLAG_TARGETED_RERUN_REQUIRED"
        else:
            verdict = "ROBUST_INSUFFICIENT_ACROSS_TESTED_GATES"

        # DQ context from all saved measurement rows: hard DQ is never relaxed;
        # JUMP_DET/OUTLIER are counted so any candidate dominated by them is visible.
        all_m = []
        if attempts_root.exists():
            for mf in attempts_root.glob("attempt_*/exposure_measurements.csv"):
                try: all_m.append(pd.read_csv(mf))
                except Exception: pass
        if all_m:
            md = pd.concat(all_m, ignore_index=True)
            n_jump = int(md.get("dq_jump_det", pd.Series(False, index=md.index)).map(truthy).sum())
            n_outlier = int(md.get("dq_outlier", pd.Series(False, index=md.index)).map(truthy).sum())
            n_do_not_use = int(md.get("dq_do_not_use", pd.Series(False, index=md.index)).map(truthy).sum())
            n_saturated = int(md.get("dq_saturated", pd.Series(False, index=md.index)).map(truthy).sum())
        else:
            n_jump = n_outlier = n_do_not_use = n_saturated = 0

        row = {
            "candidate_id": cid,
            "current_classification": current_class,
            "current_pm_status": summary.get("pm_status", ""),
            "pair_status": summary.get("pair_status", ""),
            "n_saved_pair_attempts": n_attempts_with_measurements,
            "scenarios_reaching_two_epochs": ";".join(triggered),
            "n_scenarios_reaching_two_epochs": len(triggered),
            "dq_jump_rows": n_jump,
            "dq_outlier_rows": n_outlier,
            "dq_do_not_use_rows": n_do_not_use,
            "dq_saturated_rows": n_saturated,
            "robustness_verdict": verdict,
        }
        for s in SCENARIOS:
            row[f"{s}_two_epoch"] = any_two[s]
            row[f"{s}_best_early"] = best_counts[s][0]
            row[f"{s}_best_late"] = best_counts[s][1]
        summary_rows.append(row)

    sdf = pd.DataFrame(summary_rows).sort_values("candidate_id")
    adf = pd.DataFrame(attempt_rows)
    sdf.to_csv(outdir / "FIRST26_FINAL_ROBUSTNESS.csv", index=False)
    adf.to_csv(outdir / "PAIR_SCENARIO_DETAILS.csv", index=False)

    counts = sdf["robustness_verdict"].value_counts().to_dict()
    flagged = sdf.loc[sdf["robustness_verdict"] == "SENSITIVITY_FLAG_TARGETED_RERUN_REQUIRED", "candidate_id"].tolist()
    lines = [
        "# First26 final robustness audit",
        "",
        "This audit perturbs S/N, forced-Gaussian width, axis ratio, positional uncertainty, reduced-chi2, centroid method, and searches every saved ranked pair attempt.",
        "The <=4 pixel association envelope is held fixed to prevent the demonstrated 643980 neighbour-latching failure. DO_NOT_USE and SATURATED pixels are never relaxed.",
        "",
        "## Verdict counts",
    ]
    for k, v in counts.items():
        lines.append(f"- `{k}`: {v}")
    lines += ["", f"Candidates requiring fresh targeted rerun: {flagged}", ""]
    if flagged:
        lines.append("These candidates are not declared PM detections; the saved diagnostics only show that a plausible gate perturbation can make both epochs available. They require a fresh PM rerun with registration and full uncertainty propagation.")
    else:
        lines.append("No candidate acquires two-epoch astrometric availability under any tested plausible gate perturbation. The current INSUFFICIENT_DATA classifications are therefore robust to these threshold/centroid/pair-selection choices.")
    (outdir / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
