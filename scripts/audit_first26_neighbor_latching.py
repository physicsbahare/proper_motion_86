#!/usr/bin/env python3
from __future__ import annotations

import argparse, io, json, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

DEFAULT_OLD_REF = "5426f33b0f4e92b36e4278bb0943035dd902a518"
ASSOC_LIMIT_PIX = 4.0
SEVERE_OFFSET_PIX = 8.0


def git_show(ref: str, path: str) -> str | None:
    p = subprocess.run(["git", "show", f"{ref}:{path}"], text=True, capture_output=True)
    return p.stdout if p.returncode == 0 else None


def finite_num(s):
    return pd.to_numeric(s, errors="coerce")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default="data/candidates_first26_090630.csv")
    ap.add_argument("--current-results", default="results_first26")
    ap.add_argument("--old-ref", default=DEFAULT_OLD_REF)
    ap.add_argument("--output", default="results_first26/neighbor_latching_audit")
    args = ap.parse_args()

    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    cat = pd.read_csv(args.catalog)
    rows = []
    detail = []

    for cid in cat.candidate_id.astype(int):
        base = f"results_first26/candidates/candidate_{cid}"
        old_summary_txt = git_show(args.old_ref, f"{base}/summary.json")
        old_exp_txt = git_show(args.old_ref, f"{base}/exposure_measurements.csv")
        cur_summary_path = Path(args.current_results) / "candidates" / f"candidate_{cid}" / "summary.json"
        old_summary = json.loads(old_summary_txt) if old_summary_txt else {}
        cur_summary = json.loads(cur_summary_path.read_text()) if cur_summary_path.exists() else {}

        n_large = n_severe = n_large_astrom = 0
        max_offset = np.nan
        large_epochs = set()
        if old_exp_txt:
            d = pd.read_csv(io.StringIO(old_exp_txt))
            off = finite_num(d.get("target_seed_offset_pix", pd.Series(np.nan, index=d.index)))
            rec = d.get("target_seed_recentered", pd.Series(False, index=d.index)).astype(str).str.lower().isin(["true","1"])
            large = rec & off.gt(ASSOC_LIMIT_PIX)
            severe = rec & off.gt(SEVERE_OFFSET_PIX)
            snr = finite_num(d.get("clean_snr_err", pd.Series(np.nan, index=d.index)))
            x = finite_num(d.get("x_2dg", pd.Series(np.nan, index=d.index)))
            y = finite_num(d.get("y_2dg", pd.Series(np.nan, index=d.index)))
            old_astrom_like = large & snr.ge(3.0) & x.notna() & y.notna()
            n_large, n_severe, n_large_astrom = int(large.sum()), int(severe.sum()), int(old_astrom_like.sum())
            max_offset = float(off.max()) if off.notna().any() else np.nan
            if old_astrom_like.any():
                large_epochs = set(d.loc[old_astrom_like, "epoch_label"].astype(str))
            for _, r in d.loc[large].iterrows():
                detail.append({
                    "candidate_id": cid, "epoch_label": r.get("epoch_label"), "filter": r.get("filter"),
                    "filename": r.get("filename"), "offset_pix": r.get("target_seed_offset_pix"),
                    "clean_snr_err": r.get("clean_snr_err"), "x_2dg": r.get("x_2dg"), "y_2dg": r.get("y_2dg"),
                    "old_astrometry_like": bool(pd.notna(r.get("x_2dg")) and pd.notna(r.get("y_2dg")) and pd.to_numeric(pd.Series([r.get("clean_snr_err")]), errors="coerce").iloc[0] >= 3),
                })

        old_pm = old_summary.get("pm_status") == "PM_MEASURED"
        both_epochs_large = {"early","late"}.issubset(large_epochs)
        if old_pm and both_epochs_large:
            verdict = "CONFIRMED_OLD_PM_CONTAMINATED_BY_NEIGHBOR_LATCHING"
        elif n_large_astrom > 0:
            verdict = "HISTORICAL_LARGE_OFFSET_ASTROMETRY_PRESENT_BUT_NO_PM_RESULT"
        elif n_large > 0:
            verdict = "HISTORICAL_LARGE_OFFSET_RECENTER_PRESENT_NOT_ASTROMETRIC"
        else:
            verdict = "NO_LARGE_OFFSET_LATCHING_EVIDENCE"

        rows.append({
            "candidate_id": cid,
            "old_classification": old_summary.get("classification"),
            "old_pm_status": old_summary.get("pm_status"),
            "current_classification": cur_summary.get("classification"),
            "current_pm_status": cur_summary.get("pm_status"),
            "n_old_recenter_gt4pix": n_large,
            "n_old_recenter_gt8pix": n_severe,
            "n_old_gt4pix_astrometry_like": n_large_astrom,
            "max_old_recenter_offset_pix": max_offset,
            "large_offset_astrometry_in_both_epochs": both_epochs_large,
            "audit_verdict": verdict,
        })

    summary = pd.DataFrame(rows).sort_values("candidate_id")
    summary.to_csv(out / "FIRST26_NEIGHBOR_LATCHING_AUDIT.csv", index=False)
    pd.DataFrame(detail).to_csv(out / "LARGE_OFFSET_EXPOSURES.csv", index=False)

    counts = summary.audit_verdict.value_counts().to_dict()
    critical = summary[summary.audit_verdict == "CONFIRMED_OLD_PM_CONTAMINATED_BY_NEIGHBOR_LATCHING"].candidate_id.tolist()
    review = summary[summary.audit_verdict.str.startswith("HISTORICAL_")].candidate_id.tolist()
    md = ["# First26 historical neighbour-latching audit", "", f"Historical reference: `{args.old_ref}`", f"Association limit tested: > {ASSOC_LIMIT_PIX:.1f} pix", "", "## Verdict counts"]
    md += [f"- `{k}`: {v}" for k,v in counts.items()]
    md += ["", f"Confirmed old PM contaminated by large-offset neighbour latching: {critical or 'none'}", f"Other candidates with historical >4-pix recentering to review: {review or 'none'}", "", "A historical large-offset recenter is not itself a PM detection. The critical condition requires an old PM measurement plus astrometry-like >4-pix associations in both epochs. Current results remain the science result unless this audit identifies a separate reproducible issue."]
    (out / "SUMMARY.md").write_text("\n".join(md) + "\n")
    print(summary.to_string(index=False))
    print("counts", counts)

if __name__ == "__main__":
    main()
