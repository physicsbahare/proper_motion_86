"""Scheduled, bounded JWST archive recheck of separate Ali and Dave samples.

Only the original candidate coordinates are used; previous classifications are
NOT pipeline inputs. One small deterministic batch is processed per day.
Unchanged archive inventories are skipped after the first pass to avoid
repeatedly downloading the same large CAL image cutouts without new evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import date, timezone, datetime
from pathlib import Path

import pandas as pd
from defensible.cached_exhaustive_pm import run_candidate_cached
from pm86.archive import query_candidate_inventory

START = date(2026, 10, 13)
METHOD_VERSION = "jwst-li-fynbo-covariance-diagnostic-v1"
CATALOGS = {
    "ali": ("data/candidates_ali_cosmos2025_45.csv", 45),
    # Published Dave sample: 125 sources, of which 387132 and 1651050
    # already have motion estimates. Here we recheck the other 123.
    "dave123": ("data/candidates_dave_remaining_123.csv", 123),
}


def inventory_digest(df):
    if df is None or df.empty:
        return "EMPTY"
    preferred = ["obs_id", "obsid", "filters", "instrument_name",
                 "obs_collection", "t_min", "t_max", "proposal_id", "target_name"]
    cols = [c for c in preferred if c in df.columns]
    if not cols:
        cols = sorted(df.columns)
    norm = df[cols].astype(str).sort_values(cols).reset_index(drop=True)
    return hashlib.sha256(norm.to_csv(index=False).encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", choices=tuple(CATALOGS), required=True)
    ap.add_argument("--slot", type=int, required=True)
    ap.add_argument("--slots-per-day", type=int, default=5)
    ap.add_argument("--output-root", default="round_robin_run")
    ap.add_argument("--state-root", default="results_li2026_recheck")
    a = ap.parse_args()

    path, required_n = CATALOGS[a.sample]
    cat = pd.read_csv(path)
    if len(cat) != required_n or cat["candidate_id"].nunique() != required_n:
        raise RuntimeError(f"Wrong candidate sample: {path}, expected {required_n} unique IDs")
    if not (0 <= a.slot < a.slots_per_day):
        raise ValueError("slot must be within the daily batch")

    day = max(0, (datetime.now(timezone.utc).date() - START).days)
    index = (day * a.slots_per_day + a.slot) % len(cat)
    row = cat.iloc[index]
    cid, ra, dec = int(row.candidate_id), float(row.ra_deg), float(row.dec_deg)

    out = Path(a.output_root) / a.sample
    out.mkdir(parents=True, exist_ok=True)
    previous_path = Path(a.state_root) / a.sample / f"summary_{cid}.json"
    previous = json.loads(previous_path.read_text()) if previous_path.exists() else {}

    # Query inventory first; only rerun Stage-2 extraction when it actually
    # changes or the analysis version differs. A source with no independent
    # epochs can remain unresolved indefinitely; that is scientific, not a
    # reason to manufacture a stationary label.
    inventory = query_candidate_inventory(cid, ra, dec)
    fingerprint = inventory_digest(inventory)
    unchanged = (previous.get("archive_fingerprint") == fingerprint and
                 previous.get("method_version") == METHOD_VERSION)
    if unchanged:
        result = dict(previous)
        result["recheck_outcome"] = "SKIPPED_ARCHIVE_UNCHANGED"
    else:
        result = run_candidate_cached(cid, ra, dec, out)
        result["recheck_outcome"] = "REMEASURED"
    result.update(
        sample=a.sample, candidate_id=cid, ra_deg=ra, dec_deg=dec,
        archive_fingerprint=fingerprint, method_version=METHOD_VERSION,
        last_checked_utc=datetime.now(timezone.utc).isoformat(),
        # Do not confuse theoretical chi2 significance with an empirical
        # JWST null calibration or with the Gaia threshold of D²=7.097.
        jwst_null_calibration_validated=False,
    )

    dest = out / f"summary_{cid}.json"
    dest.write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(f"{a.sample}: {index+1}/{len(cat)}, candidate {cid}: "
          f"{result.get('classification', 'ERROR')}, {result['recheck_outcome']}", flush=True)


if __name__ == "__main__":
    main()
