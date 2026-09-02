#!/usr/bin/env python3
"""Run one deterministic shard of the 45-candidate audit.

This runner adds two safeguards learned from the remote audit itself:
1. retry transient MAST inventory-query disconnects, and
2. if forced photometry at the catalog coordinate is not secure, allow the
   same bounded motion-aware recentering used by the independent first26 audit.

The recentering does not alter field controls or PM inference.  Its offset and
whether it was used are retained in each exposure-measurement table so the
result remains auditable.
"""

from __future__ import annotations

import argparse
import importlib.util
import time
from pathlib import Path

import pandas as pd

import pm86.archive as archive
import pm86.pipeline as pipeline


def _install_inventory_retry(attempts: int = 5) -> None:
    """Retry transient MAST query_region failures without hiding real errors."""
    base_query = archive.query_candidate_inventory

    def query_with_retry(candidate_id: int, ra_deg: float, dec_deg: float):
        last = None
        for attempt in range(attempts):
            try:
                return base_query(candidate_id, ra_deg, dec_deg)
            except Exception as exc:
                last = exc
                if attempt + 1 < attempts:
                    # Small exponential-ish backoff; MAST occasionally closes a
                    # connection while many shards query simultaneously.
                    time.sleep(3 * (attempt + 1))
        raise RuntimeError(
            f"MAST inventory query failed after {attempts} attempts: {last}"
        ) from last

    # pipeline imported the function directly, so patch that reference.
    pipeline.query_candidate_inventory = query_with_retry


def _install_motion_aware_recenter() -> None:
    """Load the already-audited bounded recentering helper used by first26."""
    helper = Path(__file__).resolve().parents[1] / "first26" / "recenter_measurement.py"
    if not helper.exists():
        raise FileNotFoundError(f"Missing recenter helper: {helper}")
    spec = importlib.util.spec_from_file_location("pm86_run45_recenter", helper)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load recenter helper: {helper}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.install()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", default="data/candidates_45.csv")
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--n-shards", type=int, default=5)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    _install_inventory_retry()
    _install_motion_aware_recenter()

    catalog = pd.read_csv(args.catalog).sort_values("candidate_id").reset_index(drop=True)
    selected = catalog.iloc[[i for i in range(len(catalog)) if i % args.n_shards == args.shard]].copy()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    print(f"Shard {args.shard}/{args.n_shards}: {len(selected)} candidates")
    for _, row in selected.iterrows():
        cid = int(row["candidate_id"])
        print("=" * 80)
        print(f"Candidate {cid}")
        summary = pipeline.run_candidate(row.to_dict(), out)
        summaries.append(summary)
        print(
            f"classification={summary.get('classification')} "
            f"status={summary.get('pm_status')} "
            f"reason={summary.get('reason', '')}"
        )
        pd.DataFrame(summaries).to_csv(out / "shard_summary.csv", index=False)

    print("Shard complete")


if __name__ == "__main__":
    main()
