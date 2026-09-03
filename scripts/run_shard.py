#!/usr/bin/env python3
"""Run one deterministic shard of the 45-candidate proper-motion audit.

The core pipeline now owns ranked epoch-pair fallback.  This shard driver must
call it exactly once per candidate; wrapping it in a second fallback loop would
multiply the work (up to 6 x 8 pair trials) and can hit the GitHub-hosted runner
time limit before a candidate summary is written.
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
    """Retry transient MAST inventory failures without changing science logic."""
    base_query = archive.query_candidate_inventory

    def query_with_retry(candidate_id: int, ra_deg: float, dec_deg: float):
        last = None
        for attempt in range(attempts):
            try:
                return base_query(candidate_id, ra_deg, dec_deg)
            except Exception as exc:
                last = exc
                if attempt + 1 < attempts:
                    time.sleep(3 * (attempt + 1))
        raise RuntimeError(
            f"MAST inventory query failed after {attempts} attempts: {last}"
        ) from last

    pipeline.query_candidate_inventory = query_with_retry


def _install_motion_aware_recenter() -> None:
    """Install the tested motion-aware recentering helper used by first26."""
    helper = Path(__file__).resolve().parents[1] / "first26" / "recenter_measurement.py"
    spec = importlib.util.spec_from_file_location("pm86_run45_recenter", helper)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load recenter helper: {helper}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.install()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", default="data/candidates_45.csv")
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--n-shards", type=int, default=5)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    _install_inventory_retry()
    _install_motion_aware_recenter()

    catalog = pd.read_csv(args.catalog).sort_values("candidate_id").reset_index(drop=True)
    selected = catalog.iloc[
        [i for i in range(len(catalog)) if i % args.n_shards == args.shard]
    ].copy()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    print(f"Shard {args.shard}/{args.n_shards}: {len(selected)} candidates")
    for _, row in selected.iterrows():
        cid = int(row["candidate_id"])
        print("=" * 80)
        print(f"Candidate {cid}")

        # IMPORTANT: pipeline.run_candidate already ranks and tries up to
        # MAX_PAIR_ATTEMPTS epoch pairs.  Do not add another pair loop here.
        summary = pipeline.run_candidate(row.to_dict(), out)
        summaries.append(summary)

        print(
            f"FINAL classification={summary.get('classification')} "
            f"status={summary.get('pm_status')} "
            f"reason={summary.get('reason', '')} "
            f"n_pair_attempts={summary.get('n_pair_attempts', '')}"
        )
        pd.DataFrame(summaries).to_csv(out / "shard_summary.csv", index=False)

    print("Shard complete")


if __name__ == "__main__":
    main()
