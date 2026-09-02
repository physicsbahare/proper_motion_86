#!/usr/bin/env python3
"""Run one deterministic shard of the 45-candidate audit."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pm86.pipeline import run_candidate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", default="data/candidates_45.csv")
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--n-shards", type=int, default=5)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

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
        summary = run_candidate(row.to_dict(), out)
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
