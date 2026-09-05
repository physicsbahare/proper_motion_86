"""CLI wrapper for exhaustive proper-motion audits."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from defensible.exhaustive_pm import run_candidate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--output-root", required=True)
    ap.add_argument("--index", type=int, default=None, help="Run one zero-based catalog row (for CI sharding).")
    args = ap.parse_args()

    df = pd.read_csv(args.catalog)
    required = {"candidate_id", "ra_deg", "dec_deg"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"missing required columns: {sorted(missing)}")

    Path(args.output_root).mkdir(parents=True, exist_ok=True)
    rows = df.iloc[[args.index]] if args.index is not None else df

    summaries = []
    for _, row in rows.iterrows():
        summaries.append(run_candidate(
            int(row["candidate_id"]),
            float(row["ra_deg"]),
            float(row["dec_deg"]),
            args.output_root,
        ))

    out = pd.DataFrame(summaries)
    name = f"summary_{args.index:03d}.csv" if args.index is not None else "ALL_CANDIDATE_SUMMARIES.csv"
    out.to_csv(Path(args.output_root) / name, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
