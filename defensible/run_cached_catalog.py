"""CLI for the epoch-cached exhaustive proper-motion audit."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from defensible.cached_exhaustive_pm import run_candidate_cached


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--output-root", required=True)
    ap.add_argument("--index", type=int, required=True)
    args = ap.parse_args()

    cat = pd.read_csv(args.catalog)
    row = cat.iloc[int(args.index)]
    candidate_id = int(row["candidate_id"])
    ra_deg = float(row["ra_deg"])
    dec_deg = float(row["dec_deg"])

    final = run_candidate_cached(candidate_id, ra_deg, dec_deg, Path(args.output_root))
    out = Path(args.output_root)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([final]).to_csv(out / f"summary_{candidate_id}.csv", index=False)
    print(final, flush=True)


if __name__ == "__main__":
    main()
