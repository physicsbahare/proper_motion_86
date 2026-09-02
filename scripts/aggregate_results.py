#!/usr/bin/env python3
"""Merge shard outputs into repository-level result tables."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="results")
    args = parser.parse_args()

    source = Path(args.input)
    out = Path(args.output)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "candidates").mkdir(exist_ok=True)

    summaries = []
    for summary_file in source.rglob("summary.json"):
        data = json.loads(summary_file.read_text())
        summaries.append(data)
        cid = int(data["candidate_id"])
        src_dir = summary_file.parent
        dst_dir = out / "candidates" / f"candidate_{cid}"
        shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)

    df = pd.DataFrame(summaries).sort_values("candidate_id") if summaries else pd.DataFrame()
    df.to_csv(out / "ALL45_PM_AUDIT.csv", index=False)

    if len(df):
        counts = (
            df["classification"].fillna("UNKNOWN")
            .value_counts()
            .rename_axis("classification")
            .reset_index(name="n_candidates")
        )
        counts.to_csv(out / "CLASSIFICATION_COUNTS.csv", index=False)
        lines = [
            "# Remote audit result\n",
            f"Candidates completed: **{len(df)}/45**\n",
            "## Classification counts\n",
        ]
        for _, row in counts.iterrows():
            lines.append(f"- `{row['classification']}`: {int(row['n_candidates'])}\n")
        lines += [
            "\n`CONSISTENT_WITH_ZERO` means no significant PM at the achieved precision; it is not proof of physically zero motion.\n",
            "\n`INSUFFICIENT_DATA` is used when the archive lacks a usable independent epoch pair, the target is not detected/centroidable in both epochs, or local registration cannot be supported.\n",
        ]
        (out / "SUMMARY.md").write_text("".join(lines), encoding="utf-8")

    print(df.to_string(index=False) if len(df) else "No summaries found")


if __name__ == "__main__":
    main()
