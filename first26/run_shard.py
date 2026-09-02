#!/usr/bin/env python3
"""Run an independent PM audit shard for the first 26 090630 candidates.

The source SED/model classification is preserved as metadata only and is never
used to choose epochs, centroids, registration, or the PM classification.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

from first26.recenter_measurement import install as install_recenter_patch
from pm86.pipeline import run_candidate


def json_safe(value):
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--n-shards", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # The base pipeline first measured only at the catalog coordinate.  A real
    # mover can therefore be several pixels away and be rejected before PM is
    # even estimated.  Install a first26-only local recentering layer that
    # searches up to 12 pixels around the catalog prediction and preserves the
    # search offset in the exposure evidence.
    install_recenter_patch()

    catalog = pd.read_csv(args.catalog).sort_values("source_row").reset_index(drop=True)
    selected = catalog.iloc[args.shard :: args.n_shards].copy()

    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=True)
    manifest_rows = []

    print(
        f"first26 shard={args.shard}/{args.n_shards} "
        f"n_candidates={len(selected)} ids={selected.candidate_id.astype(int).tolist()}"
    )

    for _, row in selected.iterrows():
        candidate = row.to_dict()
        cid = int(candidate["candidate_id"])
        t0 = time.time()
        print(f"\n=== candidate {cid} ===", flush=True)

        summary = run_candidate(candidate, root)

        cdir = root / f"candidate_{cid}"
        input_payload = {k: json_safe(v) for k, v in candidate.items()}
        (cdir / "source_input.json").write_text(
            json.dumps(input_payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )

        summary_path = cdir / "summary.json"
        saved = json.loads(summary_path.read_text(encoding="utf-8"))
        for key, value in input_payload.items():
            saved[f"catalog_{key}"] = value
        summary_path.write_text(
            json.dumps(saved, indent=2, sort_keys=True),
            encoding="utf-8",
        )

        elapsed = time.time() - t0
        manifest_rows.append(
            {
                "candidate_id": cid,
                "elapsed_seconds": round(elapsed, 3),
                "pair_status": saved.get("pair_status"),
                "pm_status": saved.get("pm_status"),
                "classification": saved.get("classification"),
                "reason": saved.get("reason"),
            }
        )
        print(
            f"candidate {cid}: {saved.get('classification')} / "
            f"{saved.get('pm_status')} in {elapsed:.1f}s",
            flush=True,
        )

    pd.DataFrame(manifest_rows).to_csv(root / "shard_manifest.csv", index=False)


if __name__ == "__main__":
    main()
