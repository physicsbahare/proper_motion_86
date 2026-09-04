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

from recenter_measurement import install as install_recenter_patch
from remote_io_patch import install as install_remote_io_patch
import pm86.pipeline as pipeline


# Nine ranked pairs is the deepest robust search for the current remote-I/O path.
# In run 33819456053 candidate 575676 completed attempts 1--9 with no accepted
# forced astrometric fit in either epoch, then the job spent nearly five hours in
# the attempt-10 remote product path and hit GitHub's six-hour hosted-runner hard
# limit.  That is an infrastructure/remote-I/O failure, not useful science.
# The first nine already span the high-priority F444W/F356W/F277W alternatives
# plus multiple shorter-wave cross-filter combinations.  Per-request HTTP bounds
# are installed separately so a single bad CAL product cannot consume a shard.
DEEP_PAIR_ATTEMPTS = 9


def json_safe(value):
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def install_deep_pair_search():
    """First26-only expansion of the ranked epoch-pair search."""
    original = pipeline._rank_epoch_pairs

    def deep_rank(inventory):
        return original(inventory, max_pairs=DEEP_PAIR_ATTEMPTS)

    pipeline._rank_epoch_pairs = deep_rank


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--n-shards", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    install_remote_io_patch()
    install_recenter_patch()
    install_deep_pair_search()

    catalog = pd.read_csv(args.catalog).sort_values("source_row").reset_index(drop=True)
    selected = catalog.iloc[args.shard :: args.n_shards].copy()

    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=True)
    manifest_rows = []

    print(
        f"first26 shard={args.shard}/{args.n_shards} "
        f"n_candidates={len(selected)} ids={selected.candidate_id.astype(int).tolist()} "
        f"max_pair_attempts={DEEP_PAIR_ATTEMPTS}"
    )

    for _, row in selected.iterrows():
        candidate = row.to_dict()
        cid = int(candidate["candidate_id"])
        t0 = time.time()
        print(f"\n=== candidate {cid} ===", flush=True)

        summary = pipeline.run_candidate(candidate, root)

        cdir = root / f"candidate_{cid}"
        input_payload = {k: json_safe(v) for k, v in candidate.items()}
        (cdir / "source_input.json").write_text(
            json.dumps(input_payload, indent=2, sort_keys=True), encoding="utf-8"
        )

        summary_path = cdir / "summary.json"
        saved = json.loads(summary_path.read_text(encoding="utf-8"))
        for key, value in input_payload.items():
            saved[f"catalog_{key}"] = value
        saved["first26_deep_pair_attempt_limit"] = DEEP_PAIR_ATTEMPTS
        saved["first26_faint_astrometry_fallback"] = "bounded_forced_gaussian_v1"
        saved["first26_remote_cal_io"] = "aiohttp_bounded_v1"
        summary_path.write_text(
            json.dumps(saved, indent=2, sort_keys=True), encoding="utf-8"
        )

        elapsed = time.time() - t0
        manifest_rows.append({
            "candidate_id": cid,
            "elapsed_seconds": round(elapsed, 3),
            "pair_status": saved.get("pair_status"),
            "pm_status": saved.get("pm_status"),
            "classification": saved.get("classification"),
            "reason": saved.get("reason"),
        })
        print(
            f"candidate {cid}: {saved.get('classification')} / "
            f"{saved.get('pm_status')} in {elapsed:.1f}s", flush=True
        )

    pd.DataFrame(manifest_rows).to_csv(root / "shard_manifest.csv", index=False)


if __name__ == "__main__":
    main()
