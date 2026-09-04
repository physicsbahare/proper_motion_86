#!/usr/bin/env python3
"""Run the validated first26-style PM search on the remaining 19 090630 candidates."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "first26"))

from recenter_measurement import install as install_recenter_patch
from remote_io_patch import install as install_remote_io_patch
import pm86.pipeline as pipeline

DEEP_PAIR_ATTEMPTS = 9


def json_safe(value):
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def install_deep_pair_search():
    original = pipeline._rank_epoch_pairs
    def deep_rank(inventory):
        return original(inventory, max_pairs=DEEP_PAIR_ATTEMPTS)
    pipeline._rank_epoch_pairs = deep_rank


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--catalog", required=True)
    p.add_argument("--shard", type=int, required=True)
    p.add_argument("--n-shards", type=int, required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    install_remote_io_patch()
    install_recenter_patch()
    install_deep_pair_search()

    catalog = pd.read_csv(args.catalog).sort_values("source_row").reset_index(drop=True)
    selected = catalog.iloc[args.shard :: args.n_shards].copy()
    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=True)
    manifest = []

    print(f"remaining19 shard={args.shard}/{args.n_shards} ids={selected.candidate_id.astype(int).tolist()}")
    for _, row in selected.iterrows():
        candidate = row.to_dict()
        cid = int(candidate["candidate_id"])
        t0 = time.time()
        print(f"\n=== candidate {cid} ===", flush=True)
        pipeline.run_candidate(candidate, root)

        cdir = root / f"candidate_{cid}"
        payload = {k: json_safe(v) for k, v in candidate.items()}
        (cdir / "source_input.json").write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        sp = cdir / "summary.json"
        saved = json.loads(sp.read_text(encoding="utf-8"))
        for k, v in payload.items():
            saved[f"catalog_{k}"] = v
        saved["pm_search_pair_attempt_limit"] = DEEP_PAIR_ATTEMPTS
        saved["pm_search_faint_astrometry_fallback"] = "bounded_forced_gaussian_v1"
        saved["pm_search_remote_cal_io"] = "aiohttp_bounded_v1"
        sp.write_text(json.dumps(saved, indent=2, sort_keys=True), encoding="utf-8")
        manifest.append({
            "candidate_id": cid,
            "elapsed_seconds": round(time.time()-t0, 3),
            "pair_status": saved.get("pair_status"),
            "pm_status": saved.get("pm_status"),
            "classification": saved.get("classification"),
            "reason": saved.get("reason"),
        })
        print(f"candidate {cid}: {saved.get('classification')} / {saved.get('pm_status')}", flush=True)

    pd.DataFrame(manifest).to_csv(root / "shard_manifest.csv", index=False)


if __name__ == "__main__":
    main()
