#!/usr/bin/env python3
"""Run one deterministic shard of the 45-candidate audit with robust pair fallback."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import tempfile
import time
from pathlib import Path

import pandas as pd

import pm86.archive as archive
import pm86.pipeline as pipeline
from pm86.pair_candidates import ranked_epoch_pairs


def _install_inventory_retry(attempts: int = 5) -> None:
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
        raise RuntimeError(f"MAST inventory query failed after {attempts} attempts: {last}") from last

    pipeline.query_candidate_inventory = query_with_retry


def _install_motion_aware_recenter() -> None:
    helper = Path(__file__).resolve().parents[1] / "first26" / "recenter_measurement.py"
    spec = importlib.util.spec_from_file_location("pm86_run45_recenter", helper)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load recenter helper: {helper}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.install()


def _attempt_score(summary: dict) -> tuple:
    """Prefer real PM results, then attempts with secure astrometry in both epochs."""
    cls = str(summary.get("classification", ""))
    status = str(summary.get("pm_status", ""))
    n1 = int(summary.get("n_astrometric_early") or 0)
    n2 = int(summary.get("n_astrometric_late") or 0)
    if status != "ERROR" and cls != "INSUFFICIENT_DATA":
        return (0, 0, 0)
    if status != "ERROR" and n1 > 0 and n2 > 0:
        return (1, -(n1 + n2), -min(n1, n2))
    if status != "ERROR":
        return (2, -(n1 + n2), -min(n1, n2))
    return (3, 0, 0)


def _run_with_pair_fallback(candidate: dict, out: Path) -> dict:
    cid = int(candidate["candidate_id"])
    ra = float(candidate["ra_deg"])
    dec = float(candidate["dec_deg"])

    # Discover inventory once, with retry, then reuse it across pair attempts.
    inventory = pipeline.query_candidate_inventory(cid, ra, dec)
    pairs = ranked_epoch_pairs(inventory, max_pairs=6)
    if not pairs:
        # Preserve the base pipeline's exact no-pair scientific reason.
        original_query = pipeline.query_candidate_inventory
        pipeline.query_candidate_inventory = lambda *_args, **_kwargs: inventory.copy()
        try:
            return pipeline.run_candidate(candidate, out)
        finally:
            pipeline.query_candidate_inventory = original_query

    attempts = []
    original_query = pipeline.query_candidate_inventory
    original_choose = pipeline.choose_epoch_pair

    with tempfile.TemporaryDirectory(prefix=f"pm86_{cid}_") as td:
        tmpbase = Path(td)
        for i, pair in enumerate(pairs, start=1):
            attempt_root = tmpbase / f"attempt_{i}"
            pipeline.query_candidate_inventory = lambda *_args, **_kwargs: inventory.copy()
            pipeline.choose_epoch_pair = lambda _inventory, p=pair: p
            summary = pipeline.run_candidate(candidate, attempt_root)
            attempts.append((i, pair, summary, attempt_root / f"candidate_{cid}"))
            print(
                f"  pair attempt {i}: {pair['filter_early']}->{pair['filter_late']} "
                f"baseline={pair['baseline_days_inventory']:.1f} d "
                f"classification={summary.get('classification')} "
                f"reason={summary.get('reason', '')}"
            )
            if _attempt_score(summary)[0] == 0:
                break

        pipeline.query_candidate_inventory = original_query
        pipeline.choose_epoch_pair = original_choose

        best = min(attempts, key=lambda x: _attempt_score(x[2]))
        best_i, best_pair, best_summary, best_dir = best
        final_dir = out / f"candidate_{cid}"
        if final_dir.exists():
            shutil.rmtree(final_dir)
        shutil.copytree(best_dir, final_dir)

        # Preserve alternate-pair evidence without introducing nested summary.json
        # files that would confuse repository aggregation.
        alt_root = final_dir / "pair_attempts"
        alt_root.mkdir(exist_ok=True)
        manifest = []
        for i, pair, summary, cdir in attempts:
            manifest.append({
                "attempt": i,
                "selected": i == best_i,
                "pair_type": pair["pair_type"],
                "filter_early": pair["filter_early"],
                "filter_late": pair["filter_late"],
                "baseline_days": pair["baseline_days_inventory"],
                "classification": summary.get("classification"),
                "pm_status": summary.get("pm_status"),
                "reason": summary.get("reason", ""),
                "n_astrometric_early": summary.get("n_astrometric_early"),
                "n_astrometric_late": summary.get("n_astrometric_late"),
            })
            if i == best_i:
                continue
            dst = alt_root / f"attempt_{i}_{pair['filter_early']}_{pair['filter_late']}"
            shutil.copytree(cdir, dst, dirs_exist_ok=True)
            nested_summary = dst / "summary.json"
            if nested_summary.exists():
                nested_summary.rename(dst / "attempt_summary.json")
        pd.DataFrame(manifest).to_csv(final_dir / "pair_attempt_manifest.csv", index=False)

        best_summary = dict(best_summary)
        best_summary["pair_attempt_count"] = len(attempts)
        best_summary["selected_pair_attempt"] = best_i
        best_summary["pair_fallback_used"] = best_i != 1
        (final_dir / "summary.json").write_text(
            json.dumps(best_summary, indent=2, sort_keys=True, default=str), encoding="utf-8"
        )
        return best_summary


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
        summary = _run_with_pair_fallback(row.to_dict(), out)
        summaries.append(summary)
        print(
            f"FINAL classification={summary.get('classification')} "
            f"status={summary.get('pm_status')} "
            f"reason={summary.get('reason', '')}"
        )
        pd.DataFrame(summaries).to_csv(out / "shard_summary.csv", index=False)

    print("Shard complete")


if __name__ == "__main__":
    main()
