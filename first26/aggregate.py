#!/usr/bin/env python3
"""Aggregate the first-26 PM audit while separating science insufficiency from infrastructure failure."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd


INFRA_ERROR_PREFIXES = (
    "ERROR:ImportError:",
    "ERROR:ModuleNotFoundError:",
    "ERROR:ConnectionError:",
    "ERROR:Timeout",
    "ERROR:OSError:",
)


def is_infrastructure_error(status: object) -> bool:
    text = str(status)
    if text.startswith(INFRA_ERROR_PREFIXES):
        return True
    needles = (
        "HTTPFileSystem requires",
        "aiohttp",
        "Temporary failure",
        "Connection reset",
        "timed out",
        "RemoteProtocolError",
        "SSL",
    )
    return text.startswith("ERROR:") and any(x.lower() in text.lower() for x in needles)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-count", type=int, default=26)
    args = parser.parse_args()

    source = Path(args.input)
    out = Path(args.output)
    catalog = pd.read_csv(args.catalog)
    expected_ids = set(catalog["candidate_id"].astype(int))

    if out.exists():
        shutil.rmtree(out)
    (out / "candidates").mkdir(parents=True, exist_ok=True)

    summaries = []
    source_dirs = {}
    for summary_file in source.rglob("summary.json"):
        data = json.loads(summary_file.read_text(encoding="utf-8"))
        cid = int(data["candidate_id"])
        summaries.append(data)
        source_dirs[cid] = summary_file.parent
        shutil.copytree(
            summary_file.parent,
            out / "candidates" / f"candidate_{cid}",
            dirs_exist_ok=True,
        )

    df = pd.DataFrame(summaries)
    if len(df):
        df["candidate_id"] = pd.to_numeric(df["candidate_id"], errors="coerce").astype("Int64")
        df = df.sort_values("candidate_id").reset_index(drop=True)

    # The compact final table retains catalog_* provenance fields; detailed
    # exposure-level evidence remains in results_first26/candidates/.
    df.to_csv(out / "ALL26_PM_AUDIT.csv", index=False)

    if len(df):
        counts = (
            df["classification"].fillna("UNKNOWN")
            .value_counts()
            .rename_axis("classification")
            .reset_index(name="n_candidates")
        )
        counts.to_csv(out / "CLASSIFICATION_COUNTS.csv", index=False)

        reasons = (
            df["reason"].fillna("NONE")
            .value_counts()
            .rename_axis("reason")
            .reset_index(name="n_candidates")
        )
        reasons.to_csv(out / "REASON_COUNTS.csv", index=False)
    else:
        counts = pd.DataFrame(columns=["classification", "n_candidates"])

    # Deep QC: do not allow a hidden remote-data failure to masquerade as a
    # scientifically valid INSUFFICIENT_DATA classification.
    coverage_rows = []
    hidden_infra_candidates = set()
    pair_found_ids = set()
    if len(df) and "pair_status" in df:
        pair_found_ids = set(
            pd.to_numeric(
                df.loc[df["pair_status"].astype(str).eq("PAIR_FOUND"), "candidate_id"],
                errors="coerce",
            ).dropna().astype(int)
        )

    for cid, cdir in source_dirs.items():
        audit_file = cdir / "product_coverage_audit.csv"
        if not audit_file.exists():
            continue
        try:
            audit = pd.read_csv(audit_file)
        except pd.errors.EmptyDataError:
            continue
        if "status" not in audit:
            continue
        statuses = audit["status"].astype(str)
        for status in statuses:
            coverage_rows.append(
                {
                    "candidate_id": cid,
                    "status": status,
                    "is_infrastructure_error": is_infrastructure_error(status),
                }
            )
        if cid in pair_found_ids and len(statuses):
            infra_mask = statuses.map(is_infrastructure_error)
            # If every attempted product failed for infrastructure reasons,
            # this candidate has not actually had an astrometric coverage test.
            if bool(infra_mask.all()):
                hidden_infra_candidates.add(cid)

    coverage_df = pd.DataFrame(coverage_rows)
    coverage_df.to_csv(out / "PRODUCT_ACCESS_QC.csv", index=False)

    got_ids = (
        set(pd.to_numeric(df["candidate_id"], errors="coerce").dropna().astype(int))
        if len(df) else set()
    )
    missing = sorted(expected_ids - got_ids)

    pm_error_ids = []
    if len(df) and "pm_status" in df:
        pm_error_ids = sorted(
            pd.to_numeric(
                df.loc[df["pm_status"].astype(str).eq("ERROR"), "candidate_id"],
                errors="coerce",
            ).dropna().astype(int).tolist()
        )

    problems = []
    if len(expected_ids) != args.expected_count:
        problems.append(
            f"Catalog has {len(expected_ids)} unique candidates, expected {args.expected_count}."
        )
    if missing:
        problems.append("Missing candidate summaries: " + ", ".join(map(str, missing)))
    if pm_error_ids:
        problems.append("Candidate-level pm_status=ERROR: " + ", ".join(map(str, pm_error_ids)))
    if hidden_infra_candidates:
        problems.append(
            "PAIR_FOUND candidates whose product attempts all failed for infrastructure reasons: "
            + ", ".join(map(str, sorted(hidden_infra_candidates)))
        )

    qclines = [
        "# First-26 run quality control\n\n",
        f"Candidate summaries present: **{len(got_ids)}/{len(expected_ids)}**\n\n",
        f"Candidate-level `pm_status == ERROR`: **{len(pm_error_ids)}**\n\n",
        f"Hidden all-product infrastructure failures: **{len(hidden_infra_candidates)}**\n\n",
    ]
    if problems:
        qclines.append("## Problems\n\n")
        qclines.extend(f"- {p}\n" for p in problems)
    else:
        qclines.append(
            "No missing summaries, explicit candidate errors, or hidden all-product "
            "remote-access failures were found.\n"
        )
    (out / "RUN_QC.md").write_text("".join(qclines), encoding="utf-8")

    pair_found = int(df["pair_status"].astype(str).eq("PAIR_FOUND").sum()) if len(df) and "pair_status" in df else 0
    valid_pm = int((~df["classification"].astype(str).eq("INSUFFICIENT_DATA")).sum()) if len(df) and "classification" in df else 0

    lines = [
        "# First 26 candidates — independent PM audit\n\n",
        f"Candidate summaries: **{len(got_ids)}/{len(expected_ids)}**\n\n",
        f"Archive-selected NIRCam epoch pair found: **{pair_found}**\n\n",
        f"Candidates reaching a non-INSUFFICIENT PM classification: **{valid_pm}**\n\n",
        "## Classification counts\n\n",
    ]
    for _, row in counts.iterrows():
        lines.append(f"- `{row['classification']}`: {int(row['n_candidates'])}\n")
    lines += [
        "\nThe source Galaxy/Brown-Dwarf classification is preserved only as metadata "
        "and is not used by the astrometric decision chain.\n\n",
        "`CONSISTENT_WITH_ZERO` means no significant PM at the achieved precision, "
        "not proof of physically zero motion.\n\n",
        "`INSUFFICIENT_DATA` is retained when the archive or source measurement genuinely "
        "cannot support a PM; infrastructure failures are checked separately in RUN_QC.md.\n",
    ]
    (out / "SUMMARY.md").write_text("".join(lines), encoding="utf-8")

    print(df.to_string(index=False) if len(df) else "No summaries found")
    print("\nQC problems:", problems)

    if problems:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
