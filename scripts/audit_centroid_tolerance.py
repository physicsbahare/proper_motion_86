#!/usr/bin/env python3
"""Audit whether a fixed 60 mas centroid check depends on analysis method."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pm86.tolerance import (
    DEFAULT_TOLERANCE_MAS,
    add_tolerance_flags,
    collect_tolerance_rows,
    summarize_by_filter,
    summarize_tolerance,
    tolerance_sweep,
)


def fmt(value):
    if pd.isna(value):
        return "n/a"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default="results")
    parser.add_argument("--output", default="results/centroid_tolerance")
    parser.add_argument("--tolerance-mas", type=float, default=DEFAULT_TOLERANCE_MAS)
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    detail = collect_tolerance_rows(args.results)
    detail = add_tolerance_flags(detail, args.tolerance_mas) if len(detail) else detail
    summary = summarize_tolerance(detail, args.tolerance_mas)
    by_filter = summarize_by_filter(detail, args.tolerance_mas)
    sweep = tolerance_sweep(detail)

    detail.to_csv(out / "CENTROID_TOLERANCE_DETAIL.csv", index=False)
    pd.DataFrame([summary]).to_csv(out / "CENTROID_TOLERANCE_SUMMARY.csv", index=False)
    by_filter.to_csv(out / "CENTROID_TOLERANCE_BY_FILTER.csv", index=False)
    sweep.to_csv(out / "CENTROID_TOLERANCE_SWEEP.csv", index=False)

    lines = [
        "# Centroid tolerance audit\n\n",
        f"Proposed verifier tolerance: **{args.tolerance_mas:.1f} mas**.\n\n",
        "This is an empirical method-sensitivity check, not an assumption that 60 mas is correct. "
        "It compares the two independently retained target centroid methods (2-D Gaussian and "
        "center-of-mass) and also measures how much local astrometric registration changes the "
        "raw gWCS position. Only clean target measurements with S/N >= 3 are used in the headline "
        "statistics.\n\n",
        "## Headline numbers\n\n",
    ]
    for key in [
        "eligible_snr3_rows",
        "eligible_candidates",
        "method_comparable_rows",
        "method_rows_over_tolerance",
        "method_fraction_within_tolerance",
        "method_sep_p95_mas",
        "method_sep_p99_mas",
        "method_sep_max_mas",
        "registration_comparable_rows",
        "registration_rows_over_tolerance",
        "registration_fraction_within_tolerance",
        "registration_shift_p95_mas",
        "registration_shift_p99_mas",
        "registration_shift_max_mas",
        "assessment",
    ]:
        if key in summary:
            lines.append(f"- `{key}`: {fmt(summary[key])}\n")

    lines += [
        "\n## How to interpret this\n\n",
        "A `method_rows_over_tolerance` count above zero means that the same S/N>=3 source can "
        "move by more than the proposed tolerance solely because a different legitimate centroid "
        "algorithm is used. A `registration_rows_over_tolerance` count above zero means that raw "
        "gWCS versus local-registration frame choice can also move a position by more than the "
        "tolerance. Either result is evidence that a hard 60 mas rule should be method-aware or "
        "calibrated more broadly.\n\n",
        "The audit directly tests 2DG versus COM and raw-gWCS versus local affine registration. "
        "It does not claim to replace a dedicated empirical/ePSF fit comparison; if the measured "
        "spread approaches 60 mas, an ePSF/PSF-fit robustness test should be added before signing "
        "off the fixed threshold.\n",
    ]
    (out / "CENTROID_TOLERANCE_REPORT.md").write_text("".join(lines), encoding="utf-8")
    print("\n".join(lines[:8]))
    print(pd.DataFrame([summary]).to_string(index=False))


if __name__ == "__main__":
    main()
