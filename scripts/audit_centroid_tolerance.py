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
    summarize_by_group,
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
    parser.add_argument("--reference-centroids", default=None)
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    detail = collect_tolerance_rows(args.results, args.reference_centroids)
    detail = add_tolerance_flags(detail, args.tolerance_mas) if len(detail) else detail
    summary = summarize_tolerance(detail, args.tolerance_mas)
    by_filter = summarize_by_group(detail, "filter", args.tolerance_mas)
    by_source = summarize_by_group(detail, "source", args.tolerance_mas)
    sweep = tolerance_sweep(detail)

    detail.to_csv(out / "CENTROID_TOLERANCE_DETAIL.csv", index=False)
    pd.DataFrame([summary]).to_csv(out / "CENTROID_TOLERANCE_SUMMARY.csv", index=False)
    by_filter.to_csv(out / "CENTROID_TOLERANCE_BY_FILTER.csv", index=False)
    by_source.to_csv(out / "CENTROID_TOLERANCE_BY_SOURCE.csv", index=False)
    sweep.to_csv(out / "CENTROID_TOLERANCE_SWEEP.csv", index=False)

    reference = by_source[by_source.get("source", pd.Series(dtype=str)).astype(str).eq("REFERENCE_282040")] if len(by_source) else pd.DataFrame()

    lines = [
        "# Centroid tolerance audit\n\n",
        f"Proposed verifier tolerance: **{args.tolerance_mas:.1f} mas**.\n\n",
        "This is an empirical method-sensitivity check, not an assumption that 60 mas is correct. "
        "It compares 2-D Gaussian and center-of-mass centroids measured on the same detector image "
        "and separately measures local-registration frame shifts. Only clean S/N >= 3 measurements "
        "enter the headline statistics.\n\n",
    ]

    if len(reference):
        r = reference.iloc[0]
        lines += [
            "## Independent 282040 reference pre-check\n\n",
            f"The previously completed 282040 run contributes {int(r['method_comparable_rows'])} S/N>=3 exposures. "
            f"At {args.tolerance_mas:.0f} mas, {int(r['method_rows_over_tolerance'])} exposure(s) have a 2DG-COM "
            f"separation larger than the threshold; the maximum separation is {r['method_sep_max_mas']:.1f} mas "
            f"and the 95th percentile is {r['method_sep_p95_mas']:.1f} mas.\n\n",
        ]

    lines += ["## Headline numbers including available reference/calibration rows\n\n"]
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
        "\n## Interpretation\n\n",
        "A non-zero `method_rows_over_tolerance` is direct evidence that the same real S/N>=3 source can "
        "differ by more than the proposed threshold solely because a different defensible centroid method "
        "is used. A non-zero `registration_rows_over_tolerance` indicates analogous astrometric-frame "
        "sensitivity. In either case, a hard pass/fail radius should be made method-aware or calibrated "
        "with a broader method set rather than treated as universally neutral.\n\n",
        "This audit directly tests 2DG versus COM and raw-gWCS versus local affine registration. It does "
        "not replace a dedicated empirical/ePSF comparison. Because the 282040 reference already reaches "
        "the vicinity of 60 mas, an ePSF/PSF-fit robustness comparison would be the strongest final check "
        "before approving 60 mas as a method-independent hard cutoff.\n",
    ]
    (out / "CENTROID_TOLERANCE_REPORT.md").write_text("".join(lines), encoding="utf-8")
    print(pd.DataFrame([summary]).to_string(index=False))
    if len(reference):
        print("\n282040 reference:\n", reference.to_string(index=False))


if __name__ == "__main__":
    main()
