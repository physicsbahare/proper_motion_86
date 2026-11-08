#!/usr/bin/env python3
"""Pair stationary and moving fake-source recoveries from an existing Actions ZIP.

This uses the *same* injected sky site in stationary and moving modes and
subtracts both recovered PM vectors, then compares the difference to the known
injected motion. It is NOT an independent JWST PM measurement; results are
conditional on both successful recoveries and the existing unvalidated
cross-filter template PSFs.
"""
from __future__ import annotations

import argparse
from io import BytesIO
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def summarize(trials, injected_mu):
    required = {"trial", "mode", "mu_east", "mu_north", "status"}
    if not required.issubset(trials.columns):
        raise ValueError("Missing trial columns: " + str(sorted(required - set(trials.columns))))
    if trials.duplicated(["trial", "mode"]).any():
        raise ValueError("Duplicate sky-site/mode records")
    moving = trials[trials["mode"] == "moving"].copy()
    stationary = trials[trials["mode"] == "stationary"].copy()
    if len(moving) != len(stationary):
        raise ValueError("Stationary/moving sample sizes differ")
    pairs = stationary.merge(moving, on="trial", suffixes=("_stationary", "_moving"),
                             validate="one_to_one")
    ok = (pairs.status_stationary.eq("RECOVERED")
          & pairs.status_moving.eq("RECOVERED"))
    for col in ("mu_east_stationary", "mu_north_stationary",
                "mu_east_moving", "mu_north_moving"):
        ok &= np.isfinite(pd.to_numeric(pairs[col], errors="coerce"))
    good = pairs[ok]
    if len(good) < 2:
        raise ValueError("Need two or more paired successes")
    recovered_difference = good[["mu_east_moving", "mu_north_moving"]].to_numpy(dtype=float) - good[
        ["mu_east_stationary", "mu_north_stationary"]].to_numpy(dtype=float)
    resid = recovered_difference - np.asarray(injected_mu, dtype=float)
    size = np.linalg.norm(resid, axis=1)
    return {
        "n_sites": int(len(pairs)),
        "n_both_recovered": int(len(good)),
        "fraction_both_recovered": float(len(good)/len(pairs)),
        "truth_mu_east_north_masyr": np.asarray(injected_mu, dtype=float).tolist(),
        "paired_difference_error_mean_east_north_masyr": np.mean(resid, axis=0).tolist(),
        "paired_difference_error_median_east_north_masyr": np.median(resid, axis=0).tolist(),
        "paired_difference_error_std_east_north_masyr": np.std(resid, axis=0, ddof=1).tolist(),
        "paired_residual_amplitude_p50_p90_p95_masyr": np.quantile(size, [.5, .9, .95]).tolist(),
        "paired_residual_amplitude_max_masyr": float(np.max(size)),
        "scientific_status": "DIAGNOSTIC_CONDITIONAL_ON_BOTH_RECOVERIES_ONLY",
        "limitations": [
            "Shared same-site biases cancel partially, but reference registration errors may remain.",
            "Only originally recovered site pairs are included (selection-dependent).",
            "The injected PSFs are compact-field-source stamps, not validated ultra-red ePSFs.",
            "Not an empirical Gaussian sigma, new actual-source motion, or a trigonometric parallax."
        ]
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=int, choices=[387132,1651050], required=True)
    parser.add_argument("--artifact-zip", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    with zipfile.ZipFile(args.artifact_zip) as z:
        name = f"{args.candidate}/trials.csv"
        trials = pd.read_csv(BytesIO(z.read(name)))
    frozen = pd.read_csv(ROOT / "data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    row = frozen.loc[args.candidate]
    mu = [float(row.mu_alpha_cosdec_masyr), float(row.mu_delta_masyr)]
    result = {"candidate_id": args.candidate,
              "input_artifact_zip": args.artifact_zip.name,
              **summarize(trials, mu)}
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
