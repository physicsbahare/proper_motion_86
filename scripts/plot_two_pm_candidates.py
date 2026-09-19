#!/usr/bin/env python3
"""Publication-style proper-motion vector plot for the two moving candidates
in the 125-source Step3 sample.

Run from the repository root:
    python scripts/plot_two_pm_candidates.py

Outputs:
    figures/two_pm_candidates_pm_plane.png
    figures/two_pm_candidates_pm_plane.pdf
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DATA = Path("data/two_pm_candidates_125.csv")
OUTDIR = Path("figures")


def main():
    df = pd.read_csv(DATA)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.4, 5.4))

    markers = {
        "validated_mover": "o",
        "strong_pm_candidate": "s",
    }

    for _, row in df.iterrows():
        label_id = str(int(row["sample_candidate_id"]))
        if pd.notna(row["legacy_candidate_id"]):
            label_id += f" (legacy {int(row['legacy_candidate_id'])})"

        status_text = (
            "validated mover"
            if row["status"] == "validated_mover"
            else "strong PM candidate"
        )
        legend_label = f"{label_id}: {status_text}"

        ax.errorbar(
            row["mu_alpha_cosdec_masyr"],
            row["mu_delta_masyr"],
            xerr=row["sigma_mu_alpha_cosdec_masyr"],
            yerr=row["sigma_mu_delta_masyr"],
            fmt=markers[row["status"]],
            markersize=7,
            capsize=3,
            linewidth=1.2,
            label=legend_label,
        )

        ax.annotate(
            f"{int(row['sample_candidate_id'])}\n"
            f"{row['filter_early']}-{row['filter_late']}, "
            f"{row['significance_2d']:.1f}σ",
            (row["mu_alpha_cosdec_masyr"], row["mu_delta_masyr"]),
            xytext=(7, 7),
            textcoords="offset points",
            fontsize=9,
        )

    ax.axhline(0, linestyle="--", linewidth=0.8)
    ax.axvline(0, linestyle="--", linewidth=0.8)

    ax.set_xlabel(r"$\mu_{\alpha*}$ (mas yr$^{-1}$)")
    ax.set_ylabel(r"$\mu_{\delta}$ (mas yr$^{-1}$)")
    ax.set_title("Proper-motion candidates in the 125-source sample")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    fig.tight_layout()

    png = OUTDIR / "two_pm_candidates_pm_plane.png"
    pdf = OUTDIR / "two_pm_candidates_pm_plane.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {png}")
    print(f"Wrote {pdf}")


if __name__ == "__main__":
    main()
