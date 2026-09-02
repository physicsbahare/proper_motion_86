# proper_motion_86

Independent public-archive validation of proper motion for the **45 remaining candidates** in the 86-object COSMOS sample.

## Why this repository exists

The 86-object candidate list has been split by independent validation into:

- **40 objects**: independently confirmed as same-filter F444W single-exposure artifacts/cosmic-ray-like events.
- **282040**: independently remeasured as a real source with a long-baseline JWST displacement (`mu_alpha*cos(dec) = -22.9 +/- 4.1 mas/yr`, `mu_delta = +2.4 +/- 5.8 mas/yr`).
- **45 objects in `data/candidates_45.csv`**: the objects audited here.

This repository intentionally does **not** use another team's frozen cutouts, proper-motion measurements, or stationary/moving classifications as inputs. The physical input is the original candidate RA/Dec plus public MAST imaging.

## What the remote pipeline does

GitHub Actions runs the analysis on GitHub-hosted runners; no local computer is required.

For each of the 45 candidates the pipeline:

1. queries public MAST imaging at the catalog coordinate;
2. separates same-night exposures from genuinely independent epochs;
3. prefers a **same-filter JWST/NIRCam** long baseline, with F444W first;
4. only if no same-filter baseline exists, considers a JWST cross-filter pair;
5. retrieves only individual Stage-2 `*_cal.fits` products that are needed;
6. uses embedded distortion-aware JWST **gWCS** to determine actual detector coverage;
7. reads small `SCI`, `ERR`, and `DQ` cutouts by HTTP byte range rather than downloading whole CAL files;
8. performs forced photometry, empirical local-noise checks, and DQ diagnostics;
9. measures the target with both 2-D Gaussian and center-of-mass centroid chains;
10. detects local field controls and robustly registers each exposure to a reference exposure;
11. fits proper motion only when the target is detected/centroidable in two independent epochs;
12. saves exposure measurements, control residuals, registration quality, pairwise PM diagnostics, plots, and a final classification.

## Scientific rules built into the code

- A non-detection in a different filter is **not** evidence that a source is an artifact or a mover.
- Same-night dithers are not independent PM epochs.
- Same-filter baselines have the highest evidence grade.
- Cross-filter PM is allowed only when the target is detected in both epochs and is accompanied by local color/residual and centroid-method diagnostics.
- DQ `JUMP_DET` / `OUTLIER` pixels are diagnostics for artifact-sensitive photometry, but they are masked for precision astrometry.
- The eight (or more) exposure pairings of two epochs are consistency diagnostics, not statistically independent PM measurements.
- `CONSISTENT_WITH_ZERO` means no significant motion at the achieved precision; it is not a proof that the true PM is exactly zero.
- If the data cannot support a clean fit, the answer is `INSUFFICIENT_DATA`, not `STATIONARY` by assumption.

## Result classes

- `MOVING`
- `CONSISTENT_WITH_ZERO`
- `AMBIGUOUS`
- `AMBIGUOUS_SYSTEMATICS`
- `INSUFFICIENT_DATA`

## Repository layout

```text
.
├── data/
│   ├── candidates_86_input.csv
│   ├── candidates_45.csv
│   ├── validated_artifact_ids.txt
│   └── PROVENANCE.md
├── docs/
│   ├── METHOD.md
│   └── OUTPUTS.md
├── scripts/
│   ├── run_shard.py
│   └── aggregate_results.py
├── src/pm86/
│   ├── archive.py
│   ├── astrometry.py
│   ├── config.py
│   ├── measurement.py
│   └── pipeline.py
├── tests/
└── .github/workflows/run_45.yml
```

## Running remotely

The workflow runs automatically when the pipeline/data/workflow is changed, and can also be launched from **Actions -> Independent PM audit of 45 candidates -> Run workflow**.

The 45 candidates are split across five parallel jobs. A final aggregation job collects all shard outputs and commits the compact scientific results back under `results/`. Full JWST FITS files and local NPZ cutouts are deliberately not committed.

## Main outputs

After a successful run:

```text
results/
├── ALL45_PM_AUDIT.csv
├── CLASSIFICATION_COUNTS.csv
├── SUMMARY.md
└── candidates/
    └── candidate_<ID>/
        ├── archive_inventory.csv
        ├── selected_epoch_pair.json
        ├── products_early.csv
        ├── products_late.csv
        ├── product_coverage_audit.csv
        ├── exposure_measurements.csv
        ├── field_controls.csv
        ├── registered_positions.csv
        ├── registration_quality.csv
        ├── registration_control_residuals.csv
        ├── pairwise_pm.csv
        ├── diagnostic.png
        └── summary.json
```

See `docs/METHOD.md` for the detailed methodology and `docs/OUTPUTS.md` for column/interpretation notes.
