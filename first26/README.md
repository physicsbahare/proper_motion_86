# Proper-motion audit: first 26 rows of `galaxy and brown dwarf candidates_new_090630`

This is a separate run inside the same repository. It does **not** overwrite the existing 45-candidate audit.

## Input

`data/candidates_first26_090630.csv` contains the first 26 numbered rows of the supplied workbook. The original SED/model columns and the `Classification_3σ` value are retained for provenance, but **none of them are used to decide proper motion**. The PM analysis uses only the candidate sky position plus public archival imaging.

## Astrometric hierarchy

The shared `pm86` pipeline is reused so this run benefits from the lessons learned during the 45-candidate audit:

1. Prefer a long-baseline **same-filter JWST/NIRCam** epoch pair.
2. If unavailable, allow **cross-filter NIRCam** only as an explicitly labeled lower-evidence case.
3. Do not silently mix HST and JWST. HST coverage is preserved in the archive inventory and can be used later in a dedicated cross-instrument registration test.
4. Verify detector coverage with the **full embedded JWST gWCS**, not approximate header WCS.
5. Measure every retained exposure with forced photometry, empirical blank-aperture noise, DQ-aware masks, 2-D Gaussian and COM centroids, Monte-Carlo centroid errors, and local pixel scale.
6. Register epochs locally using field controls and preserve registration residuals.
7. Save pairwise PM diagnostics and classify conservatively as `MOVING`, `CONSISTENT_WITH_ZERO`, `AMBIGUOUS`, `AMBIGUOUS_SYSTEMATICS`, or `INSUFFICIENT_DATA`.
8. Run the empirical **60 mas method-sensitivity audit** on any usable exposure measurements.

## Speed without throwing information away

The 26 candidates are distributed over **9 independent GitHub-hosted shards** (about 3 candidates per runner), which reduces wall-clock time without changing the per-candidate science. Each candidate still keeps the full evidence directory:

- archive inventory
- selected epoch pair
- MAST product tables
- product coverage audit
- exposure measurements
- field controls
- registered positions and residuals
- pairwise PM table
- diagnostic figure
- complete source input row
- final summary JSON

The aggregation step also has fail-visible QC. In particular, a run cannot be treated as scientifically successful if a `PAIR_FOUND` candidate merely became `INSUFFICIENT_DATA` because every remote FITS attempt failed for an infrastructure/dependency reason. This catches the failure mode discovered in the earlier runs.

## Outputs

Compact outputs are committed to `results_first26/`, while detailed candidate evidence is under `results_first26/candidates/`.

The main files are:

- `ALL26_PM_AUDIT.csv`
- `CLASSIFICATION_COUNTS.csv`
- `REASON_COUNTS.csv`
- `PRODUCT_ACCESS_QC.csv`
- `RUN_QC.md`
- `SUMMARY.md`
- `centroid_tolerance/`
