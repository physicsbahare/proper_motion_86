# Outputs

## `results/ALL45_PM_AUDIT.csv`

One row per input candidate. Important fields include:

- `candidate_id`, `ra_deg`, `dec_deg`
- `pair_type`: `SAME_FILTER_JWST` or `CROSS_FILTER_JWST`
- `filter_early`, `filter_late`
- `baseline_days`, `baseline_year`
- `mu_alpha_cosdec_masyr`, `mu_delta_masyr`
- corresponding uncertainty columns
- `significance_2d`
- `pair_direction_fraction`
- `centroid_method_disagreement`
- `cross_filter_color_systematic_flag`
- `classification`
- `reason` when no measurement can be made

## Per-candidate files

`archive_inventory.csv` preserves the public-archive observations found around the catalog coordinate.

`selected_epoch_pair.json` records the epoch pair selected before any target PM is measured.

`product_coverage_audit.csv` records every individual CAL product tested with full gWCS and whether the target was actually on that detector.

`exposure_measurements.csv` contains forced photometry, empirical noise, DQ diagnostics, and the two centroid chains.

`field_controls.csv`, `registration_quality.csv`, and `registration_control_residuals.csv` document the local astrometric frame solution.

`registered_positions.csv` contains the locally registered target position in each usable exposure.

`pairwise_pm.csv` contains all earliest x latest exposure PM estimates. These rows are consistency diagnostics and are not independent measurements.

`diagnostic.png` shows local stamps, S/N, and registered positions for fast human review.
