# Centroid tolerance audit

Proposed verifier tolerance: **60.0 mas**.

This is an empirical method-sensitivity check, not an assumption that 60 mas is correct. It compares 2-D Gaussian and center-of-mass centroids measured on the same detector image and separately measures local-registration frame shifts. Only clean S/N >= 3 measurements enter the headline statistics.

## Independent 282040 reference pre-check

The previously completed 282040 run contributes 6 S/N>=3 exposures. At 60 mas, 1 exposure(s) have a 2DG-COM separation larger than the threshold; the maximum separation is 83.0 mas and the 95th percentile is 74.9 mas.

## Headline numbers including available reference/calibration rows

- `eligible_snr3_rows`: 6
- `eligible_candidates`: 1
- `method_comparable_rows`: 6
- `method_rows_over_tolerance`: 1
- `method_fraction_within_tolerance`: 0.833
- `method_sep_p95_mas`: 74.882
- `method_sep_p99_mas`: 81.399
- `method_sep_max_mas`: 83.029
- `registration_comparable_rows`: 0
- `registration_rows_over_tolerance`: 0
- `registration_fraction_within_tolerance`: n/a
- `registration_shift_p95_mas`: n/a
- `registration_shift_p99_mas`: n/a
- `registration_shift_max_mas`: n/a
- `assessment`: 60MAS_METHOD_OR_FRAME_SENSITIVITY_DETECTED

## Interpretation

A non-zero `method_rows_over_tolerance` is direct evidence that the same real S/N>=3 source can differ by more than the proposed threshold solely because a different defensible centroid method is used. A non-zero `registration_rows_over_tolerance` indicates analogous astrometric-frame sensitivity. In either case, a hard pass/fail radius should be made method-aware or calibrated with a broader method set rather than treated as universally neutral.

This audit directly tests 2DG versus COM and raw-gWCS versus local affine registration. It does not replace a dedicated empirical/ePSF comparison. Because the 282040 reference already reaches the vicinity of 60 mas, an ePSF/PSF-fit robustness comparison would be the strongest final check before approving 60 mas as a method-independent hard cutoff.
