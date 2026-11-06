# Original JWST F356W pixel inspection for candidate 1651050

The original September 17, 2026 astrometric audit detected 1651050 in only **one of four covering F356W CAL exposures**. Its accepted F356W 2-D Gaussian and center-of-mass centroids differed by **2.18 pixels (~138 mas)**. A `JUMP_DET` pixel occurred within four pixels of the fixed catalog position (the adopted centroid pixel itself was clean). This result must be inspected before calling the F444W/F356W displacement a confirmed proper motion.

## New diagnostic outputs

The dedicated workflow reads only six source-verified individual JWST Stage-2 CAL products from public MAST: all four covering F356W exposures plus the two astrometrically usable F444W reference exposures, and writes:

- `panels/*_review.png`: SCI, error-weighted S/N, DQ, Gaussian residual, wide-field neighbors and 1-D profiles for each exposure, with catalog position, Gaussian, and COM marked.
- `pixel_stamps/*.npz`: original small SCI, ERR, DQ arrays, source coordinates and offsets, retaining full-detector origin metadata for reproducible independent inspection.
- `centroid_quality_summary.csv`: numerical centroid-method disagreement, DQ flags, measured SNR, residual fitting metrics.
- `MANUAL_REVIEW_CHECKLIST.csv`: blank reviewer signoff fields for all six exposures.
- `READ_ME_FIRST.md`: instructions for independent visual checks.

The workflow does **not** overwrite any published proper-motion measurement or reclassify the source. Its Gaia-like/catalog location cross is a sky coordinate, not an epoch-propagated PM model. Some rejected exposures can display flux without providing usable precision astrometry.

[Review workflow](https://github.com/physicsbahare/proper_motion_86/actions/workflows/visual_review_1651050.yml)

### Manual review is required

Focus on (a) whether the accepted F356W Gaussian centroid lands on the correct faint source, (b) whether the center-of-mass shift arises from a neighbor or structure, and (c) whether nearby masked `JUMP_DET` pixels or registration effects could bias the fitted position. A third independent epoch and an exact JWST L2 parallax model remain advisable.
