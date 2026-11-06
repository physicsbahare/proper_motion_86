# Step3 42-candidate proper-motion audit

Method: epoch-cached exhaustive DQ-aware JWST/NIRCam audit.

Candidates completed: **42/42**

## Classification counts

- `INSUFFICIENT_DATA`: 38
- `AMBIGUOUS_SYSTEMATICS`: 2
- `AMBIGUOUS`: 1
- `CONSISTENT_WITH_ZERO`: 0
- `MOVING`: 1

## PM status counts

- `INSUFFICIENT_DATA`: 38
- `PM_MEASURED`: 4

These counts match `ALL42_PM_SUMMARY.csv` and `CLASSIFICATION_COUNTS.csv`
after the conservative systematic-precedence correction (candidate 1038406).
The older text summary had become stale and incorrectly reported one
`CONSISTENT_WITH_ZERO` object. No measurements or classifications are
altered by this documentation correction.

The one `MOVING` classification (1651050) uses a **cross-filter**
F444W/F356W epoch pair. Its significance is the magnitude of the 2-D
displacement in units of coordinate uncertainty; it is not a
parallax-corrected 1-D Gaussian detection statistic.
