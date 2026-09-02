# Testing the verifier's 60 mas centroid tolerance

The question is not whether 60 mas sounds reasonable. The useful question is whether a submission can pass or fail only because it used a different defensible astrometric method.

The remote pipeline therefore performs a method-sensitivity audit after the candidate analysis finishes. For every clean target measurement with S/N >= 3 it measures:

1. **2DG vs COM centroid separation** on the same original detector exposure.
2. **Raw gWCS vs locally registered position**, to quantify astrometric-frame sensitivity.

The independent six-exposure 282040 measurement is also stored in `data/reference_282040_centroids.csv` and is included as a calibration/reference sample.

## What the 282040 pre-check already shows

For 282040, the six 2DG-vs-COM same-exposure separations are approximately:

```text
14.6, 50.1, 15.8, 7.8, 50.4, 83.0 mas
```

All six exposures have S/N > 8. At a hard 60 mas threshold, 5/6 pass but one legitimate centroid-method comparison differs by about 83 mas. The 95th percentile is about 74.9 mas. Therefore 60 mas is **not demonstrably method-neutral** even in the already-validated reference mover; there is direct evidence that method choice can matter.

This does not mean that 60 mas must be discarded. It means the threshold should either be calibrated across methods, allow a method-aware uncertainty/tolerance, or be supported by an additional PSF/ePSF comparison.

## Remote outputs

The audit writes:

```text
results/centroid_tolerance/CENTROID_TOLERANCE_REPORT.md
results/centroid_tolerance/CENTROID_TOLERANCE_DETAIL.csv
results/centroid_tolerance/CENTROID_TOLERANCE_SUMMARY.csv
results/centroid_tolerance/CENTROID_TOLERANCE_BY_FILTER.csv
results/centroid_tolerance/CENTROID_TOLERANCE_BY_SOURCE.csv
results/centroid_tolerance/CENTROID_TOLERANCE_SWEEP.csv
```

The tolerance sweep evaluates 20, 30, 40, 50, 60, 80, and 100 mas rather than tuning the answer to one pre-selected value.

## Limitation

This directly compares the two centroid chains already retained by the analysis (2-D Gaussian and center of mass), plus local-registration frame effects. It is not a dedicated empirical/ePSF fit. Because the 282040 pre-check already reaches and exceeds 60 mas, a true PSF/ePSF comparison is the strongest remaining test before calling 60 mas method-independent.
