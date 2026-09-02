# Testing the verifier's 60 mas centroid tolerance

The question is not whether 60 mas sounds reasonable. The useful question is whether a submission can pass or fail only because it used a different defensible astrometric method.

The remote pipeline therefore performs a method-sensitivity audit after the candidate analysis finishes.

For every clean target measurement with S/N >= 3 it measures:

1. **2DG vs COM centroid separation.** This directly asks how far two independently retained centroid algorithms place the same source on the same original detector exposure.
2. **Raw gWCS vs locally registered position.** This asks how much a valid local affine registration changes the target's astrometric coordinate relative to the unregistered gWCS solution.

The audit writes a row-level table, filter-level summaries, and a tolerance sweep at 20, 30, 40, 50, 60, 80, and 100 mas. A source with `method_sensitive_at_tolerance=True` is a concrete case where a 60 mas verifier could depend on the centroid algorithm. A source with `registration_sensitive_at_tolerance=True` shows analogous frame/registration sensitivity.

The headline report is saved as:

```text
results/centroid_tolerance/CENTROID_TOLERANCE_REPORT.md
```

with detailed evidence in:

```text
CENTROID_TOLERANCE_DETAIL.csv
CENTROID_TOLERANCE_SUMMARY.csv
CENTROID_TOLERANCE_BY_FILTER.csv
CENTROID_TOLERANCE_SWEEP.csv
```

## Important limitation

This is a direct comparison of the two centroid chains already retained by the pipeline (2-D Gaussian and center of mass), plus local-registration frame effects. It is not a substitute for a dedicated empirical/ePSF fit. If the observed method spread is comfortably below 60 mas, that is evidence against strong method dependence in this sample. If it approaches or exceeds 60 mas, a true PSF/ePSF comparison should be added before approving a hard threshold.
