# Cross-filter centroiding systematics for Dave's two moving candidates

## Purpose and frozen measurement pairs

This is a **preliminary empirical field-control null** addressing Priyanka's
request to quantify cross-filter centroiding systematics, not an automatic
amendment to the paper's proper-motion results.

The runner uses public MAST Stage-2 NIRCam detector exposures (SCI/ERR/DQ,
embedded gWCS, local cutouts only). Candidate coordinates and the published
measurements are taken from **frozen repository inputs**. It selects the
observation epochs whose filter combination and temporal separation match:

| COSMOS2025 ID | Original filters | Baseline (days) | Nominal PM score |
| --- | --- | ---: | ---: |
| 387132 (legacy 282040) | F444W → F410M | 727.4425 | 5.60 |
| 1651050 | F444W → F356W | 497.6405 | 5.06 |

The paper calls these nominal scores 5.60 and 5.06 sigma. The original code
defines these as the quadrature of the two coordinate S/N values, which is **not**
a properly calibrated one-dimensional Gaussian sigma statistic.

## What is actually measured

1. Query the MAST inventory and find the *specified* filter/epoch pair within
   5 days of the frozen baseline; fail closed rather than silently switch to
   unrelated observations.
2. Load small local Stage-2 CAL cutouts, preserve archive coverage errors,
   extract DQ-aware catalog source and reference-source centroids, and attempt
   an independent two-epoch target refit for comparison only.
3. Identify matched, S/N >= 10, isolated local reference sources in a
   representative target-centroidable early/late exposure pair.
4. For **each** matched reference source, leave it **out** of a local robust
   affine registration fitted to the *other* references. Record its apparent
   astrometric displacement in east/north (mas). Reject fewer-than-six-control
   fits. One source is never allowed to help fit its own null residual.
5. Estimate the robust median residual vector, 2x2 covariance, outliers and
   instrumental flux-ratio/color-proxy correlations. The raw flux ratio is NOT
   a calibrated filter color.
6. Compute a **conservative scenario** using the frozen published candidate
   vector/error covariance plus the control residual covariance. This
   illustrative quadrature *can double-count the registration noise* already
   present in the original reported errors. It is NOT a revised paper
   significance or a validated PSF correction.
7. Preserve all data and failures as workflow artifacts. The job cannot silently
   promote a source to a secure moving object.

### Run

**GitHub Actions → Dave cross-filter astrometric systematics audit → Run workflow**
on the feature branch, or push a change matching the workflow path filters.

Locally, if MAST is accessible:

```bash
pip install -r requirements.txt
pip install -e .
pytest -q tests/test_cross_filter_null.py
python -u scripts/run_cross_filter_null.py --out results_cross_filter_null
```

Files produced for each candidate:

- `product_coverage_audit.csv`: exactly which original detector files loaded
- `exposure_measurements.csv`: target centroids, S/N, DQ, filter and MJD
- `registered_target_positions.csv` and `independent_target_refit.json`
  when a defensible refit is possible
- `leave_one_out_residuals.csv`: held-out source offsets, color proxies and S/N
- `summary.json`: per-candidate diagnostics and conservative sensitivity

The root `summary_table.csv` is a quick comparison of both candidates.

## Caveats, important for a defensible response to Priyanka

**This preliminary run is not yet the full requested calibration.**

- A catalog of automatically detected controls contains a mixture of
  stars and galaxies. Proper motions of foreground stars can broaden the
  null, and control objects are not yet morphologically classified.
- Control fluxes are segmentation measurements on detector cutouts, not
  calibrated AB photometry. The very red ultracool-dwarf spectra may fall
  outside the field-source color range. Do not extrapolate a chromatic PSF
  correction to either target without supporting data.
- A single representative exposure pair is not equivalent to the full
  multi-exposure/epoch target fit. A final null must use the **same combination
  of registered exposures and correlation model** as the actual candidate fit.
- Empirical control scatter includes registration and reference-centroid
  noise. Adding it wholesale to the current error budget may **double count**.
  A defensible final covariance requires source bootstraps/injection and a
  variance-component analysis.
- With N reference controls, empirical tail resolution is only about
  1/(N+1), so a few controls cannot establish 5-sigma significance.
- An F444W/F356W systematic is not interchangeable with the F444W/F410M
  systematic. Both objects must be tested separately.
- A two-epoch astrometric displacement is not a parallax-subtracted
  intrinsic proper motion.

### Decision table

| Outcome | Interpretation / manuscript action |
| --- | --- |
| Null cannot be measured (archive failure, insufficient controls) | Keep the nominal measurements explicitly **provisional**, report the limitation, do not claim independently validated 5 sigma |
| Null scatter small, enough sources, no color trend | Encouraging; still verify galaxy stationarity, target-color compatibility, and multi-exposure covariance before adopting a revised significance |
| Null broad or strong color dependence | Re-estimate systematics using appropriate source colors, PSF fits, and control bootstraps; revise uncertainty and detection category if needed |
| Systematic-corrected motion no longer significant | Retain measured displacement but classify as **tentative**, do not state astrometric confirmation; SED classification remains separate |

A successful GitHub run is **not** by itself enough to label the tests
publication-ready. The next phase, if the cutouts can be obtained, is to
validate background-galaxy selection, perform morphology/PSF tests, and
propagate the empirical per-epoch covariance into a true source-level null.
No number-density prediction is changed by this audit; PM completeness and
photometric selection are distinct.
