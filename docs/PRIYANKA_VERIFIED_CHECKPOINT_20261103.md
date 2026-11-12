# Priyanka cross-filter audit: verified real-data checkpoint (2026-11-03)

**Source:** successful GitHub Actions run [41605446965](https://github.com/physicsbahare/proper_motion_86/actions/runs/41605446965), artifact `dave-cross-filter-astrometry-audit` (ID 13218202106), head SHA `acc8483b9475ca385cde1632b05bc27ca440c156`. These values were read from both candidate `summary.json` and `frozen_refit_reconciliation.json`, not inferred from the workflow's green check. This is a **diagnostic checkpoint**, not an adopted manuscript measurement.

| Property | 387132 | 1651050 |
|---|---:|---:|
| Filters | F444W/F410M | F444W/F356W |
| Baseline, days | 727.4425 | 497.6405 |
| Raw matched controls | 12 | 9 |
| Robust retained controls | 10 | 9 |
| Median east residual, mas | -0.854 | +6.980 |
| Median north residual, mas | +4.024 | +0.154 |
| East scatter, mas | 9.979 | 15.580 |
| North scatter, mas | 10.849 | 12.233 |
| Empirical residual covariance (mas²) | [[99.581,-25.390],[-25.390,117.701]] | [[242.751,138.655],[138.655,149.635]] |
| Nominal published 2D quadrature score | 5.60 | 5.06 |
| Conservative *diagnostic* Mahalanobis score | 3.49 | 4.56 |
| Original μ_east ± err, mas/yr | -22.9 ± 4.1 | -55.2303 ± 11.5336 |
| Original μ_north ± err, mas/yr | +2.4 ± 5.8 | +16.4082 ± 9.9982 |
| Independent refit μ_east ± err, mas/yr | -19.983 ± 4.101 | -61.741 ± 15.397 |
| Independent refit μ_north ± err, mas/yr | +2.056 ± 6.320 | +17.150 ± 10.369 |
| Registered exposures early / late | 4 / 2 | 4 / **1** |
| Independent refit classification | AMBIGUOUS | AMBIGUOUS |

## What the newer color robustness test actually says

For **387132**, full-sample east residual versus instrumental-color Pearson r = 0.756, p = 0.00449 (N=12), but retained-sample Pearson r = 0.535, p = 0.111 (N=10), and retained-sample Spearman rho = 0.091, p = 0.803. The retained-sample leave-one-out Pearson r range is **-0.088 to +0.665**. Thus a significant and stable chromatic correlation **has not been established**; the earlier apparent correlation is sensitive to influential observations. For **1651050**, retained east Pearson r = -0.105 (p=0.789) and north Pearson r = -0.458 (p=0.215), N=9. None of these small-N tests excludes a color-dependent offset at the brown-dwarf color.

## New important audit finding

The independent refit for **1651050 has only ONE registered late-epoch exposure** (four early, one late). For **387132**, four early and two late exposures were registered. This asymmetry is critical: a one-exposure late epoch cannot provide an internal late-epoch repeatability test, so a statistically precise cross-filter shift is especially vulnerable to centroid misassociation, PSF variation, and epoch-specific astrometric error. The original manuscript exposure selection/centroid algorithm must be matched exactly before interpreting the independent refit discrepancy as physical. The refit differences are +2.917 and -0.344 mas/yr (east, north) for 387132; -6.511 and +0.742 mas/yr for 1651050. These differences are **not independent** of the frozen solution; no naive significance-of-difference test is valid.

## Decision

**Neither candidate has a publication-ready systematic-corrected Gaussian significance.** The empirical field-control null is useful but (a) N=9–10 is insufficient for a calibrated extreme-tail null; (b) the controls are not individually verified stationary background galaxies; (c) their instrumental flux ratios are not calibrated stellar/galaxy colors and may not represent brown-dwarf PSFs; (d) only a representative detector pair is used for the null, whereas the target uses multiple early exposures; (e) adding the entire null covariance to published registration errors can double-count noise; and (f) the original/refit centroid associations and exposure subsets have not been reconciled. **Do not publish 3.49 or 4.56 as sigma.** The photometric number-density prediction is unchanged; the number of *astrometrically confirmed* objects must not be asserted without the additional validation.

## Specific next tests

1. Trace the original published candidate pair filenames and centroid methods and perform a like-for-like frozen/refit comparison with the **same exposures and weights**. Audit the single late exposure for 1651050 and visually inspect source association/fit residuals.
2. Select stationary galaxies using morphology plus independent astrometric quality checks, cross-match control identities across epochs, and report per-control position, flux ratio, isolation and fit flags. Do not assume stars are stationary.
3. Run a source-level leave-one-out *full epoch* fit for every control, using the exact target estimator; use robust bootstraps of **sources**, not dither pairs, to assess correlated covariance. Separate centroid noise already in published errors from additional systematic components.
4. Test filter-dependent PSF and target-color extrapolation using empirical/ePSF fits and injection/recovery. Only then derive an uncertainty budget and a calibrated null detection probability.

## Provisional English reply to Priyanka (do not send as final result)

> Thank you for highlighting the possible cross-filter astrometric systematics. We have now performed an initial leave-one-out local-reference test on the original JWST CAL exposures, separately for F444W/F410M and F444W/F356W. The preliminary retained-reference scatter is approximately 10–11 mas per coordinate for the former and 12–16 mas for the latter, with only 10 and 9 usable controls, respectively. A full-sample color correlation in the F444W/F410M field is not robust after outlier handling. We are validating stationary reference-galaxy selection and reconciling the independent epoch fits with the original source-centroid measurements before adopting a systematic covariance or revised detection significance. We will report these checks and avoid interpreting the preliminary quadrature sensitivity scores as Gaussian sigma values.

**Manuscript placeholder (not yet a validated result):** “We assessed possible filter-dependent centroid shifts by applying leave-one-out local affine registrations to nearby high-S/N reference sources in the two cross-filter epoch pairs. The preliminary reference residual distributions have coordinate scatters of order 10–16 mas. These include reference centroid noise and potential foreground motions and therefore do not yet isolate a chromatic PSF systematic. We are verifying stationary controls and propagating a full-epoch covariance before assigning a revised astrometric detection significance.”
