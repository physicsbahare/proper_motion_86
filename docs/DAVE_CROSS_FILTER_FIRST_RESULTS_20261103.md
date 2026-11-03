# Real MAST cross-filter null: first empirical results (2026-11-03)

**Status: preliminary, NOT a validated correction or revised detection
significance.** Reproduced from GitHub Actions run
[41593637302](https://github.com/physicsbahare/proper_motion_86/actions/runs/41593637302),
artifact `dave-cross-filter-astrometry-audit` (ID 13211494828).
The workflow completed successfully, including synthetic tests and public
MAST cutout retrieval for **both** candidates.

| COSMOS2025 ID | filters | baseline (days) | held-out matches | retained after clipping | median east/north null (mas) | east/north scatter (mas) | illustrative conservative 2D score | frozen nominal 2D score |
| --- | --- | ---: | ---: | ---: | --- | --- | ---: | ---: |
| 387132 | F444W/F410M | 727.4425 | 12 | 10 | -0.85 / +4.02 | 9.98 / 10.85 | 3.49 | 5.60 |
| 1651050 | F444W/F356W | 497.6405 | 9 | 9 | +6.98 / +0.15 | 15.58 / 12.23 | 4.56 | 5.06 |

The null covariance matrices (mas², east/north order) are:

- 387132: `[[99.58, -25.39], [-25.39, 117.70]]`
- 1651050: `[[242.75, 138.65], [138.65, 149.63]]`

**Neither source is established at a revised Gaussian 5σ significance.**
The sensitivity scores above add total held-out residual covariance to the
frozen PM errors and subtract the median null displacement. They are **not**
calibrated detection significances; they may double-count registration noise.
The 2D quadrature nominal scores are likewise not Gaussian z-values.

### Crucial limitations uncovered by the real-data run

1. **Only 10 and 9 controls survived**. The code's diagnostic threshold is
   15; both results explicitly report `sufficient_for_diagnostic=false`.
   Individual controls are not classified as stationary galaxies. Empirical
   null tail probabilities cannot resolve 5σ at these sample sizes.
2. **Color dependence is not calibrated.** 387132 has an apparent
   east-residual versus instrumental flux-ratio correlation
   `r=0.756, p=0.0045` in the **12 pre-clipping controls**. The most extreme
   negative-color residual is -96.7 mas east, -83.7 mas north, and is clipped
   from the covariance estimate. This association must be re-tested with
   robust/leave-one-out statistics and vetted sources; it is NOT evidence
   of a calibrated chromatic PSF shift. 1651050 has `r=-0.105, p=0.789`
   with only 9 sources; absence of a significant correlation is not evidence
   of no chromatic systematic. The 1651050 target's *instrumental* color
   proxy is outside the control range; that proxy is not AB color.
3. **The independent all-exposure refit disagrees with the frozen published
   measurements**: 387132 gives (`mu_east=-19.98`, `mu_north=+2.06`)
   mas/yr and 1651050 gives (`-61.74`, `+17.15`) mas/yr. Both are
   labeled `classification=AMBIGUOUS` by the existing refit pipeline.
   Before adoption, identify why the frozen values and current refit
   differ, and validate centroid consistency; the 387132 COM and Gaussian
   PM solutions differ substantially.
4. The representative-pair null is **not the full epoch-combined target
   estimator**. Multi-exposure correlations and source-level bootstrap
   calibration are needed. A finite stationary-galaxy reference set, matched
   target color and S/N, and PSF/injection tests remain outstanding.
5. Do not modify the published number-density prediction from this test.
   PM confirmation confidence and photometric selection completeness are
   distinct; any inference about the Galactic fraction of F444W dropouts
   needs the separately established PM-capable footprint.

### Next implementation priorities

- Reconcile independent refits against the frozen target measurement,
  exposure-selection and centroid-association decisions.
- Aggregate held-out nulls at the **source/epoch level**, not by treating
  multiple dither observations of the same source as independent.
- Vet stationary galaxies by morphology and multi-epoch consistency;
  remove moving stars, blends and unstable Gaussian fits.
- Use calibrated multi-filter fluxes and/or ePSF color-dependent injection
  experiments to assess ultracool dwarf versus reference SED mismatch.
- Separate registration variance already included in nominal covariance
  from any *additional* filter-dependent variance before recomputing PM
  uncertainty. Calibrate 2D null p-values rather than label raw scores
  'sigma'.

### Safe provisional response to Priyanka (not final manuscript wording)

> Thank you for raising this important point. We have now performed an
> initial empirical cross-filter test using the original JWST detector
> exposures. Leave-one-out astrometric registration of local reference
> sources indicates residual scatters of approximately 10–11 mas per
> coordinate for the F444W/F410M pair and 12–16 mas for F444W/F356W.
> These are preliminary control-source dispersions, not yet isolated
> chromatic PSF systematics. The local samples are small and not yet
> restricted to stationary galaxies, so we are validating the reference
> population and the epoch-level covariance before adopting revised
> proper-motion significances. We will report the validated systematic
> budget and adjust the astrometric interpretation accordingly.

This response must not claim the requested complete validation has
already been achieved.


### Follow-up validation, 2026-11-03

The original color-correlation diagnostic computed Pearson correlations on **all**
held-out controls even when the covariance estimate had clipped extreme
residuals. This is inconsistent and can overstate apparent chromatic evidence,
especially for 387132 (one source had a very large residual). The branch now
computes correlations separately on **all** and **retained** controls, reports
Pearson and Spearman coefficients, and includes leave-one-out Pearson ranges.
A synthetic regression test injects a high-leverage outlier. The new workflow
must finish before interpreting updated color results; do not quote the old
r=0.756 as a validated chromatic detection.

The full scientific validation remains incomplete. The most important
remaining question is the source/epoch association difference between the
frozen PM estimates and the independent all-exposure refits. It must be
resolved before revised candidate PMs or detection classifications are adopted.
