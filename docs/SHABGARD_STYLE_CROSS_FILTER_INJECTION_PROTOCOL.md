# Shabgard-style artificial-source injection and recovery for Dave's two PM sources

## Provenance and difference from Shabgard

The Shabgard (COSMOS2025 ID 639188) manuscript used **500 artificial sources per epoch pair** with an empirical representation of the target PSF and comparable brightness, injected into clean regions of the genuine JWST images. It tested both stationary fakes and fakes displaced by the measured motion. Main *same-filter* PRIMER/COSMOS-3D F115W, F200W and F356W pairs yielded stationary-source 95th-percentile false PM amplitudes about **3.2, 5.0, 27.1 mas/yr**, versus measured **88.62, 88.89, 90.86 mas/yr**, respectively. In the Shabgard analysis these percentiles characterized detection robustness rather than calibrating a formal 5-sigma tail. The older correction notebook `proper_motion_fig4_corrected_baseline.ipynb` and `proper_motion_baseline_fig4_patch.py` reference `ID639188_fake_source_PM_test/ID639188_fake_source_trials.csv`, but the **original source injection notebook/code has not been verified in this GitHub repository**. Do not claim code reuse without source inspection.

Shabgard's principal epoch pairs used the *same* filter at both epochs; Dave's pairs use different filters, making chromaticity, PSF variation, different S/N and source association central rather than negligible.

## Concrete proposed test, using existing MAST CAL infrastructure

Targets and *frozen* epoch pairs:

- COSMOS2025 387132 (legacy 282040): early F444W / late F410M, baseline 727.4425 days, nominal two-coordinate quadrature score 5.60 (not a 1D Gaussian sigma).
- COSMOS2025 1651050: early F444W / late F356W, baseline 497.6405 days, nominal two-coordinate quadrature score 5.06.

1. Retrieve the **same original detector CAL exposures** actually used in each target result from MAST, using the existing `pm86.archive` bounded remote cutout loader. Save filenames, filters, exposure times, actual WCS mapping and DQ. Compare against the frozen paper exposure selection, particularly the single accepted late F356W image for 1651050. Do not change the science result while fitting mocks.
2. For each filter/epoch create a *filter-specific, spatially representative PSF*. Validate against real point sources. For 1651050 the injected F356W source must reproduce its actual much fainter F356W flux and S/N, while F444W must match F444W. A constant-amplitude, same-PSF fake across both filters is scientifically misleading. Consider source-color-dependent effective PSFs/ePSF and subpixel centroid phase; mark unverifiable approximations.
3. Preselect clean sky coordinates on the **common tangent plane**, before any injection and independent of the recovery outcome. Require sky positions map inside every chosen exposure, respect DQ masks and avoid catalogued real sources and candidate. Pick a range of detector locations and background/crowding environments. Seed random trials for reproducibility; target at least 500 *independent sky positions* per object/pair, if coverage permits. Treat repeated dithers of the same sky position as correlated, not independent samples.
4. Inject PSF sources *into genuine SCI cutouts*, accounting for the exposure's pixel sampling, local background, flux/ERR and DQ, using sky-to-detector gWCS for every image. Generate the source on each native detector pixel grid; do **not** overwrite original CAL products. Document whether injected source Poisson noise is represented by the ERR model. Do not reuse the identical PSF as a perfect recovery template without an independent mismatched-PSF check.
5. Run exactly the target pipeline: source association, full target centroid method/forced fit, quality/SNR gate, local registration against real controls, same exposure weighting and epoch combination, resulting 2D PM. Exclude the injected source from reference-control selection. For each trial store success/rejection and *why*. Do not measure only recovered centroids using a simpler method if the paper method differs.
6. Two hypotheses in the SAME selected positions and noise environment:
   - **Null:** source is fixed at one sky position at both epochs. Record the recovered `mu_east`, `mu_north`, amplitude, target-statistic score and false detection rate (including method failures).
   - **Moving:** source is displaced by the observed 2D vector times the actual per-exposure baseline. Record recovery bias, dispersion, containment and completeness. If needed consider parallax **separately** from the injected proper-motion displacement.
7. Include color/PSF stress cases: real source's measured flux ratio, redder/bluer plausible SEDs, PSF width/shape perturbations, subpixel phase and local environment. Compare null and moving trials. Reproduce the Shabgard **95th-percentile stationary false-motion amplitude** and ratio of observed amplitude to that limit, but do not call this 95% result a 5-sigma confirmation.
8. Track errors that synthetic sources **cannot** constrain: survey-common WCS/registration biases, true chromatic PSF terms absent from the assumed injection PSFs, selection of foreground moving control stars, common epoch registration covariance and parallax. Use independent stationary galaxies, cross-validated control registration, and inject after/before re-solving registration as appropriate. Need a systematic uncertainty that is incremental, NOT the full control scatter added twice.

## Required reporting

For each source produce a trial-level CSV with sky position, filter fluxes/SNR, type `stationary|moving`, true and recovered 2D displacement/PM, original and recovered estimator covariance, recovered/rejected flag, exposure IDs, method and PSF case. Aggregate:

- 500 trial target (report actual N), accepted fractions, centroiding failure and faint-filter dropout frequencies.
- Stationary recovered PM amplitude p50/p68/p95/p99, and 2D axis bias/covariance.
- Moving recovery bias and 68% interval for both axes, with red/PSF stress-case range.
- Comparison of observed motion amplitude against stationary p95, specifically at each candidate's **actual brightness, filter ratio, and epoch selection**.
- Empirical null tail fraction for the observed *same* 2D statistic, with explicit Monte Carlo resolution: 500 trials only resolve p down to approximately 1/(500+1); **cannot validate 5-sigma Gaussian tails**.
- Clear separation of observational proper-motion evidence from photometric BD classification and number-density forecast.

## Scientific go/no-go

**Good diagnostic:** injecting a fixed red source rarely yields an apparent displacement comparable to the observed candidate, while a simulated moving red source is recovered without material bias and retains adequate late-filter detection completeness.

**Failure or inconclusive:** systematic false motion near the observed amplitude, color-sensitive bias, extensive false association, or insufficient S/N in late F356W. Treat detection as tentative and revise wording.

**Status 2026-11-04:** This is an implementation-ready **protocol, not results**. The present PR's empirical field-control null has 10 retained reference sources for 387132 and 9 for 1651050, with ~10–16 mas coordinate scatters. Actual cross-filter injection/recovery has **not yet run**; its results must not be inferred from the Shabgard thresholds. The original Shabgard injection code has not yet been independently retrieved, although the manuscript and the correction notebook document its design and output paths.
