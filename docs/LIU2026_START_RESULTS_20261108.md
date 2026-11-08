# Execution checkpoint: Liu 2026 bounded astrometric checks (8 Nov 2026)

**Branch**: `plan/liu2026-minimal-astrometry-followups-20261108`.  
**Science scope**: COSMOS2025 IDs 387132 and 1651050 only.  
**Base**: `review/cross-filter-systematics-20261103`.  
**Status:** Previous-data audit completed, one *new* small detector-level fit-window sensitivity run started; final newly measured fit-window results are not yet claimed in this checkpoint.

## 1. Data availability and inputs

- Current steps use **public MAST Stage-2 NIRCam CAL original detector exposures**; no user-local computer data are required. The previous local COSMOS2025/LePHARE crossmatch already yielded [17 photometric galaxy candidates](PRIYANKA_STATIONARY_GALAXY_CONTROL_MATCH.md) (9 by 387132, 8 by 1651050), committed as `data/field_controls_priyanka/cosmos2025_photometric_galaxy_controls_17.csv`. Those are **not spectroscopically verified** stationary galaxies.
- Original frozen two-epoch pairs: **387132 F444W/F410M** (727.4425 days); **1651050 F444W/F356W** (497.6405 days).
- Original target-centroidable CALs for multi-exposure refits: **4+2** for 387132 and **4+1** for 1651050.
- No third independently validated epoch exists for 1651050 in the archived HST/MIRI searches; a measured trigonometric parallax is **not available** for either candidate.

## 2. Verified existing full-epoch injection results, NOT a new rerun

**Exact source:** [real-CAL full-epoch run 41814514100](https://github.com/physicsbahare/proper_motion_86/actions/runs/41814514100), including successful jobs 126521779534 (387132) and 126521777933 (1651050) and artifact ZIPs 13314054065 / 13314805087. Values below were verified directly from their job logs (`AUDIT_SUMMARY` / terminal `result.json` printout), not inferred from a green status.

| Metric | 387132 | 1651050 |
| --- | ---: | ---: |
| Independent injected sky locations | 32 | 32 |
| Stationary sources recovered / attempted | 30/32 = 93.75% | 24/32 = 75% |
| Known-moving sources recovered / attempted | 31/32 = 96.88% | 24/32 = 75% |
| Stationary recovered PM amplitude ≥ observed candidate PM | 0/30 | 3/24 |
| **Diagnostic empirical exceedance**, (k+1)/(n+1) | **0.03226** | **0.1600** |
| Stationary-source 95th-percentile recovered PM amplitude (mas/yr) | **20.11** | **63.41** |
| Stationary east/north mean recovered PM (mas/yr) | (-5.80, +11.99) | (+6.92, -13.43) |
| Moving recovery east/north bias (mas/yr) | (-6.22, +10.91) | (+7.52, -21.68) |

**Interpretation and limitations:** 387132 is more robust than 1651050 under this specific diagnostic. The source-site sample is small; the empirical exceedance is conditional on successful recovery, and the injected PSF was derived from local compact sources **without a validated ultra-red source SED**. This does **not** calibrate a 5-sigma significance, a full PSF systematic, a revised number-density estimate or the intrinsic PM/parallax decomposition. The 1651050 late epoch still has only one accepted faint F356W detector image. Do not conflate 32 independent sky positions with multiple dithered exposures.

## 3. Verified full-epoch held-out galaxy results from prior code

[Original mixed-reference frame + held-out galaxies](PRIYANKA_ORIGINAL_FRAME_GALAXY_RESULTS_20261106.md):

| Metric | 387132 | 1651050 |
| --- | ---: | ---: |
| Full-CAL independently held-out photometric-galaxy controls | 6 | 7 |
| Median null displacement east/north (mas) | (-0.92, +2.22) | (-1.76, +0.16) |
| Empirical control scatter east/north (mas) | (10.98, 8.80) | (18.82, 15.35) |
| Full-epoch target independent PM east/north (mas/yr) | (-19.98, +2.06) | (-61.74, +17.15) |
| Frozen reported PM east/north (mas/yr) | (-22.90, +2.40) | (-55.23, +16.41) |
| Refit algorithmic classification | AMBIGUOUS | AMBIGUOUS |

Both displacements persist under this reference-frame selection; the *significances are not recalibrated*. Seven small-N galaxy controls cannot prove stationary-tail p-values at 5 sigma. The galaxy-only alternative has NO fully-complete six-CAL held-out null for 387132; do not substitute its five-CAL null covariance.

## 4. Live new step: limited fit-window sensitivity using the original CAL frames

**Source:** `scripts/liu2026_real_cal_window_audit.py`, triggered through `.github/workflows/liu2026_real_cal_window_audit.yml` on this plan branch; [run 42444259609](https://github.com/physicsbahare/proper_motion_86/actions/runs/42444259609).

Per candidate the new diagnostic:
1. Reads existing original frozen filter/epoch pair via the original MAST retrieval functions.
2. Retains the same **target usable** CALs as the original recentered fitting gate and writes all original rejected exposures to `original_target_quality.csv`.
3. Compares 2D-Gaussian and center-of-mass centroids on detector cutouts across fitting-window half-sizes **4,5,7,9 pixels**, using the same DQ masking and raw detector WCS.
4. Writes per-exposure `window_centroids.csv`, `original_target_quality.csv`, `summary.json` and `coverage.csv` as 30-day workflow artifacts, even if one candidate fails.

**Important:** this is a target-specific methodological sensitivity test, **NOT** a color-weighted hybrid PSF, independent Gaia-galaxy reference solution, final proper-motion refit, or measurement of parallax. It can flag centroid instability without establishing the cause. Report only observed diagnostics once jobs finish.

## 5. What to do only if the new test is informative

- Inspect **1651050's single F356W** image against the existing PR #1 ePSF and SCI/ERR/DQ diagnostics. Decide whether weak source association or fit-window dependence makes its displacement fragile; do **not** claim to fix it through more synthetic samples alone.
- If, and only if, this gives a scientifically validated reason and an independent adequate PSF, consider *one* small targeted color-aware centroid/injection comparison using existing code; otherwise stop and document the limit.
- No local FITS is needed for the above. User-local COSMOS2025 images/cat become relevant **only** if needing to inspect additional local mosaics or new locally held astrometric control data not in public MAST or committed catalog summaries.

**Do not automatically merge or modify** paper values, confirmed-object counts, or number-density forecasts.

## 6. Additional NEW paired-site analysis of EXISTING real-CAL injection data (8 Nov)

I downloaded and inspected the actual full-epoch ZIP artifacts (not just workflow logs). At each trial ID, I paired its recovered **moving** fake with its recovered **stationary** fake from the **same sky location**. The residual diagnostic is:

`(recovered moving PM - recovered stationary PM) - known injected PM`

This can cancel some site-specific bias without claiming that the pipeline's full cross-filter or WCS errors have disappeared. Only sites where **both** trials were recovered enter the summary.

| Metric | 387132 | 1651050 |
| --- | ---: | ---: |
| Both stationary + moving recovered at the same site | **30/32** | **23/32** |
| Paired residual mean, east/north (mas/yr) | (-0.14, -1.13) | (+1.16, -9.93) |
| Paired residual coordinate standard deviations (mas/yr) | (4.69, 3.76) | (37.37, 28.55) |
| **Paired residual 2D amplitude median (mas/yr)** | **2.34** | **25.61** |
| Paired residual amplitude 90th percentile (mas/yr) | 4.05 | 51.52 |
| Paired residual amplitude 95th percentile (mas/yr) | 12.30 | 67.64 |

**Interpretation:** The injection pipeline recovers the *difference* between known moving and stationary fakes much more consistently for 387132 than for 1651050. This makes the 1651050 single-late-exposure/cross-filter sensitivity concern more concrete. **Not a new independent observed PM measurement**; same small injected-site sample, recovery conditioning, template PSF and noise limitations remain. Do not treat the numbers as statistical significance or convert to Gaussian sigma.

**Reproduction:** `scripts/summarize_paired_injection_artifact.py --candidate 387132 --artifact-zip github-actions-artifact-13314054065.zip` and `--candidate 1651050 --artifact-zip github-actions-artifact-13314805087.zip`. The ZIP files are existing [run 41814514100](https://github.com/physicsbahare/proper_motion_86/actions/runs/41814514100) downloads; their CSVs are not asserted to live on `main`.

## 7. New original-CAL fitting-window audit: verified partial result

Workflow [run 42444259609](https://github.com/physicsbahare/proper_motion_86/actions/runs/42444259609).

**387132 job completed successfully** (artifact `liu2026-target-centroid-window-387132`, ID 13602946807). All **six** usable original CAL exposures passed the 2D Gaussian fit at four tested window half-sizes (4,5,7,9 pixels). The within-exposure Gaussian centroid span across the tested windows was **at most 0.0156 mas east and 0.0139 mas north** (the maxima occur on one early F444W frame). These tiny differences are **numerical stability of this particular target-anchored Gaussian fit**, NOT astrometric error bars, independent registration, or independent validation of the real PSF. The center-of-mass centroids were much more window/background sensitive and could differ from the adopted Gaussian centroids by **several to ~83 mas** in the per-frame comparison. COM is not a validated substitute for the chosen Gaussian fit.

**1651050:** the separate real-CAL fitting-window job was still being processed when this partial checkpoint was written. No new centroid-window result is claimed for it here. The previous original-pixel audit remains applicable: only one late F356W CAL target has an accepted forced Gaussian centroid, while other covering exposures fail that astrometric gate.

No manuscript tables, significance thresholds, or source classifications have been changed.


## 8. Verified archived original-refit provenance (8 November)

Independently inspected artifact 13218202106 from [run 41605446965](https://github.com/physicsbahare/proper_motion_86/actions/runs/41605446965), including both frozen_refit_reconciliation.json files and full per-exposure diagnostics.

- **387132:** Four accepted F444W frames plus two accepted F410M frames ALL use catalog_position_2dg centroiding. The independent-minus-frozen PM difference is **(+2.91696, -0.34358) mas/yr** (east,north). The fitted baseline differs by just +0.0000157 day, so the baseline is not the explanation.
- **1651050:** Four accepted early F444W frames use catalog_position_2dg, but the **ONLY accepted late F356W frame** (jw05893016008_05101_00001_nrcblong_cal.fits) uses **forced_gaussian**, S/N 5.053. This is a filter change AND centroid-estimator change between epochs. The independent-minus-frozen PM difference is **(-6.51065, +0.74229) mas/yr**; the baseline changes by only -0.000390 day. The other three F356W covering frames fail the accepted finite-centroid gate despite diagnostic target S/N **2.255, 5.086, 3.220**, respectively; S/N alone is not an astrometric centroid.
- The archived original-versus-refit provenance reconciliation remains marked UNRESOLVED. We identified a plausible source of method sensitivity for 1651050, **not a proven cause of the PM discrepancy**.

**Scientific priority:** vet the only late F356W forced-centroid association/PSF and DQ conditions; do not silently include rejected frames or claim high-sigma intrinsic PM. No published classification or numerical result was changed.


## 9. COMPLETED result from both original-CAL fit-window workflow jobs (8 November 2026)

The [workflow run 42444259609](https://github.com/physicsbahare/proper_motion_86/actions/runs/42444259609) is now **completed: success**, both target jobs **success**, with archived per-exposure output artifacts:
- 387132: `liu2026-target-centroid-window-387132` (artifact 13602946807)
- 1651050: `liu2026-target-centroid-window-1651050` (artifact 13603012339)

The previously written section 7 was an incomplete live-job checkpoint. The **final 1651050 detector-centroid result** from `1651050/window_centroids.csv` and `summary.json` is:

- 4/4 F444W early frames remain usable with originally adopted `catalog_position_2dg` centroids. For each, Gaussian fits succeed at half-window sizes **4, 5 and 7 pixels** but **fail at 9 pixels**, meaning a test of 4 through 9 pixels is **NOT globally stable** for this candidate.
- Among the four covering late F356W CAL frames, only **one** original **forced_gaussian** frame is usable at all (S/N **5.053**). Three other F356W frames fail the original accepted-centroid quality gate.
- On the sole usable late F356W CAL, the new, separately refitted centroid succeeds **only** at half-window sizes **4 and 5 pixels** and **fails at 7 and 9**. The two successful re-fits agree closely **with each other** (east span 0.00456 mas; north span 0.00533 mas), but are **offset from the original accepted forced Gaussian centroid** by about **(-1.95, -6.32) mas**, or **6.61 mas 2D**. Do not misinterpret the tiny between-window span as evidence of no cross-method centroid bias; the originally adopted *forced* estimator differs and two larger-window attempts do not converge.
- For that late F356W, the default center-of-mass versus accepted Gaussian detector offset is approximately **137.39 mas**, computed with the audit's fixed approximate 63 mas/pixel scale. Center-of-mass is not a validated alternative PSF estimator: this size is a *warning flag about sensitivity*, not a measured systematic correction.
- **387132 comparison**: 4 early F444W + 2 late F410M target centroids all accepted with `catalog_position_2dg`; all six succeed at **all four** window sizes (4,5,7,9). Largest within-exposure Gaussian span 0.01556 mas east / 0.01384 mas north. These are internal re-fit stability values, not absolute astrometric accuracies.

**Bottom line:** 387132's chosen Gaussian centroids show better test-window stability on every original accepted detector CAL. 1651050's sole faint cross-filter late-epoch centroid shows *limited convergence* and a cross-method difference from its accepted forced-fit coordinate. Coupled to existing full-epoch injection false positives (3/24 recovered stationary fakes at or above observed PM amplitude), the 1651050 displacement remains **not yet robustly attributable to intrinsic proper motion**; parallax versus motion cannot be separated using only two independent epochs.

**Do not** infer formal covariance, updated PM significance, a certified ePSF correction, or a revised paper/source classification from this test. No original paper numbers were changed.


## 10. Targeted 1651050 rejected-F356W follow-up — new successful REAL-CAL checks (8 Nov 2026)

**Full, reproducible details:** [LIU2026_1651050_F356W_TARGETED_AUDIT.md](LIU2026_1651050_F356W_TARGETED_AUDIT.md). New code `scripts/audit_1651050_f356w_joint.py` and `scripts/audit_1651050_f356w_register.py` executes on the same isolated plan branch.

- [Real-CAL shared-location model run 42486214684](https://github.com/physicsbahare/proper_motion_86/actions/runs/42486214684): **success** (artifact 13619821892).
- [Independent F356W-only field-registration run 42487059542](https://github.com/physicsbahare/proper_motion_86/actions/runs/42487059542): **success** (artifact 13620167168).
- P0 original rejection audit: **only one of four F356W exposures remains accepted**. A failed F356W has DO_NOT_USE on modeled target source core; two others have lower-than-required original forced-fit SNR (2.996 and 3.860); first rejected even weaker. Do not promote frames based on joint fitting.
- All four F356W CAL images exhibit positive, **jointly selected/conditional** source-model amplitudes near the target sky location. Using an analytic Gaussian and initial raw gWCS, the joint fit excluding the accepted image returns a position 12–22 mas from that original accepted centroid. This does **not** establish independent single-frame detections.
- New **same-F356W-filter local affine control registration** succeeds with **11, 8 and 11** field controls in the three rejected exposures (reference has 13); per-coordinate robust registration scatters are about **4–8 mas**. After this reference registration, four-image Gaussian-width-model offsets relative to the accepted F356W center are **(+2.94,+8.64)**, **(+0.30,+10.67)**, **(-2.40,+13.95)** mas at fixed Gaussian widths 0.8, 1.05 and 1.3 pixels. Excluding the accepted image, the three rejected frames give **(+3.87,+11.03)**, **(+0.62,+14.64)**, **(-3.14,+19.50)** mas, respectively.
- The new same-filter fit is promising as a **feasibility/sanity diagnostic**, not publication-ready astrometry: simplified Gaussian instead of independent color/SED-matched PSF, source location chosen near original accepted source, neighboring bright source/background contamination, DQ flags, conditional positive-amplitude bias, very small local-control N, no quantified absolute centroid covariance, and no cross-filter or parallax resolution.
- **Gate:** Hold any fresh injection campaign until a truly validated source-color/PSF model is available. Reusing the old compact-source empirical stamp would NOT eliminate the known color-dependent systematic. Keep original `main` and manuscript results unchanged.

**Recommended description of 1651050 today:** candidate with observed cross-filter two-epoch positional displacement and additional suggestive, registered same-filter F356W weak-frame signal; not a newly confirmed independent high-significance intrinsic-proper-motion brown dwarf. A third independent and suitable epoch would be most decisive.
