# 1651050 real-CAL F444W physical-PSF plus compact control check

**Date:** 2026-11-08  
**Branch:** `plan/liu2026-minimal-astrometry-followups-20261108` (not `main`)  
**Source:** Public JWST MAST Stage-2 CAL detector frames, unchanged SCI/ERR/DQ  
**Workflow:** [GitHub Actions run 42491583765](https://github.com/physicsbahare/proper_motion_86/actions/runs/42491583765)  
**Code:** [audit_1651050_f444w_psf_controls.py](../scripts/audit_1651050_f444w_psf_controls.py)

## Question

Does the *early* four-frame F444W target position depend strongly on the assumed PSF compared to the accepted late F356W target centroid? Are analogous Gaussian-versus-STPSF shifts observed among nearby, independently selected high-S/N isolated **compact field objects**?

The F356W physical two-source sensitivity is documented in [LIU2026_1651050_STPSF_COLOR_NEIGHBOR_20261108.md](LIU2026_1651050_STPSF_COLOR_NEIGHBOR_20261108.md). That test alone did not address the reported large and inconsistent original empirical-ePSF versus Gaussian offsets for F444W in [PR #1 publication audit](https://github.com/physicsbahare/proper_motion_86/blob/publication-pm-parallax-audit-20261008/docs/PUBLICATION_ASTROMETRY_FINAL_FINDINGS.md).

## Defined method and boundaries

1. Reuse the *same originally target-accepted original CALs* only: four early F444W and one late F356W for the ID 1651050. No originally rejected source frame is promoted.
2. For each source and image, fit a target-centered **fixed-width Gaussian (σ=1.05 pixels)** and **detector-/filter-/position-specific STPSF** to exactly the same 13×13 original SCI/ERR/DQ patch, both with per-image signed linear source amplitude and a planar local background.
3. STPSF target source SEDs are hypothetical 500 K/900 K blackbodies; the 5700 K blackbody is an ordinary-star comparison for field compact objects. These are **not atmospheric models or SED measurements**, especially at molecular absorption bands.
4. Select up to three bright (field-segmentation S/N≥12), isolated (≥10 px from other detected controls, ≥16 px from target), defect-clean-core compact **control candidates** per original exposure. Exclude them from target-source fitting. **Being compact is not equivalent to being an independently verified star**; do not claim spectral or astrometric stationarity.
5. Record per-object Gaussian–STPSF centroid shifts in tangent-plane east/north mas from per-exposure gWCS, model fit χ² and amplitude, target position-boundary status, and per-exposure controls. Report if no credible controls exist; do not invent a PSF calibration from them.
6. Only as an **illustrative, unregistered raw-WCS method-sensitivity calculation**: annualize the difference between the F356W detector method offset and the median of the four F444W detector method offsets. This is **NOT a full two-epoch astrometric re-registration**, original F444W `catalog_position_2dg` is not necessarily the comparator's fixed-width Gaussian, source-color astrometric systematics are not calibrated, and no formal sigma or revised physical PM is inferred.

**Scientific acceptance gate:** A credible independent PSF calibration requires spectroscopically or morphologically well-vetted stars and a validated detector subpixel interpolation/phase residual test. Model convergence on a handful of compact objects does not satisfy this gate. If the model-to-model centroid shifts or their field-control scatter are large, do **not** run new injection or revise the manuscript PM. The two-epoch parallax degeneracy persists irrespective of PSF fit precision.

**Status:** new real-CAL workflow submitted; numerical outcomes added only after the workflow artifact is checked.


## Results checked from actual completed Actions artifact (8 November 2026)

- [GitHub run 42491583765](https://github.com/physicsbahare/proper_motion_86/actions/runs/42491583765): **success**. Raw artifact **13622405043**, `1651050-F444W-STPSF-independent-field-control-audit`.
- Four originally accepted F444W CALs have original target astrometric S/N **16.33, 21.05, 16.71, 19.49**; the one accepted late F356W has S/N **5.053**. Twelve isolated compact bright *field control candidates* were fitted on F444W (three per CAL) and three on F356W (one CAL). These objects are **not confirmed stars**, and counts contain repeat field sources across exposures.
- Every control passed the numerical STPSF fit, but **STPSF had higher weighted pixel χ² than the fixed-width Gaussian in all 15/15 controls**. For these controls, the **median 2D STPSF-minus-fixed-Gaussian displacement** is **1.818 mas (F444W, 12 fits)** and **1.779 mas (F356W, 3 fits)**. Coordinate median shifts (east,north) are (+0.128,+1.186) mas for F444W and (-0.456,-0.391) mas for F356W. No independent proof of accurate STPSF is obtained. These are fit-to-fit differences and the control sample may contain repeated objects, not 15 independent stars.
- **Target fitted fixed-width Gaussian versus STPSF**: under hypothetical 900 K spectra, the difference in an individual accepted F444W exposure is approximately (-1.01,+1.42), (-0.47,+1.69), (+0.96,-0.83), (+0.09,+1.76) mas (raw-WCS east,north). The accepted F356W exposure is (+1.31,+4.19) mas. The indicative late-minus-median-early method shift annualizes to (+1.10,+1.93) mas/yr over the frozen 497.6405-day baseline. This is **NOT** a re-registered or PSF-calibrated proper-motion correction.
- **Critical distinction: ORIGINAL paper-like `catalog_position_2dg` centroid versus the NEW fixed-width Gaussian/STPSF centroid.** The original F444W `...00004` has original pixel coordinate (199.6429,198.7402) versus STPSF 900 K (199.8143,199.3662), **0.649 pixel (~40.9 mas using nominal 63 mas/pixel)** difference; the new fixed-width Gaussian also differs by 0.631 pixel. This is **much larger than the STPSF-versus-fixed-width-Gaussian method difference of ~1–2 mas** in that same CAL. DO NOT ignore or conflate these estimators. Other accepted F444W raw detector original-versus-STPSF separations are **~1.1, 5.8, 1.0 mas** (files `...00003`, `...00001`, `...00002` respectively), nominal scale. F356W original forced Gaussian versus new STPSF differs **~4.8 mas**.
- The fitted 3×3 target detector cores include **5 astrometrically masked pixels in F444W `...00004`** and **4 in F444W `...00001`**, zero in F444W `...00003` and `...00002`; accepted late F356W has one. Some mask counts can differ with the exact fitted center selected. This is a warning about DQ-dependent centroid sampling, not proof the original measurement was incorrect.
- Only **four of five target STPSF trials** (900 K each) have lower weighted χ² than the fitted Gaussian; **0/15 independent compact-control trials** do. As controls are not spectroscopically verified and no pixel-phase validation exists, this test **fails the predeclared independent PSF calibration gate**.

### Critical next check, launched on the existing original estimator only

A tightly bounded **original pipeline leave-one-early-exposure-out test** is therefore necessary *before* interpreting apparent physical-PSF astrometry, particularly for F444W `...00004`. New script `scripts/audit_1651050_f444w_frame_loo.py` and [workflow run 42492250731](https://github.com/physicsbahare/proper_motion_86/actions/runs/42492250731) repeat full-CAL original `measure_exposure_recentered` + local control `register_exposures` + `infer_proper_motion` with all 4 early+1 late accepted, each early image omitted in turn, and both originally accepted early exposures with bad 3×3 target DQ pixels omitted together. The test explicitly preserves original registered controls and source quality; no STPSF centroid is used as publication substitute.

**Interim scientific status:** a stable *within-physical-PSF-model* F356W position does not settle the more consequential early F444W original-estimator discrepancy. Do NOT interpret `~1 mas/yr` raw fixed-Gaussian/STPSF difference as the whole astrometric systematic, promote 1651050 to a confirmed intrinsic PM object or modify frozen tables.


## Completed source-exposure leave-one-out check on original PM estimator

[Run 42492250731](https://github.com/physicsbahare/proper_motion_86/actions/runs/42492250731) **SUCCESS**, archived actual output **artifact 13622077969** `1651050-original-F444W-detector-LOO-robustness`. I inspected its `summary.json`, `original_estimator_f444w_leaveout_summary.csv`, and original-target DQ core quality. The analysis uses actual Stage-2 CAL detector images and the original code's `measure_exposure_recentered` + local-control `register_exposures` + `infer_proper_motion`; no new STPSF centroid has been substituted for the original PM.

| Scenario | Accepted F444W early frames | East PM (mas/yr) | North PM (mas/yr) | Nominal original-pipeline 2D score | Internal class |
| --- | ---: | ---: | ---: | ---: | --- |
| Original independent refit, all four early | 4 | **-61.74095** | **+17.15049** | 4.3377 | AMBIGUOUS |
| Drop `...00004` (5 hardbad core pixels) | 3 | **-55.40382** | **+17.79839** | 4.8554 | AMBIGUOUS |
| Drop `...00003` (0 hardbad core pixels) | 3 | **-62.98206** | **+17.93729** | 4.0461 | AMBIGUOUS |
| Drop `...00001` (4 hardbad core pixels) | 3 | **-61.91000** | **+16.03356** | 3.7903 | AMBIGUOUS |
| Drop `...00002` (0 hardbad core pixels) | 3 | **-64.37397** | **+18.07865** | 4.1053 | AMBIGUOUS |
| Drop both early F444W `...00004` and `...00001`, retain clean `...00003` + `...00002` | **2** | **-55.230296** | **+16.408212** | 5.0621 | MOVING |

All cases keep the **one originally accepted** late F356W CAL. The last *two-clean-early plus one-late* refit reproduces the frozen pair **(-55.2303,+16.4082) mas/yr** to about **0.00002 mas/yr** in each component; the original frozen baseline differs from refitted baseline by only ~0.00039 day. Such exceptionally close numerical agreement strongly suggests the earlier frozen pipeline derived its result using these cleaner F444W CAL exposures or an equivalent rejection. **It does not by itself certify historical exposure provenance**, which would require the exact archived original run's selected-exposure list and code/parameters.

**Important robust finding:** even with each individual originally accepted early F444W detector exposure excluded in turn, the source's measured **apparent two-epoch displacement retains the same east-negative/north-positive direction**. The east PM across LOO cases ranges -64.37 to -55.40 mas/yr; north +16.03 to +18.08 mas/yr. The result is not driven solely by the F444W `...00004` exposure whose unvalidated new STPSF center diverged by ~0.65 detector pixel from the original estimator.

**Still not confirmed intrinsic PM:** The model-based 5.062 original-pipeline score for the two-clean-frame case is **not** calibrated using enough independent ultra-red matched stationary sources or reliable across-filter PSFs. Field-control small N and the existing 3/24 recovered stationary injected false-motion diagnostic remain pertinent. JWST parallax can mimic some or all of the two-epoch apparent displacement under different assumed distances. The physically modeled F356W PSF and F444W compact-control fit exercises did not deliver independently certified ~mas absolute PSF errors.

**Decision:** Save these real-data sensitivity results and optionally present the **stable observed displacement and original-astrometry numerical reconciliation** to collaborators, while NOT updating `main`, the paper table, 5-sigma label or source classification without additional independent evidence.
