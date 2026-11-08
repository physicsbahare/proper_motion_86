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
