# 1651050: final two-epoch sensitivity and publication decision (2026-10-10)

## Inputs and provenance

This is a **re-analysis of existing archived numbers**, not a new original-pixel measurement. Source: `docs/PUBLICATION_ASTROMETRY_FINAL_FINDINGS.md`, `docs/1651050_NEXT_ASTROMETRIC_DECISION_20261008.md`, and the previously successful JWST-L2 Horizons run [37771840466](https://github.com/physicsbahare/proper_motion_86/actions/runs/37771840466).

- Early F444W observations: MJD 60315.151156 and 60315.157742 (same visit/epoch).
- Later F356W observation: MJD 60812.794560. Two **independent epochs**, despite three accepted exposures.
- Archived annualized apparent displacement: east -55.2303, north +16.4082 mas/yr; formal component errors 11.5336, 9.9982 mas/yr. Formal errors **exclude** calibrated cross-filter PSF/registration covariance.
- JWST-L2 differential parallax factors: east -1.5987386, north +0.6994593. Baseline approximately 497.64 days (1.362 years).
- Photometric distance assumptions: 7.2, 11.2, 17.0 pc. These are **assumptions, not astrometric distances**.

## Distance-dependent motion (existing Horizons calculation)

| Assumed distance (pc) | Conditional east motion (mas/yr) | Conditional north motion (mas/yr) | Interpretation |
|---:|---:|---:|---|
| 7.2 | +107.74 | -54.89 | Requires substantial motion under fixed distance and archived centroids |
| 11.2 | +49.54 | -29.43 | Requires motion under fixed distance and archived centroids |
| 17.0 | +13.79 | -13.79 | Much smaller residual, sensitive to centroid systematics |

The existing forced-zero-PM parallax-only diagnostic gives 44.00 mas (22.73 pc) and residual chi-square 0.48 for 1 formal degree of freedom, but **22.73 pc lies outside the adopted photometric interval**. This is a degeneracy illustration, not an astrometric distance measurement and not evidence that zero motion is physically preferred.

## Centroid-systematics stress test

- The original F356W target has only one accepted exposure among four covering CAL images, with S/N 5.05, Gaussian-versus-center-of-mass separation 2.18 detector pixels (~138 mas), and a nearby JUMP_DET pixel.
- Independent local ePSF models for two F444W dithers shift the Gaussian north centroids by -56.31 and +22.09 mas, respectively: **78.40 mas disagreement between those ePSF-Gaussian offsets**, even within the early visit.
- F444W ePSF reduced chi-square values 4.29 and 2.56; F356W 1.06. These are not calibrated uncertainty estimates.
- The previously reported two-method annualized motion discrepancy is ~59.41 mas/yr, larger than the nominal conditional residual amplitude at 17 pc (~19.5 mas/yr). This comparison is a **robustness warning, not a Gaussian significance test**.
- Do not add 78.40 mas or 59.41 mas/yr as if they were statistically calibrated 1-sigma errors. They demonstrate the need for model/registration uncertainty calibration.

## Third-epoch gate

MAST inventories, original MIRI CAL images and HST FLT/FLC images were checked in prior runs. MIRI has 18/152 WCS-covering exposures, maximum conservative S/N 2.16, **no usable independent centroid**. HST has eight unique exposures; apparent F300W signal is not repeatable and has DQ/fit issues. **WCS coverage is not a detection.** There is currently no validated third independent epoch suitable for a joint two-component proper-motion plus JWST-parallax fit.

## What the current data support

1. **Report:** an apparent two-epoch NIRCam displacement with clearly labeled *formal* errors, the sensitivity of inferred intrinsic motion to assumed distance, and the method-dependent PSF/centroid systematics.
2. **Do not report:** a secure 5-sigma intrinsic proper motion, a measured trigonometric parallax or a uniquely determined Galactic tangential velocity for 1651050.
3. **Classification:** *astrometrically inconclusive*. Keep independent photometric/SED classification separate. Do not include 1651050 in a confirmed-proper-motion count.
4. **Next useful analysis, if revisiting:** same-filter/matched-color detector-level centroid and local reference-star registration with bootstrap/leave-one-out controls and an explicit cross-filter uncertainty budget. This may improve the *apparent displacement*, but **cannot alone resolve the two-epoch parallax/motion degeneracy**.
5. **For a true solution:** a new distinct-visit high-S/N NIRCam detection, preferably F444W or matched red filter at a different JWST-L2 parallax phase, with clean DQ, stable PSF and adequate astrometric field controls. Verify original GWCS, actual detection, and registration before fitting.

No main-branch classifications or tables are modified.

## Suggested manuscript wording

> Candidate 1651050 exhibits an apparent displacement between two NIRCam observing epochs. However, the inferred motion depends on the assumed photometric distance and on cross-filter centroid methodology. Local empirical-PSF and Gaussian fits exhibit substantial, exposure-dependent offsets, and archival HST and MIRI data do not provide a reliable independent third astrometric centroid. We therefore classify its astrometric status as inconclusive and do not interpret the two-epoch displacement as independently confirmed intrinsic proper motion or derive a joint astrometric proper-motion/parallax solution.
