# Publication audit of the two putative Galactic UCD movers

**Date:** 2026-10-08. **Branch:** `publication-pm-parallax-audit-20261008`. **Status:** REVIEW REQUIRED, scientific classifications on `main` are untouched.

## 1. JWST observer parallax, not Earth's orbit

The new Horizons scripts query the actual JWST ephemeris (**target `-170`**, reference center `@0`, equatorial ICRF vectors, observation times converted UTC MJD to TDB JD). The resulting annual-parallax component is subtracted *conditional on adopted photometric distances*, with a full provenance record and small input manifests.

### 1651050 (1.362-year F444W/F356W baseline)

Original apparent annualized displacement: `(-55.2303,+16.4082)` mas/yr; formal 1D component errors `(11.5336,9.9982)` mas/yr.

Differential JWST parallax factors `(-1.5987386,+0.6994593)`.

| Photometric distance assumed | Conditional intrinsic mu_east | Conditional intrinsic mu_north | 2D radial residual statistic |
|---|---:|---:|---:|
| 7.2 pc | +107.74 mas/yr | -54.89 mas/yr | 10.84 |
| 11.2 pc | +49.54 mas/yr | -29.43 mas/yr | 5.21 |
| 17.0 pc | +13.79 mas/yr | -13.79 mas/yr | 1.83 |

The conditional *zero-proper-motion* parallax-only least-squares model gives `pi=44.00 mas` (distance `22.73 pc`) with residual chi-square `0.48` for 1 remaining degree of freedom (`p=0.49`). The 22.73 pc value lies outside the quoted 7.2–17 pc photometric interval, so the parallax-only fit **is not a measured distance** and does not supersede the atmosphere-fit distance. This demonstrates two-epoch degeneracy: **two-epoch astrometry alone cannot establish nonzero intrinsic proper motion**.

All p-values and radial statistics are conditional on fixed photometric distance, uncorrelated coordinate errors, and the archived formal astrometric error model; they do not account for cross-filter systematic covariance or uncertainty in the assumed distance. Do not write "5-sigma confirmed intrinsic PM".

[Successful JWST L2 calculation](https://github.com/physicsbahare/proper_motion_86/actions/runs/37771840466).

### 387132 / legacy 282040 (1.992-year F444W/F410M baseline)

Six original reference CAL centroids were used to sample JWST observer vectors. Under the nominal photometric distance 19.1 pc, differential parallax contributes only `(+1.513,-0.831) mas`. The conditional intrinsic motion remains `(-23.66,+2.82) mas/yr`, consistent with the original `(-22.9,+2.4) mas/yr`. The previously validated mover classification survives this parallax check. This is a conditional parallax sensitivity calculation, not a direct trigonometric parallax determination.

[Successful six-exposure JWST L2 calculation](https://github.com/physicsbahare/proper_motion_86/actions/runs/37772665878).

## 2. Original detector pixels and empirical effective PSF

The September 2026 original DQ-aware audit finds **one accepted F356W target centroid among four covering CAL exposures**, at `jw05893016008_05101_00001_nrcblong_cal.fits`: S/N `5.05`, original accepted 2D Gaussian/center-of-mass separation `2.18 pixels` (`~138 mas`). The accepted centroid's *pixel itself* is DQ-clean, but `JUMP_DET` is present within four pixels. Faint neighbors and background structure are visible.

[Original six CAL SCI/ERR/DQ panels and reviewer checklist](https://github.com/physicsbahare/proper_motion_86/actions/runs/37769062034).

The local empirical ePSF fits were run **independently for each of the three accepted CAL images** using photutils EPSFBuilder, selecting 5–6 isolated compact high-S/N local comparison sources. The builder converged at 7, 10, and 12 iterations (tolerance 0.01 detector pixels). **These reference sources are selected for compactness but not spectroscopically known stars.**

| CAL exposure | Filter | ePSF ref sources | ePSF minus Gaussian, east | ePSF minus Gaussian, north | ePSF reduced chi-square |
|---|---|---:|---:|---:|---:|
| jw01727065001_04101_00002_nrcblong | F444W | 5 | +3.62 mas | -56.31 mas | 4.29 |
| jw01727065001_04101_00003_nrcblong | F444W | 6 | +5.24 mas | +22.09 mas | 2.56 |
| jw05893016008_05101_00001_nrcblong | F356W | 6 | -10.67 mas | -16.80 mas | 1.06 |

The F356W empirical-PSF and Gaussian fits are closer than the F356W center-of-mass position, but F444W shifts are large and inconsistent between dithered exposures. The F444W fits also have elevated reduced chi-square, indicating model inadequacy relative to pixel errors. Consequently **these ePSF fits do not produce a robust, independent cross-filter PM measurement**. The model/centroid-method spread must not be ignored or treated as calibrated Gaussian error. Also, the early F356W 4-iteration model differed from the 12-iteration model, indicating model dependence even within the nominally accepted F356W frame.

[Successful higher-iteration empirical-ePSF analysis](https://github.com/physicsbahare/proper_motion_86/actions/runs/37772444851).

## 3. Publication decision

- **387132:** Astrometric mover supported by existing checks; JWST parallax contribution small at its nominal photometric distance. Keep as the stronger astrometric detection.
- **1651050:** Two-epoch **apparent displacement**, but *intrinsic proper motion is not independently confirmed*. Photometric-distance assumptions, parallactic displacement, weak F356W coverage, model-dependent ePSF fits and cross-filter registration remain limiting. Describe as an astrometrically interesting **candidate pending confirmation**, not a secure 5-sigma PM brown dwarf.
- A third genuinely independent epoch, preferably a comparable red NIRCam filter and enough field controls, permits a joint position/PM/parallax fit and would be the most decisive next step. A careful WebbPSF/STPSF color-dependent model and explicit ePSF reference-star validation would further test the single F356W centroid and cross-filter biases.
- The predicted ~2 Galactic UCDs selected by Dave's photometric cuts is a **population expectation**, not an independent validation of the classification of either particular object.

No original astrometric measurements, classifications, or main-branch result tables were changed by this audit. Review all downloaded original-pixel panels and local ePSF model files before merging this draft PR.
