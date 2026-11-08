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
