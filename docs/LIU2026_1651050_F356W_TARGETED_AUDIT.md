# 1651050: original F356W CAL rejection diagnosis and bounded joint-centroid test

**Started:** 2026-11-08. **Branch:** `plan/liu2026-minimal-astrometry-followups-20261108`. **No paper astrometry has been modified.**

## Original exposure-level baseline, independently re-inspected from archived real-CAL CSV

Source: [existing cross-filter original-CAL run 41605446965](https://github.com/physicsbahare/proper_motion_86/actions/runs/41605446965), artifact `dave-cross-filter-astrometry-audit` id **13218202106**, raw `1651050/exposure_measurements.csv` (not invented or freshly measured). Filter F356W, epoch 2:

| CAL filename | Original aperture/astrometric S/N | Original forced Gaussian fit S/N | Forced width sigma x/y (pix) | Rejection / status | DQ at catalog center |
| --- | ---: | ---: | --- | --- | ---: |
| `jw05893016008_05101_00001_nrcblong_cal.fits` | 5.053 | **5.053** | 0.728 / 0.726 | **ACCEPTED forced Gaussian**, reduced χ²=1.334 | 0 |
| `jw05893016005_05101_00002_nrcblong_cal.fits` | 2.255 | 0.448 | 0.550 / 0.648 | **REJECTED:** not an accepted reliable target centroid; catalog-center DQ=1 (DO_NOT_USE) | **1** |
| `jw05893016005_05101_00001_nrcblong_cal.fits` | 5.086 | 2.996 | 0.976 / **0.550** | **REJECTED:** forced-fit S/N below required 4; fitted width at minimum boundary (0.55 pix) | 0 |
| `jw05893016005_03101_00001_nrcblong_cal.fits` | 3.220 | 3.860 | 0.645 / 0.815 | **REJECTED:** forced-fit S/N < required 4; although fit shape approximately compact, no accepted 2D centroid | 0 |

- Original fit thresholds are defined in [first26/forced_astrometry.py](../first26/forced_astrometry.py), notably S/N ≥4, 0.55<σ<1.8 pix, reduced χ²≤3, and astrometric positional errors ≤1.2 pixels.
- Do **not** infer that original aperture S/N 5.086 for the rejected second row is a source-centroid detection; its forced-fit S/N 2.996 and width-on-boundary violate the intended quality cuts.
- The one accepted F356W image has a JUMP_DET DQ-flag within four detector pixels, but not at the catalog center. It is still source-model sensitive and is only one independent detector image.
- The four early F444W exposures remain the accepted astrometric baseline; this new targeted step only revisits F356W feasibility.

## New real-CAL targeted audit now launched

- Script: [audit_1651050_f356w_joint.py](../scripts/audit_1651050_f356w_joint.py).
- Workflow: [audit_1651050_f356w_joint.yml](../.github/workflows/audit_1651050_f356w_joint.yml); [run 42486214684](https://github.com/physicsbahare/proper_motion_86/actions/runs/42486214684).
- First stage exports per-frame original acceptance, forced-fit diagnostic and a figure of original unmodified SCI, ERR, DQ stamps (with expected location derived from accepted F356W target coordinate).
- Strictly bounded **feasibility-only** source-centroid modeling: one shared sky position for four F356W exposure WCS locations, **independent nonnegative amplitudes and background** in each, *circular fixed-width Gaussian* with widths 0.8, 1.05, 1.3 pixels. It uses the original DQ masks and uncertainty map; no rejected frame is automatically promoted to a paper detection.
- Compare all four exposures with leave-one-frame-out solutions, **including removal of the sole accepted frame**. A solution dependent on that frame or highly width-sensitive **does not rescue the other three observations**.
- This is *not* a color-weighted STPSF/WebbPSF fit. A validated cross-filter physical PSF model would require a separate focused check after seeing the DQ/fit diagnostics. **No new color-aware injection is authorized yet**, because an unvalidated model would only repeat the earlier problem.

## Science decision rules

1. Do not change original exposure quality labels just because a shared-position fit converges. Demand association consistency, positive per-frame amplitudes, and stability *without* the originally accepted frame.
2. An original catalog-center DO_NOT_USE frame cannot become reliable merely through a joint fit; core pixel masking and WCS registration need careful assessment.
3. If the result is unstable across reasonable widths or exclusion of the accepted frame, label the attempt `NOT_RESCUED`, and stop before any larger injection experiment.
4. Only if stable, physically credible, and independent of the originally accepted frame: evaluate an appropriate target-color PSF and verify calibration controls before running even a small additional injection.
5. Two independent epochs cannot independently separate parallax and proper motion; no claim of measured parallax or five-sigma intrinsic PM.

**Status:** run submitted; no new real-CAL joint-fit measurements are asserted before the workflow completes and artifacts are inspected.


## New result 1: four-CAL *unregistered* shared-centroid feasibility (COMPLETED 8 November)

[GitHub Actions 42486214684](https://github.com/physicsbahare/proper_motion_86/actions/runs/42486214684) **success** with self-test success and raw artifact **13619821892**, `1651050-f356w-original-cal-joint-feasibility`. I downloaded and read the raw `summary.json`, `joint_loo_width_summary.csv` and `conditional_per_frame_amplitudes.csv`, and inspected the SCI/ERR/DQ panels. This is an **additional diagnostic**, not validation of previously rejected centroids.

- Four original F356W CAL frames all enter a **fixed common-position analytic Gaussian fit** with detector DQ/ERR masks. All four individually show *positive, conditional fitted amplitudes* at the **jointly selected location** for widths 0.8–1.3 pixels. Do **not** interpret the conditional amplitude SNR as an independent detection probability because the common location was chosen from the accepted frame and fit to the same data. All three originally rejected frames **remain rejected**.
- Relative to the *original accepted F356W forced-Gaussian sky coordinate*, the all-four fit moves by (east,north) **(+3.08,+8.88) mas at width 0.8**, **(+0.11,+12.53) mas at width 1.05**, and **(-2.58,+16.25) mas at width 1.3**. These shifts are model-dependent astrometric *sensitivity*, not error bars.
- Repeating the joint fit **without the sole accepted F356W frame**, the three originally rejected images return offsets **(+3.98,+11.36)**, **(+0.36,+16.78)**, and **(-3.30,+21.92) mas** for the same width sequence. The optimization converges without hitting the imposed ±160-mas coordinate bounds. This shows possible *shared-location information* in those weak frames, **not three separately validated centroids**, and does not calibrate field registration.
- There is a nearby much brighter source in the original frames, and the source stamps exhibit nonuniform noise/flags, so a constant-background circular Gaussian may be biased. The image with original catalog-center DO_NOT_USE remains unsuitable without independent DQ checks.
- Joint residual χ² across the four frames is **1242–1327**, with model width; these are *raw* weighted pixel χ² values, not a reliable goodness-of-fit probability given correlated CAL noise/neighbor residuals and simplified PSFs.

**Interim decision:** A catastrophic total nondetection in the three rejected frames is **not** established, but equally **no independent rescue** has been demonstrated. Before any claim, the same-filter CAL field registration must be tested (the above fit used initial per-exposure archive gWCS only).

## New result 2: same-filter field registration initiated

- Code: `scripts/audit_1651050_f356w_register.py`; workflow `.github/workflows/audit_1651050_f356w_register.yml`.
- [Run 42487059542](https://github.com/physicsbahare/proper_motion_86/actions/runs/42487059542).
- Each F356W original exposure is locally affine-registered to the accepted F356W exposure **using independently detected, high-S/N field sources, with the candidate excluded**. The analysis saves the number of matched field references, robust registration scatter, and fixed-Gaussian joint-centroid shifts **before and after registration** and **with/without the original accepted frame**, for widths 0.8, 1.05, 1.3 pixels.
- This directly tests whether the apparent agreement of rejected frames is merely uncorrected WCS alignment (a potential important systematic).
- **Gate:** no color-aware PSF or new injection until reference registration, source association, and PSF/background sensitivity have been checked. A good numerical fit alone cannot promote a DQ-rejected frame.
- No new observations are fabricated; no journal PM, error bars, parallax, or source classification has been changed.
