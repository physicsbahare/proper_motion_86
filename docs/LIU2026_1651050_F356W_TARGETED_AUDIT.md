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


## New result 2: same-filter registration test COMPLETED (8 November 2026)

[GitHub Actions run 42487059542](https://github.com/physicsbahare/proper_motion_86/actions/runs/42487059542): **success**, and raw artifact **13620167168** (`1651050-F356W-same-filter-registration-feasibility`). I independently inspected the actual archived `summary.json`, `registered_vs_raw_summary.csv`, `per_frame_registration_quality.json` and `registered_vs_raw_per_frame_amplitudes.csv`.

### Same-band reference field quality

Local affine registration targets the originally accepted F356W exposure. Excluding the science target, the **reference** has 13 locally detected field sources; the three other detector images have **11, 8, 11** matched and retained controls, respectively. Per-axis robust residual scatter of field-source matches:

- Rejected `...05101_00002`: east **4.03 mas**, north **7.96 mas** (11 matched/used);
- Rejected `...05101_00001`: east **5.95 mas**, north **5.12 mas** (8 matched/used);
- Rejected `...03101_00001`: east **4.72 mas**, north **6.56 mas** (11 matched/used).

These values are *empirical same-filter registration residuals* from only 8–11 field matches, **not an error bar on the source centroid**, a stationarity-certified galaxy null, or a cross-filter calibration. The source itself was excluded from frame registration.

### Location agreement and PSF-width dependence, relative to accepted F356W centroid

| Fixed Gaussian width (pixels) | Original four-frame WCS east/north (mas) | Registered four-frame east/north (mas) | Registered three rejected frames alone east/north (mas) |
| ---: | --- | --- | --- |
| 0.80 | (+3.08,+8.88) | **(+2.94,+8.64)** | **(+3.87,+11.03)** |
| 1.05 | (+0.11,+12.53) | **(+0.30,+10.67)** | **(+0.62,+14.64)** |
| 1.30 | (-2.58,+16.25) | **(-2.40,+13.95)** | **(-3.14,+19.50)** |

The fixed shared-position model converged in all these checks **without the source hitting the ±160 mas center-position bounds**. In this model, the rejection of the sole original accepted image does **not** erase the tentative cross-image common location. The within-F356W local affine registration changes the best-fitting location by a few mas rather than introducing a catastrophic disagreement.

However: this **does not constitute three validated independent source detections**. Every fit is seeded near the original accepted position and assumes a fixed circular Gaussian, nonnegative target amplitudes and local constant backgrounds. Faint close neighbors/background, CAL correlated noise, DQ flags, selection conditioned on a shared position, and uncertain PSF can bias location and detection evidence. The original weak-frame quality decisions remain in force.

### DQ-sensitive rejected frame

I also re-inspected the archived original SCI/ERR/DQ pixel arrays (artifact **11545878714**, run 37769062034) at the new full-four-frame width 1.05 modeled source coordinate:
- Rejected F356W `...05101_00002`: **DO_NOT_USE bit 1 at the modeled target center**, and **3/9** pixels in the central 3x3 contain astrometric hard-bad DQ flags. This exposure is particularly unsuitable for rescue, regardless of its positive conditional flux.
- Rejected `...05101_00001`: modeled center has DQ=0 and **0/9** bad 3x3 core pixels; nonetheless its independently attempted forced fit had S/N 2.996 and one Gaussian width hit the 0.55-pix lower bound.
- Rejected `...03101_00001`: modeled center has DQ=0 and **0/9** bad 3x3 core pixels; its independently attempted forced fit S/N was 3.860 (<4).
- The single already accepted `...16008_05101_00001` also had a central 3x3 bad flag within one pixel despite a clean modeled center, highlighting the need for rigorous PSF validation.

These are **rechecks of archived pixel stamps with source coordinates derived from the unregistered common Gaussian model**; they are not proof that a same-filter re-registration unflags pixels or that noise/scene PSF are controlled.

## Current scientific decision / stop gate

1. **P0 complete:** all four F356W original CAL images, the rejected-frame reason, and real weighted-pixel comparisons were checked. All original rejection labels remain unchanged.
2. **P1 partial but scientifically informative:** a fixed-position, same-filter control-registered joint fitting experiment is feasible with 8–11 local references per rejected frame. Its shifts are small, stable within the *specified simple Gaussian family*, including when the accepted frame is excluded. This motivates **a single bounded physical F356W PSF / neighbor-background validation**, but **does not itself certify the late centroid** or chromatic F444W/F356W PM.
3. **No P2 injection rerun yet:** at present there is no **independently validated color/SED-weighted PSF for this ultra-red target**, which the pre-specified P2 gate explicitly requires. Repeating artificial-source injections with the old unvalidated field-source stamps would not resolve the systematic and would risk overconfidence.
4. **Do not change paper values or call 1651050 a confirmed high-significance intrinsic mover.** Separating proper motion from JWST parallax still requires another independent epoch (or an externally validated distance).
5. Data needed for the next *physical-PSF* stage: an appropriate target spectral prior or uncertainty-bracketed SED, verified NIRCam F356W/F444W PSF and detector sampling at each source position, and independent control-source/neighbor-residual model. Original CAL images are public at MAST; **user-local FITS are not required for this completed stage**.

**Conclusion:** a promising common faint-source position across the late F356W CAL frames was found, including after local field registration and removal of the previously accepted exposure, but **three individual late astrometric centroids were not independently validated or recovered**. The evidence is stronger than "all three rejected frames contain no information", but not strong enough to change the published astrometry.


## New result 3: physically simulated STPSF + neighbor sensitivity (COMPLETED)

[STPSF run 42490151272](https://github.com/physicsbahare/proper_motion_86/actions/runs/42490151272) **success**; independent real-CAL artifact **13621172168**. All 24 fits over 3 hypothetical blackbody source spectra (500,900,5700 K), target-only versus explicit neighboring bright source, flat versus planar background, and four-frame versus original-accepted-frame-excluded cases converged. The source-model amplitudes were positive in each included frame, conditional on a source coordinate anchored to the only accepted F356W detection. The strongest neighboring peak lies 12.10 detector pixels away.

At 900 K with explicit neighbor and planar background, the **all-F356W** registered centroid is offset by **(+1.97,+5.04) mas** from the accepted original F356W coordinate; the **three rejected frames alone** give **(+2.78,+6.66) mas**, just ~1.8 mas farther. At fixed target+neighbor+planar background, changing the hypothetical source SED from 500 K to 5700 K moves the all-frame centroid by only 0.27 mas east and 0.25 mas north. This is a bounded **within-model sensitivity**, not a verified color calibration or astrometric systematic error. The bright neighbor substantially alters the pixel χ²; its identity and PSF residuals have not been externally validated.

**Full scientific report:** [LIU2026_1651050_STPSF_COLOR_NEIGHBOR_20261108.md](LIU2026_1651050_STPSF_COLOR_NEIGHBOR_20261108.md).

**Current status:** the *requested physical-PSF and neighbor-model sensitivity step is done*. The PSF accuracy, original cross-filter F444W astrometry and parallax split remain unvalidated. Do not call the three originally rejected exposures independently confirmed, do not modify `main` or manuscript PM and do not rerun injection with unvalidated PSFs.
