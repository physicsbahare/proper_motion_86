# 1651050 F356W physical-PSF and bright-neighbor sensitivity: bounded experiment

**Date:** 2026-11-08.  
**Branch:** `plan/liu2026-minimal-astrometry-followups-20261108`.  
**Only candidate:** 1651050.  
**Workflow:** [run 42490151272](https://github.com/physicsbahare/proper_motion_86/actions/runs/42490151272).  
**Script:** `scripts/audit_1651050_stpsf_neighbor.py`.  
**Source:** original JWST NIRCam detector `*_nrcblong_cal.fits` Stage-2 data served by public MAST.

## Scientific motivation and hard stop gates

Prior [same-filter registered four-frame Gaussian test](LIU2026_1651050_F356W_TARGETED_AUDIT.md) found evidence that three initially rejected F356W detector frames have a common faint-source location, but did not validate individual centroids or the extreme-red target's PSF. The accepted F356W frame has a bright neighboring field source around 11 detector pixels away in the original unmodified SCI stamp; field-source PSF wings and background structure can contaminate target fits. This experiment tests whether a more *physically motivated* PSF shape or a simple explicit neighbor model materially alters the common late-epoch centroid.

## Experimental design

- Compute **instrument + filter + detector-position-specific** NIRCam F356W simulated PSFs with **STPSF** and **synphot**, including physical detector distortion PSF extension `OVERDIST`; native detector `NRCBLONG` (alias `NRCB5`). Original detector positions and native pixel sampling are used. References: [STPSF documentation](https://stpsf.readthedocs.io/en/latest/usage.html) and [STPSF installation](https://stpsf.readthedocs.io/en/latest/installation.html).
- Use three **hypothetical** source SEDs, **500K, 900K and 5700K blackbodies**, synthesized through the real F356W instrument bandpass by STPSF/synphot. The **500/900K models are bracketing thought experiments, NOT measured Teff and NOT a validated atmospheric spectrum** for this object. The 5700K model is only a normal-color comparison; source colors may differ by molecular absorption even within F356W.
- Use original DQ masks, ERR weights and per-exposure local affine registration to the **sole originally accepted** F356W CAL (11/8/11 held-out field-control matches on three rejected images). The science target is excluded from registration.
- Fit one common sky centroid across four frames, **unconstrained signed target amplitudes** for each image (avoids silently forcing fake positives), optional explicitly modeled bright neighbor with free per-frame amplitude, and either constant or planar per-frame background.
- Repeat every PSF/scene combination with the **sole originally accepted F356W exposure excluded**. Inspect fitted centroid displacement, residual chi-square, per-frame amplitude signs, convergence, optimizer boundary flags, and sensitivity to source SED and scene assumptions. Signed source amplitudes mean one cannot treat this as a guaranteed target recovery.
- Compare results only as **source-model sensitivity**. Detector PSF pixel-phase interpolation is NOT independently calibrated to mas accuracy; compact reference control objects have not been confirmed spectroscopic stars. The PSF model does not encode archival observation-specific telescope wavefront evolution, exact atmospheric absorption, or a full deblended background scene.

## Prespecified acceptance criteria

The new physical-PSF results may motivate a more rigorous calibrated fit only if:
1. The source centroid is stable *across hypothesized SEDs* and changes little when the bright-neighbor term or a background gradient is included.
2. The three originally rejected F356W frames show positive (unconstrained) independent profile amplitudes *even after excluding the originally accepted frame*.
3. The residual model is plausible on SCI/ERR/DQ and local target/neighbor pixel stamps, with no extreme reliance on saturated or flagged pixels.
4. Before any paper PM update or new injection trial, PSF morphology/centroid bias must be independently checked against isolated **verified** point-like field controls and the early F444W frame, with target-color uncertainties. Numerical convergence or a successful workflow alone is not validation.

**No changes to `main`, published source counts, PM, parallax, or significance.**

**Initial execution status:** script + workflow committed and workflow running. Do not quote new physical-PSF numbers until inspecting the real-run result artifact.


## FINAL: physical STPSF real-CAL results (8 November 2026)

[Run 42490151272](https://github.com/physicsbahare/proper_motion_86/actions/runs/42490151272) completed **SUCCESS**, including successful synthetic solver test. I downloaded and independently inspected raw artifact **13621172168** (`1651050-STPSF-spectral-neighbor-sensitivity`): `physical_psf_model_comparisons.csv`, `source_amplitudes_by_frame.csv`, `neighbor_selection.json`, `registration_and_masks.json` and `model_provenance.json`.

**PSF provenance:** STPSF **2.2.0**, NIRCam **NRCB5/F356W**, detector positions (1912,483), (506,459), (381,424), (129,459); source *hypotheses* **500/900/5700 K blackbodies** through synphot bandpass and 7 simulated wavelengths/PSF. Exactly 12 independently simulated instrument-position/spectrum templates and **24** joint scene/leave-out fits. These hypotheses do **not** constitute a measured temperature/atmosphere spectrum.

**Real neighbor:** the most prominent local accepted-F356W source is **12.10 detector pixels** from the target reference position, smoothed-S/N proxy 11.09, at relative east/north **(+638.67,+422.83) mas**. This is an image-detected nearby object with no independently vetted classification.

### Centroid offsets relative to *originally accepted F356W forced Gaussian coordinate* (mas)

For **two STPSF-modeled point sources (target and neighbor) and per-image planar background**, all four images vs three *originally rejected images only*:

| Assumed blackbody | All four: east,north | Drop only accepted: east,north |
| --- | --- | --- |
| **500 K** | **(+1.83,+5.18)** | **(+2.59,+6.85)** |
| **900 K** | **(+1.97,+5.04)** | **(+2.78,+6.66)** |
| **5700 K** | **(+2.09,+4.93)** | **(+2.95,+6.50)** |

- Across blackbody scenarios 500–5700 K, the 4-frame centroid coordinate range is merely **0.27 mas east / 0.25 mas north**, conditional on this source and interpolation model. This is **NOT** a calibrated systematic error on true ultracool-dwarf SEDs. The 3-rejected-only coordinate range is 0.36/0.35 mas.
- At 900 K, including explicit neighboring-source flux rather than omitting it shifts the 4-frame centroid **(+0.51,+0.44) mas**, when both models use a plane background. This is a model comparison, NOT proof of zero true contamination.
- At 900 K, removing the sole originally accepted image changes the target+neighbor+plane centroid **(+1.97,+5.04) → (+2.78,+6.66) mas**, a vector change of **~1.8 mas**. Therefore three weaker exposures still support a **conditional** common location after matched-control F356W registration, but NOT 3 independent secure astrometric detections.
- All **24 model fits converged**, with no optimizer hitting its ±85 mas search boundary. The target-model amplitude was **positive in all included exposures** for all models; amplitudes are **signed/unconstrained** in each frame. This positivity remains **selection-conditioned** on the accepted-frame source seed.
- At 900 K, using all four images, the weighted pixel χ² declines **11948 → 7518** when adding the luminous neighbor to the planar-background target fit. The large χ² gain mainly models the obvious neighbor in the 31×31 stamp, **not** evidence of an intrinsic target PM. The χ² remains substantial; detector-correlated errors, unresolved sources and scene fidelity are not calibrated.

### Scientific evaluation

**Accomplished:** new physical instrument/filter/position-dependent and source-color-bracketed STPSF sensitivity comparison on all original late F356W CALs, with and without an explicitly modeled bright neighbor, with the accepted image separately excluded. It supports within-tested-model consistency of the **late-epoch F356W source location**.

**Still unproven:** real molecular-feature SED (blackbodies are crude proxies), field-star-calibrated STPSF PSF fidelity and detector subpixel interpolation, neighbor identity/shape and background residual structure, F444W early-epoch chromatic/ePSF centroids, absolute cross-filter systematic and parallax-vs-intrinsic-PM. Two independent epochs cannot uniquely separate parallax and true PM. Original 3/24 recovered stationary false-motion injection tail remains relevant.

**Stop gate:** Do **not** rerun arbitrary injections with unvalidated source PSFs or promote originally rejected frames into publication astrometry. Further substantive work requires isolated verified point-star PSF control and a consistent early F444W cross-filter source comparison. Until then, **1651050 remains an interesting apparent two-epoch displacement candidate, not a newly certified 5σ intrinsic mover**.

No change to the frozen published PM, number-density estimates, source counts, or `main`.
