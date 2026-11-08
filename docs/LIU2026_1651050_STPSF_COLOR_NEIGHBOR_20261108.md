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
