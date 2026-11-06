# Candidate 1651050: publication decision and next astrometric gate (2026-10-08)

## Current evidence
- Public MAST third-epoch inventory: [run 37774486721](https://github.com/physicsbahare/proper_motion_86/actions/runs/37774486721). WCS coverage is **not** a detection.
- Repaired original MIRI CAL audit: [run 37792840090](https://github.com/physicsbahare/proper_motion_86/actions/runs/37792840090); [verified report](1651050_MIRI_THIRD_EPOCH_VERIFIED_20261008.md). 18/152 CAL images cover the target; 14 yield forced-aperture S/N, none an independent centroid. Highest per-exposure conservative S/N 2.16; later F1000W at MJD 61025.482705 has S/N 1.8703.
- Earlier empirical NIRCam ePSF fits show inconsistent F444W Gaussian-ePSF north offsets (-56.3 and +22.1 mas), F444W reduced chi-square 4.29 and 2.56, and a problematic F356W target centroid. Cross-filter systematic uncertainty is not established.
- Two NIRCam epochs cannot independently separate intrinsic proper motion from JWST-L2 parallax without a distance constraint; a zero-PM parallax-only interpretation remains possible conditional on distance.

## Required gates before any three-epoch fit
1. Search for **new public NIRCam Stage-2 CAL exposures** at the target coordinate, recording observation ID, exact EXPSTART, filter, detector, data rights and product URI. Verify target falls on each detector using original GWCS (not only observation footprint). Group exposures by independent visits, not exposure count. Prefer F444W/F356W and a distinct parallax phase.
2. On original SCI/ERR/DQ pixels, measure a **real source detection** (not a forced aperture at catalog WCS). Inspect neighbor contamination, PSF morphology, residuals, DQ flags and centroid stability versus fitting window/background. Report S/N, centroid covariance and goodness of fit; no centroid from a non-detection.
3. Independently register each visit to compact unsaturated field references with robust outlier rejection; report number of controls, leave-one-out residuals, transformation complexity, local distortion and filter/color systematics. Compare Gaussian vs empirical PSF centroids and propagate disagreements rather than selecting the favorable result.
4. Only with **three independent, reliable position epochs** and sufficient parallax-factor leverage, jointly fit tangent-plane position, mu_alpha*, mu_delta and parallax using **JWST-L2 ephemerides**, full astrometric covariance and reference-frame uncertainties. Check design-matrix rank/condition and compare with a zero-PM fit; do not equate numerical convergence with identifiability.
5. If no usable third centroid exists, report apparent two-epoch displacement, distance-dependent conditional PM, and the parallax-only alternative. Do **not** classify 1651050 as independently confirmed intrinsic PM.

## Publication-ready conservative conclusion
Candidate 1651050 shows a two-epoch NIRCam apparent displacement. Existing cross-filter centroid/registration systematics and uncertain distance prevent an independent attribution to intrinsic proper motion. The public MIRI observations include later WCS coverage but no reliable third-epoch astrometric centroid. A joint proper-motion/parallax solution is therefore **not justified by the currently validated data**. This does not negate the separately validated candidate 387132 result.

This document changes neither main nor the original classification tables.
