# Candidate 1651050: what would actually confirm intrinsic proper motion?

**Status (2026-10-08):** The two previously identified NIRCam visits permit measurement of an *apparent relative displacement* but do not independently separate linear proper motion and annual parallax. The three accepted CAL detector images are only two independent observing epochs (two F444W dithers during the same early visit and one F356W later).

## In-progress archive follow-up

- [Current public MAST original-WCS inspection](https://github.com/physicsbahare/proper_motion_86/actions/runs/37774486721) checks possible third-epoch NIRCam, MIRI and HST images, including actual CAL coverage and timestamps.
- [Original MIRI F770W/F1000W/F2100W visual/forced-aperture audit](https://github.com/physicsbahare/proper_motion_86/actions/runs/37774367059) checks whether the red object is actually visible in MIRI. A high forced-aperture flux is **not** a secure new astrometric position without identification of the same object, isolated PSF model, local registration, and systematics.
- Before the fresh search, the archived Sep 2026 original JWST NIRCam data comprised two independent observing periods (around MJD 60315 and MJD 60812); most bluer-band exposures lack centroidable target signal. In the archived inventory, MIRI and HST observations occur at other times; each needs a source-detection and precision-astrometry check.

## Scientific acceptance criteria for an archival third position

1. Independently timestamped actual exposure(s), not merely more dithers within a previously used visit.
2. Unambiguous positional association to the same red source, no nearby-neighbor latching or gross PSF mismatch.
3. Original detector SCI/ERR/DQ acceptable within the PSF core and enough high-S/N detection for precise centroiding.
4. Empirical reference control sources in the original field establish an astrometric transform relative to the previously used frame. Test null-source shifts, leave-one-control-out, and orientation/filter-dependent residuals.
5. Quantify method-dependent offsets (Gauss vs empirical PSF, PSF/reference-color dependence, structured background) and propagate uncertainties.
6. Fit joint 2D position, linear 2D motion, and JWST observer parallax only after ≥3 genuinely independent and non-degenerate epochs. With exactly 3 epochs there are only six coordinate measurements for five astrometric parameters; a fourth epoch at a different parallax phase would be much stronger. MIRI/NIRCam cross-instrument offsets must be accounted for explicitly.
7. Compare true-motion and parallax-only models; report model dependence, not just a 2D displacement signal mislabeled '5-sigma proper motion.'

## If no suitable archival centroid exists

Obtain a new JWST/NIRCam epoch, preferably F444W with multiple dithers and target S/N well above ten, or another comparably red filter after validating wavelength-dependent centroid bias. A follow-up at a different annual parallax phase from the two existing observing dates is preferable; two newly separated phases are better for robustly fitting parallax and PM. Use multiple local field stars, explicit bad-pixel screening, and consistent PSF-fitting in *all* epochs rather than updating only the faint F356W source.

Until an independent epoch exists, manuscript language should remain **'candidate exhibiting a significant apparent displacement between F444W and F356W images'** rather than a confirmed intrinsic proper-motion brown dwarf. This conclusion does not invalidate the number-count prediction of roughly two selected Galactic UCDs; it means the predicted *population count* is not object-specific astrometric confirmation.
