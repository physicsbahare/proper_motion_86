# 1651050 MIRI original-CAL third-epoch audit: validated rerun (2026-10-08)

Source: [GitHub Actions run 37792840090](https://github.com/physicsbahare/proper_motion_86/actions/runs/37792840090), artifact `1651050-MIRI-original-CAL-third-epoch-exploration` (ID 11558236857). The repaired script separates `coverage_status` from the forced-aperture `status`.

## Verified results from the exported CSV (152 products)

- **18/152** original Stage-2 MIRI CAL products have WCS coverage at the catalog coordinate. The other **134** are outside the detector. WCS coverage is not source detection.
- **14/18** covering products have a computed forced-aperture statistic; **4/18** have insufficient locally valid pixels. No independent MIRI centroid was measured.
- F770W (MJD ~60314.975-60315.013): **8 covering**, 4 measured; conservative aperture S/N = 0.41, 1.03, 1.16, 0.93.
- F2100W (MJD ~60813.917-60813.941): **6 covering**, 6 measured; conservative aperture S/N = 0.39, -0.69, 0.75, 1.22, -0.16, 0.64.
- F1000W (MJD ~60813.948-60813.960): **3 covering**, 3 measured; conservative aperture S/N = 2.16, -0.07, 0.78.
- **Later F1000W** (MJD **61025.482705**): **1 covering** original CAL image; conservative aperture **S/N = 1.8703**.

The maximum single-exposure conservative aperture S/N is **2.16**. These are exploratory forced-aperture values at the catalog WCS coordinate, not PSF-fitting detections, aperture-corrected photometry, or astrometric measurements. Correlated noise, structured backgrounds, possible catalog-to-epoch displacement, and filter-dependent PSFs remain relevant. The results cannot be interpreted as formal nondetection limits or combined-epoch S/N without a validated combination and noise model.

## Astrometric interpretation

**No credible third astrometric position is obtained from these MIRI CAL products.** The later F1000W exposure is a *third observing epoch*, but its forced S/N of 1.87 is too low to justify a reliable centroid and local reference-frame registration. Consequently it cannot be added to a joint 5-parameter position/proper-motion/parallax solution.

The two-epoch NIRCam displacement for 1651050 remains conditional on assumed distance and cross-filter PSF/registration systematics. It is **not independently confirmed intrinsic proper motion**. Do not interpret 18 covering exposures as 18 detections or use the low-S/N MIRI values as astrometric constraints.

## Next decisive observation/analysis

Seek an independently measured higher-S/N NIRCam centroid at a distinct parallax phase (ideally in the same filter as an earlier epoch), with exposure-level PSF/ePSF residual checks, color-dependent centroid offsets, DQ inspection, and independent local-field registration. Only then jointly fit proper motion and JWST-L2 parallax and compare a zero-PM model.

This note is documentation of existing measured output only; it does not change candidate classification or main-branch science tables.
