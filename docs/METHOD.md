# Method

## 1. Archive discovery

For every catalog RA/Dec we query public MAST imaging within 3 arcsec. The raw inventory is saved before any PM decision. HST coverage is recorded, but the primary automated PM fit does not silently mix HST and JWST because cross-instrument PSF/distortion systematics deserve a dedicated validation layer.

## 2. Independent epochs

Observation rows are grouped into visits when their archive times are separated by <=1 day. Dithers from one visit provide repeat measurements and uncertainty diagnostics; they are not treated as a multi-year PM baseline.

A PM fit requires a baseline of at least 30 days.

## 3. Pair selection

Priority is:

1. same-filter JWST/NIRCam, in this order: F444W, F410M, F356W, F277W, F200W, F150W, F115W, F090W;
2. if unavailable, the cross-filter JWST pair closest to F444W with an adequate time baseline.

A cross-filter pair is only useful if the target is independently detected/centroidable in both epochs. A blue-filter non-detection is not interpreted as motion.

## 4. Original detector pixels

For selected epochs, individual Stage-2 NIRCam `*_cal.fits` products are retrieved through MAST. The code opens them remotely with byte-range access and reads the embedded full JWST gWCS. Products that do not actually put the candidate on the detector are rejected.

For covering products the code reads only a small local SCI/ERR/DQ cutout. Whole CAL products are not stored in the repository.

## 5. Photometry and DQ

At the exact gWCS candidate coordinate:

- circular aperture radius: 3 pixels;
- local background annulus: 6-10 pixels;
- S/N is recorded from the ERR plane;
- a local empirical blank-aperture noise estimate is also measured.

For raw/artifact-sensitive photometry only DO_NOT_USE and SATURATED are hard-masked. JUMP_DET and OUTLIER remain diagnostic. For astrometry, DO_NOT_USE, SATURATED, JUMP_DET and OUTLIER are masked.

## 6. Target centroiding

Two target centroid chains are retained:

- 2-D Gaussian (`centroid_2dg`), primary;
- center of mass, independent robustness chain.

The 2-D Gaussian centroid uncertainty is estimated with Monte-Carlo perturbations from the ERR stamp.

## 7. Local controls and registration

Field sources are detected independently in each cutout. Controls must be away from the target and detector edge and have S/N >=10. Their positions are refined and converted through each exposure's gWCS.

Controls are matched to a reference exposure and a robust 2-D affine transform is fit with iterative clipping. At least six matched controls are required. Registration residual scatter is propagated into the target position uncertainty.

## 8. Cross-filter systematics

When the selected epochs use different filters, matched controls provide a rough local color proxy from their flux ratio. Residual astrometry is tested for correlation with this color proxy. A strong significant residual-color correlation is flagged as a cross-filter systematic.

The COM and 2DG PM solutions are also compared. Large centroid-method disagreement downgrades the result.

## 9. Proper motion

Registered target positions are combined separately within the early and late epochs. The displacement divided by the epoch baseline gives:

- `mu_alpha_cosdec_masyr`
- `mu_delta_masyr`

All early x late exposure pairings are saved as consistency diagnostics. Because the pairings share exposures, they are not counted as independent significance measurements.

## 10. Classification

Primary 2DG criteria:

- `MOVING`: 2-D significance >=5 sigma, pairwise direction consistency >=75% when available, and no strong cross-filter/centroid-method systematic flag.
- `CONSISTENT_WITH_ZERO`: significance <3 sigma and neither PM-component uncertainty exceeds 20 mas/yr.
- `AMBIGUOUS_SYSTEMATICS`: a fit exists but the cross-filter/color or centroid-method diagnostics are concerning.
- `AMBIGUOUS`: a fit exists but it is neither a secure mover nor precise enough to call consistent with zero.
- `INSUFFICIENT_DATA`: no independent epoch pair, no target detection in both epochs, inadequate local controls, or another failure that prevents a defensible PM fit.

These thresholds are deliberately conservative. The output tables retain the continuous measurements so the classification can be re-evaluated without repeating the archive work.
