# Publication astrometry follow-up: 387132 and 1651050

**Aim:** Assess whether two-epoch apparent displacements can be labelled "proper motion" without a parallax term. This follow-up changes no existing result or classification.

## Why needed

The manuscript reports photometric distances of roughly 19.1 pc (387132)
and 11.2 pc (1651050). Those values would imply trigonometric parallaxes
of roughly 52 and 89 mas, respectively, if accurate. The previous
proper-motion estimates divide the measured two-epoch displacement by the
time separation; they do not fit parallax.

The 387132 baseline is ~1.992 years, close to an integer number of years,
so annual parallax is expected to cancel substantially. The 1651050 baseline
is ~1.362 years, so the two epochs are at different annual parallax phases.
Its apparent ~5.06 radial S/N displacement is *not* equivalent to 5.06
one-dimensional Gaussian sigmas.

## Already in repository

- Original exposure-level gWCS positions; DQ masking of compromised pixels
- Local affine registration with field controls
- 2D Gaussian / center-of-mass centroid comparison
- Exhaustive legal independent epoch pairs for the 42 new Step3 objects
- Neighbour-latching and centroid sensitivity audits
- 387132 associated with previously validated legacy ID 282040
- Fresh September 17, 2026 four-target astrometric verification available as
  GitHub Actions artifacts for run 35180229231, expiring October 17, 2026.

## New workflow

The on-branch workflow reruns only candidate 1651050, retains all exposure
diagnostics as an Actions artifact, and outputs a distance-conditional
Earth-barycentric parallax sensitivity table for both candidates.

Earth is used as a diagnostic approximation for the observer, **not** an
exact JWST L2 ephemeris. The table must not be mistaken for a measured
parallax or a joint PM+parallax fit. Obtaining another genuinely independent
epoch is the clean way to constrain both astrometric terms.

## Further high-priority checks

1. Reconfirm the 1651050 centroid in its late F356W exposure with
   empirical/effective PSF fitting and a nearby-neighbour null/control test.
2. Obtain a genuinely independent new epoch, preferably in F444W or a
   well-calibrated neighboring red filter; fit position, 2D PM and parallax.
3. Use a JWST-specific barycentric orbit for publication-grade parallactic
   corrections, propagate realistic covariance and systematic uncertainties.
4. Ensure the 125-object completeness statements track all insufficient and
   ambiguous cases; no usable PM fit does **not** mean stationary.

**Do not use a 2-object count as a proven complete UCD census.**
