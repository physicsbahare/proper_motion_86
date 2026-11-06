# Original-registration full-epoch held-out galaxy tests (real JWST, 2026-11-06)

## Independent provenance

Live GitHub Actions workflow [42183463607](https://github.com/physicsbahare/proper_motion_86/actions/runs/42183463607) **SUCCESS for both candidates**, using real MAST NIRCam Stage-2 CAL exposures.

- [387132 artifact](https://github.com/physicsbahare/proper_motion_86/actions/runs/42183463607/artifacts/13486716696)
- [1651050 artifact](https://github.com/physicsbahare/proper_motion_86/actions/runs/42183463607/artifacts/13486608495)
- Regression [3 synthetic tests passed](https://github.com/physicsbahare/proper_motion_86/actions/runs/42183348454)
- Original code: `scripts/run_galaxy_null_original_frame.py`, `src/pm86/galaxy_full_epoch.py`.
- Source identities: 17 user-COSMOS2025-matched, LePHARE TYPE=0 *photometric galaxy candidates*, in `data/field_controls_priyanka/cosmos2025_photometric_galaxy_controls_17.csv`.

**Method:** retain the original mixed high-S/N field controls for robust astrometric affine registration, but for every candidate reference galaxy, exclude its matched detector-source record from **every** exposure's registration before measuring the held-out galaxy's apparent PM. Register and combine all original target-centroidable exposures using the repository's original target estimator. This is **NOT galaxy-only registration**; some remaining registration sources may be stars. The alternative all-galaxy-reference fit is documented separately in `PRIYANKA_FULL_EPOCH_GALAXY_REAL_RESULTS_20261106.md`.

## Final verified measured diagnostics

| Quantity | 387132 (F444W/F410M) | 1651050 (F444W/F356W) |
| --- | ---: | ---: |
| Original target-usable early / late CALs | **4 / 2** | **4 / 1** |
| Frozen manuscript PM east/north, mas/yr | -22.900 / +2.400 | -55.230 / +16.408 |
| Original-reference independent full-epoch target refit east/north, mas/yr | **-19.983 / +2.056** | **-61.741 / +17.150** |
| Refitted formal errors east/north, mas/yr | 4.101 / 6.320 | 15.397 / 10.369 |
| Automatic classification in independent refit | **AMBIGUOUS** | **AMBIGUOUS** |
| Nominal 2D quadrature score, **not a calibrated Gaussian sigma** | 4.883 | 4.338 |
| Individual source-held-out galaxy tests spanning both epochs | 9 | 8 |
| Galaxy tests with **all original target-usable CALs** | **6** | **7** |
| Median null displacement east/north, mas | **-0.921 / +2.219** | **-1.760 / +0.156** |
| Full-epoch control null scatter east/north, mas | **10.983 / 8.796** | **18.817 / 15.347** |
| Full-epoch galaxy displacement covariance (mas²) | [[120.632,-33.567],[-33.567,77.367]] | [[354.092,203.385],[203.385,235.542]] |
| Smallest empirical source-level tail resolution | 1/7=0.1429 | 1/8=0.125 |

**The two target displacements persist** on real images, with substantially similar directions under both reference-population selections. They are not established to new 5-sigma significance by these tests. The mixed reference frame does not independently isolate filter-dependent PSF shifts.

### How sensitive are these estimates to individual controls?

Using the full-exposure held-out galaxy samples, the raw coordinate scatters' **leave-one-galaxy-out** ranges are:

- **387132 (N=6):** east **7.97–12.18 mas**; north **6.74–9.82 mas**, maximum fractional change **27.5%**.
- **1651050 (N=7):** east **10.20–20.59 mas**; north **6.24–16.78 mas**, maximum fractional change **59.4%**.

**1651050 source 726373** has a large apparent residual of `(-43.76, -35.71) mas` in this mixed-registration test. This one source strongly affects the empirical covariance; its morphology, isolation, color, detector association and centroiding need visual inspection. It must not be dismissed solely because it increases scatter.

**1651050 source 726167** had a larger catalog-position match offset (~0.124 arcsec), and should also be visually vetted before adopting its galaxy label as a secure cross-band astrometric control.

### Candidate leave-one-galaxy-out fit robustness

Under original mixed registration with each of the 9 / 8 vetted galaxy sources excluded in turn:

- **387132:** all 9 leave-one-galaxy-out target fits retained 4+2 exposures. East PM range **-21.31 to -18.90 mas/yr**, north remains positive (**+2.15 to +3.63 mas/yr**). The nominal classification is generally ambiguous rather than independently secured.
- **1651050:** all 8 target fits retained 4+1 exposures. East PM range **-73.19 to -56.70 mas/yr**, north positive (**+9.00 to +17.99 mas/yr**). One nominal result crosses the 2D algorithm's `MOVING` threshold; the other seven remain `AMBIGUOUS`. No source-level probability calibration supports treating the 5-sigma threshold as a Gaussian z-value.

### What is and is not resolved

**Supported**: these source-selection changes do not erase either target's apparent displacement; the nearby galaxy-labelled reference objects provide empirical held-out frame residuals. 387132's direction is relatively insensitive to individual galaxy omissions.

**Still not resolved**: uncalibrated astrometric significance, very small source N, cross-filter PSF color extrapolation to ultra-red targets, source-specific centroid/re-association biases, independent spectroscopic stationarity confirmation, parallax, and target/epoch-correlated noise. For 1651050, **only one usable late F356W detector image remains**, and earlier real CAL artificial-source injection had a sizable static false-motion recovery fraction. The held-out galaxy null cannot rule out that source-specific failure mode.

**Science decision**: Do not modify the frozen published PM table or infer additional Galactic count fractions from this work alone. The result favors a persistent measured *positional displacement* for both targets. It does not certify both as astrometrically confirmed high-significance stellar motions, especially 1651050. No systematic covariance should simply be added to reported registration uncertainties without separating double-counted contributions.

## Suggested concise collaborator reply

> I rechecked both candidates with the original multi-exposure JWST astrometric pipeline, holding out COSMOS2025-selected galaxies from the registrations used to measure each galaxy's displacement. The target motions remain similar to our original values. The six and seven fully covered galaxy controls give coordinate scatters of roughly 9–11 mas around 387132 and 15–19 mas around 1651050. The result supports the presence of measured positional shifts, but the small control samples, cross-filter PSF effects and the single usable late F356W exposure for 1651050 mean that a revised detection significance is not yet justified.
