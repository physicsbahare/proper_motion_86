# Real JWST full-epoch galaxy-control astrometry (2026-11-06)

## Provenance and status

**Real MAST Stage-2 NIRCam CAL exposures**, processed by the original `measure_exposure_recentered` target centroid chain, with the *original frozen target-centroidable epoch-pair selection*. Analysis: `src/pm86/galaxy_full_epoch.py` + `scripts/run_full_epoch_galaxy_null.py` using 17 nearby locally COSMOS2025-matched `LEPHARE TYPE=0` sources in `data/field_controls_priyanka/cosmos2025_photometric_galaxy_controls_17.csv`. Run: [42182049157](https://github.com/physicsbahare/proper_motion_86/actions/runs/42182049157), with [387132 artifact](https://github.com/physicsbahare/proper_motion_86/actions/runs/42182049157/artifacts/13485373580) and [1651050 artifact](https://github.com/physicsbahare/proper_motion_86/actions/runs/42182049157/artifacts/13485713433).

Code regression tests: [4 passed in run 42181946166](https://github.com/physicsbahare/proper_motion_86/actions/runs/42181946166). 1651050 live job concluded **success** (diagnostic available); 387132 live job concluded **failure by design** because no *fully complete* six-exposure held-out galaxy null could be obtained. **This failure is a scientific coverage limitation, not a failed catalog match or source-zero-motion result.** Both jobs uploaded their complete diagnostics. Neither candidate's frozen paper values nor its photometric number-density model were modified.

## Target measurements in same original detector epoch sets

| Property | 387132 | 1651050 |
| --- | ---: | ---: |
| Filters | F444W / F410M | F444W / F356W |
| Target-usable CALs early / late | 4 / 2 | 4 / 1 |
| Original frozen east / north proper motion (mas/yr) | -22.900 / +2.400 | -55.230 / +16.408 |
| Full-epoch galaxy-only reference fit east / north (mas/yr) | **-21.103 / +1.542** | **-62.486 / +17.530** |
| Galaxy-reference-fit formal PM errors east / north (mas/yr) | 4.128 / 4.440 | 14.344 / 9.225 |
| Full-epoch code classification | MOVING | AMBIGUOUS |
| Uncalibrated 2D quadrature score (NOT Gaussian sigma) | 5.124 | 4.753 |
| Source-level held-out photometric-galaxy candidates tested | 9 | 8 |
| Held-out galaxy nulls with measurements across both epochs | 9 | 8 |
| Held-out galaxies with **all** target-usable exposures included | **0** | **7** |

The galaxy-reference fitted displacements are broadly similar to the frozen target displacement directions, especially in east. They are not an exact like-for-like reproduction of the paper's original estimator, nor a final measurement.

## 1651050 full-epoch galaxy null (seven complete controls)

Among the **seven** held-out photometric galaxies with all four early and the one original target-usable late CAL, each excluded from its own affine registration:

- Median apparent displacement (east/north): **-5.361 / +0.119 mas**.
- Coordinate scatter (east/north): **11.620 / 13.862 mas**.
- Empirical residual covariance: `[[135.021,75.392],[75.392,192.152]] mas²`.
- Smallest possible source-level empirical tail resolution even for seven independent galaxies: `1/(7+1)=0.125`, **not** a 5-sigma calibration.
- The eighth galaxy (726328) was missing one early exposure and therefore excluded from the strict seven-galaxy covariance.
- Omitting any one of the eight galaxy controls leaves a negative east motion: `-69.20` to `-56.15 mas/yr`; north remains positive, `+10.09` to `+18.83 mas/yr`. The code still yields `AMBIGUOUS` in seven of the eight held-out target fits, and `MOVING` in one, based on nominal internal criteria only.
- **Critical**: there is only one usable late F356W exposure. The galaxy null does not validate that one faint ultra-red source's F356W centroid identification, source-specific chromatic PSF error, or a true 5-sigma tail. Earlier real-CAL injection-recovery controls still warrant caution.

## 387132 full-epoch target fit and strict null limitation

Using all nine candidate reference galaxies, the target fits **all six usable CAL exposures** (four early plus two late), giving `(-21.103,+1.542) mas/yr`, close to frozen `(-22.9,+2.4) mas/yr`, with no code-detected 2DG/COM disagreement.

Six of the nine galaxies were detected in all six CALs, but **excluding each held-out galaxy from every registration** leaves fewer than the six independent fitting controls required for the second late exposure in the held-out tests. The other three reference galaxies were themselves detected in only five CALs. As a result **none** of the nine held-out galaxy nulls has the same six exposures as the target. This is not equivalent to zero source motion, and no six-CAL full-epoch null covariance can be quoted.

The nine partial held-out galaxy nulls are *four early + one late* measurements, with **diagnostic-only** median displacement `(+0.615,+0.904) mas` and coordinate sample scatter `(10.334,14.333) mas`. They are not comparable to the six-CAL full-epoch target fit as a calibrated stationary distribution. Under nine leave-one-galaxy-out fits, the target's east PM remains negative (`-21.08` to `-17.90 mas/yr`), but the fits themselves frequently drop the second late CAL and often flag 2DG/COM disagreement. **Do not quote those mixed 5-/6-CAL jackknifes as a formal detection significance.**

## Scientific decision

1. Both measured **displacements persist when registration is restricted to photometric galaxy candidates**. This is encouraging about the reference-source population and means an apparent shift is not obviously caused solely by some moving reference stars.
2. It does **not** establish a final cross-filter ePSF systematic floor or a recalibrated 5-sigma detection for either object. Small galaxy N, non-independent calibration, galaxy-versus-ultra-red-source color/PSF mismatch, and per-exposure noise remain.
3. For 1651050, one F356W exposure still dominates the late epoch and prior injection/recovery sensitivity tests must be reconciled before asserting astrometric confirmation.
4. For 387132, enlarge/validate the local galaxy reference catalog (or test independent alternative registration with enough controls) so a held-out source remains measurable across **both** late CAL exposures without relaxing the validated minimum six-control affine fit.
5. No proper-motion-calibrated Galactic fraction or change to photometric number density is justified from this result.

### Reply to Priyanka, accurate as of this audit

> I repeated the astrometric analysis with COSMOS2025-selected galaxy controls using the original JWST epochs. The candidate displacements remain similar: approximately (-21.1,+1.5) mas/yr for 387132 and (-62.5,+17.5) mas/yr for 1651050. For 1651050, seven held-out galaxies measured over the available four-plus-one exposures show about 12–14 mas coordinate scatter. For 387132, the full six-exposure target fit is stable, but the held-out galaxy sample is too small to calibrate a matching six-exposure null. These tests are encouraging, although we still need to account for source-specific cross-filter centroiding effects, especially the single usable F356W exposure for 1651050.
