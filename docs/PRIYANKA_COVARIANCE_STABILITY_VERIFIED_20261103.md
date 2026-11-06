# Priyanka audit: latest covariance-stability verification (2026-11-03)

## Provenance

Successful real-data GitHub Actions run [41633534458](https://github.com/physicsbahare/proper_motion_86/actions/runs/41633534458), branch commit `5637f5c9c29a2a94f1bf6f4490b555c42ac2fd83`, artifact `dave-cross-filter-astrometry-audit` (ID 13231110402). Both per-candidate `summary.json` files and the independent refit reconciliations were examined. This report is **diagnostic only**.

| Diagnostic | 387132 | 1651050 |
|---|---:|---:|
| Filter pair | F444W/F410M | F444W/F356W |
| Epoch baseline (days) | 727.4425 | 497.6405 |
| Retained reference controls | 10 | 9 |
| East residual median (mas) | -0.854 | +6.980 |
| North residual median (mas) | +4.024 | +0.154 |
| East residual scatter (mas) | 9.979 | 15.580 |
| North residual scatter (mas) | 10.849 | 12.233 |
| Full empirical residual covariance (mas²) | [[99.581,-25.390],[-25.390,117.701]] | [[242.751,138.655],[138.655,149.635]] |
| Leave-one-source-out east scatter range (mas) | 6.936–10.581 | 13.965–16.524 |
| Leave-one-source-out north scatter range (mas) | 8.259–11.506 | 8.690–13.061 |
| Maximum fractional scatter change | 30.50% | 28.96% |
| Reference position minor/major singular-value ratio | 0.633 | 0.781 |
| Reference positions nearly collinear? | No | No |
| Target instrumental color proxy within retained controls? | Yes | **No** |
| Registered early / late exposures | 4 / 2 | 4 / **1** |
| Frozen mu east ± sigma (mas/yr) | -22.900 ± 4.100 | -55.230 ± 11.534 |
| Frozen mu north ± sigma (mas/yr) | +2.400 ± 5.800 | +16.408 ± 9.998 |
| Independent refit mu east ± sigma (mas/yr) | -19.983 ± 4.101 | -61.741 ± 15.397 |
| Independent refit mu north ± sigma (mas/yr) | +2.056 ± 6.320 | +17.150 ± 10.369 |
| Frozen 2D quadrature score (not Gaussian sigma) | 5.60 | 5.06 |
| Quadrature sensitivity scenario (not Gaussian sigma) | 3.49 | 4.56 |

**Interpretation:** The reference geometry is two-dimensional, so the covariance instability is not simply caused by collinear reference positions. Both empirical covariance estimates change by >20% when one retained source is omitted. With only 9–10 controls, this scatter is not a calibrated systematic floor. The F356W candidate's target instrumental color proxy (1.781) lies outside its controls' range (-1.211 to +0.200), so even a null color trend among controls cannot validate extrapolation to this target. The proxy is uncalibrated and must not be called an AB magnitude color.

**Scientific hold:** Neither candidate has a defensible updated Gaussian significance. A representative-exposure-pair control null cannot simply be added to the frozen full-epoch PM covariance because of potential double counting of centroid/registration noise. The frozen and refit reductions use overlapping observations, so their difference cannot be assigned an independent-samples p-value. Neither control population is morphologically or kinematically verified as stationary.

## Next actions in scientific priority order

1. **1651050 F356W source association:** inspect the one registered late CAL exposure (`jw05893016008_05101_00001_nrcblong_cal.fits`), the centroid/forced-Gaussian residual and all available alternate F356W covering exposures; document why they were rejected. Cross-check against the frozen manuscript exposure selection.
2. **Frozen estimator replication:** trace original detector filenames, fit mode, centroid windows, registration weights and DQ handling. Reproduce the frozen estimates exactly before interpreting the refit discrepancy.
3. **Stationary-galaxy vetting:** use morphology, external photometry, isolation and consistency across independent epochs; distinguish field stars from background galaxies. Do not treat small measured motion alone as proof of stationarity.
4. **Full-epoch source-level null and covariance:** leave out each stationary control from the exact target estimator; bootstrap independent sources, preserve within-epoch dither correlations, and isolate incremental cross-filter variance already absent from nominal uncertainties.
5. **Extreme-color PSF check:** measure filter-dependent PSF centroid behavior with appropriate empirical PSFs and injected red-source SEDs; the F356W candidate is currently outside the control-color range.
6. Only after these tests derive a revised two-dimensional null probability with a defensible tail calibration; otherwise retain a provisional classification and explicitly state the limitation.

No number-density model prediction changes as a result of this astrometric audit. Any reported count of *astrometrically confirmed* objects must be reconsidered separately.

## Provisional response to Priyanka, not publication-ready

> Thank you for raising the cross-filter astrometric concern. We tested local reference-source displacements separately in the F444W/F410M and F444W/F356W JWST pairs, using leave-one-out affine registration on public CAL exposures. The retained reference-source scatters are 10–11 mas per coordinate for F444W/F410M (10 sources) and 12–16 mas for F444W/F356W (9 sources). The estimated covariance is sensitive to individual controls, and the extremely red F356W target lies outside the control instrumental-color range. We are verifying stationary galaxy controls and reconciling the original multi-exposure centroiding with independent fits before adopting an additional systematic covariance or quoting a revised proper-motion significance.

**Manuscript language until validation:** “An initial cross-filter astrometric diagnostic using held-out local reference sources shows residual scatter of order 10–16 mas per coordinate. These residuals mix centroid measurement noise, registration uncertainty, and possible source motion; they cannot yet be interpreted as an isolated chromatic PSF error. The proper-motion significance therefore remains provisional pending stationary-source and full-epoch covariance validation.”
