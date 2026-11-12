# COSMOS2025 control-source vetting for Priyanka's astrometry comment

## Ready in the review branch

The 21 original matched field-control records are **committed**, with original measured residuals and sky-relative positions:

- `data/field_controls_priyanka/387132/leave_one_out_residuals.csv` (12 sources, 10 retained for covariance)
- `data/field_controls_priyanka/1651050/leave_one_out_residuals.csv` (9 sources, 9 retained)
- Provenance: GitHub Actions original detector-frame audit run [41716345445](https://github.com/physicsbahare/proper_motion_86/actions/runs/41716345445).

The catalog is **only on the user's Seagate drive**:
`/run/media/bahareh/Seagate Hub/COSMOSWeb_mastercatalog_v1.fits`.
GitHub-hosted CI cannot read that path. The public CDS/TAP backup failed on 2026-11-06 (run 42168384835, remote connection timeouts), so **do not claim the FITS cross-match completed remotely**.

## On the Linux machine with the catalog mounted

```bash
git fetch origin
git switch review/cross-filter-systematics-20261103
git pull
python -m pip install numpy pandas scipy astropy
python scripts/match_stationary_galaxy_controls.py
```

The default file path already matches the Seagate location, and the checked-in 21 control records are used automatically. Neither GitHub CLI nor any data download is needed.

Outputs:
- `results_stationary_galaxy_controls/cosmos2025_matched_field_controls.csv`
- `results_stationary_galaxy_controls/summary.json`

The cross-matcher reconstructs each control's actual sky position from measured local tangent offsets, matches within 0.25 arcsec to the COSMOS2025 photometry extension, checks LePHARE TYPE=0 using the official COSMOS2025 same-row-order convention when the full extension lengths agree (otherwise a unique shared-ID join when possible), and records any available spectroscopic redshift, blend flags and angular residuals. Every control is preserved, including originally clipped outliers. It **does not overwrite** the paper's PM values.

## How to interpret

- `galaxy_photometric_candidate=True`: supported by available LePHARE classification and uncontaminated blend flag; **not spectroscopically confirmed**.
- `spectroscopic_redshift_support=True`: a positive spec-z column was found; requires validating its provenance and quality before claiming spectroscopic confirmation.
- `stationarity_status=NOT_INDEPENDENTLY_CONFIRMED`: no independent, full-epoch stationarity verification has been done. The script deliberately never calls any object confirmed stationary.
- `catalog_position_match=False`: object not reliably matched within 0.25 arcsec.
- Distinct star/galaxy classification, filter-specific PSF behavior, intrinsic PM, common registration errors, and centroid quality still require inspection.

The first user run matched 19/21 controls to the photometric table but classified none because the initial script refused to read the ID-less LePHARE extension. The revised script uses the catalog's documented row-order convention, and reports its join method and LePHARE rows explicitly. The real FITS crossmatch requires a new run after pulling the fix. Do not infer the galaxy count until that run completes.

The Gaussian significance is **not** recalculated by this catalog crossmatch, and no change is made to the photometric number-density prediction.

## Validation

The synthetic-FITS regression tests passed in [workflow 42169048841](https://github.com/physicsbahare/proper_motion_86/actions/runs/42169048841). This validates the matching logic, **not** the real catalog classifications.


## Real COSMOS2025 crossmatch and galaxy-only LOO results (2026-11-06)

**Provenance:** User locally ran `scripts/match_stationary_galaxy_controls.py` on their 784,016-row COSMOS2025 master catalog, then `scripts/evaluate_galaxy_only_controls.py`. The two corresponding JSON/CSV result sets were provided in the research conversation. The original FITS remains on the user's machine, and neither the original candidate PMs nor their significances have been modified.

- Photometry positions matched at 0.25 arcsec: 19/21 controls, 10/12 near 387132, 9/9 near 1651050.
- LePHARE row-order alignment: `COSMOS2025_OFFICIAL_ROW_ORDER_EQUAL_LENGTH`; both HDUs have 784,016 rows.
- LePHARE TYPE=0 photometric galaxy candidates: 9 near 387132 and 8 near 1651050. Neither field has spectroscopically confirmed stationary sources from this check.
- Galaxy-only leave-one-out affine registration, original representative detector exposure pair, all 9/8 successfully fitted:

| Diagnostic | 387132 (F444W/F410M) | 1651050 (F444W/F356W) |
| --- | ---: | ---: |
| Galaxy candidate controls | 9 | 8 |
| Galaxy-only median east/north residual (mas) | -0.247 / +3.072 | +6.196 / +0.231 |
| Galaxy-only east/north scatter (mas) | 10.141 / 10.911 | 12.927 / 13.618 |
| Original mixed-control scatter (mas) | 9.979 / 10.849 | 15.580 / 12.233 |
| Galaxy-only covariance (mas²) | [[102.837,-2.709],[-2.709,119.056]] | [[167.098,72.352],[72.352,185.446]] |

**Interpretation:** Excluding likely non-galaxies does **not** radically change the observed reference scatter. This is reassuring about the *selection sensitivity* of the pair-level residuals, not evidence of a validated chromatic systematic or confirmed motion significance. Small N and shared registration remain limiting; F356W target 1651050's instrumental color is outside the reference color range.

**Source-level warnings:** LePHARE `zfinal=-99` occurs for galaxy-labeled control COSMOS2025 ID 333691 near 387132 and is **not** a physical redshift. One of the 387132 leave-one-out registrations uses only six fitting galaxies after clipping; residual covariance may be source sensitive. The 1651050 galaxy control matched to COSMOS2025 ID 726167 at 0.1236 arcsec, less secure than the other controls (mostly <0.04 arcsec), and merits visual verification. No independent morphology, robust spectral redshift, all-epoch stationary fit, or color-matched ePSF validation has yet been performed.

**Decision:** The current data support a statement that *nearby photometrically selected galaxy controls show roughly 10–14 mas positional scatter per coordinate across the relevant filter pairs*. It would be inaccurate to claim these values are an additional, independently measured PSF bias, to add the entire covariance to existing registration errors, or to quote new 5σ detection values. The 128-position real CAL injection-recovery results are a separate serious quality check, especially for 1651050. Manuscript PM measurements remain frozen until full-epoch validation.
