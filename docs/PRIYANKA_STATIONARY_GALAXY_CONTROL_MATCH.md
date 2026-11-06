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
