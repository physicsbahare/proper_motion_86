# COSMOS2025 source vetting for Priyanka's cross-filter PM question

The physical catalog is on the user's computer at:

\`/run/media/bahareh/Seagate Hub/COSMOSWeb_mastercatalog_v1.fits\`

This is **not** present on a GitHub runner or in this repository. Do not commit the multi-GB proprietary/local file. The script only reads the FITS.

1. On the user's Linux computer, \`git checkout review/cross-filter-systematics-20261103\` (or \`git pull\` if already on it).
2. \`pip install -r requirements.txt\`
3. Download actual previously measured field-control positions with \`gh run download 41716345445 --repo physicsbahare/proper_motion_86 --name dave-cross-filter-astrometry-audit --dir results_cross_filter_null\` (requires GitHub CLI login); or pass \`--download-controls\` to the script.
4. Run:

\`\`\`bash
python scripts/match_stationary_galaxy_controls.py \
    --catalog "/run/media/bahareh/Seagate Hub/COSMOSWeb_mastercatalog_v1.fits" \
    --controls-dir results_cross_filter_null \
    --out results_stationary_galaxy_controls
\`\`\`

The script converts published local astrometric tangent offsets to sky coordinates using actual candidate RA/Dec, matches to the FITS photometry extension within 0.25 arcsec, and marks whether LePHARE TYPE=0 supports a galaxy interpretation. It reads a zspec only if a recognized column actually exists. It refuses to assume that the LEPHARE extension rows align to photometry unless both have IDs and the IDs match in sampled sections. It exports all 21 preclipping control rows rather than silently selecting only already-robust residuals. Inspect \`cosmos2025_matched_field_controls.csv\` and \`summary.json\`.

**DO NOT call photometric classifications or low residuals "confirmed stationary galaxies".** Photometric morphology/SED class is insufficient, spec-z requires an actual credible value and source-specific verification; stationarity then needs an independent multi-epoch fit and careful control/source-color match. This code therefore reports **zero confirmed stationary galaxies by design** and **does not recalculate any PM significance**. The aim is to identify which controls are sufficiently promising for independent spectroscopic and astrometric verification. Use same exposure pair as baseline, then full epoch fit and source bootstrap.

Known limitations: COSMOS2025 v1 catalog may have different column names and no spectroscopic redshift. If the LePHARE extension row alignment cannot be authenticated, statuses remain UNKNOWN_OR_REJECTED; identify common IDs before accepting TYPE. NIRCam field controls may lie outside the COSMOS2025 footprint or be blended.
