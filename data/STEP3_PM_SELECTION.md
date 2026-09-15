# Step3 magnitude-cut proper-motion selection

Input catalog: `candidates_step3_magnitude_cut_strict_125.csv` (**125 rows**).

The new catalog was positionally cross-matched to all candidate sets whose proper motion has already been audited. A **0.30 arcsec** positional radius is used for catalog-version/deblend bookkeeping; the PM measurement itself retains the much stricter source-association rules in the DQ-aware pipeline.

- Original Dave 86-candidate set: **60 distinct previous sources** are recovered, represented by **62 rows** in the new catalog because two previous positions are split into close catalog components.
- Ali 45-candidate set: **0** positional overlaps with this 125-row catalog.
- Dave new F444W 23-candidate set: all **23** exact candidate IDs are present. Two additional very close catalog components lie within 0.30 arcsec of that set; these are already included in the original-86 positional overlap bookkeeping.
- Union of previously audited rows in the current 125-row catalog: **83**.
- Remaining candidates requiring a new PM audit: **42**.

The 42-candidate run list is `candidates_step3_pm_unchecked_42.csv`.

For the new run, use the epoch-cached exhaustive DQ-aware JWST/NIRCam method: enumerate every legal independent epoch pair; reuse each measured epoch rather than recomputing it for every pair; reject adopted centroids with DO_NOT_USE, SATURATED, JUMP_DET, or OUTLIER DQ flags; keep the validated <=4-pixel local association bound; and report `INSUFFICIENT_DATA` only after all legal pairs have been exhausted.
