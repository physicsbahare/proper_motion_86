# Remote run history and validity

This file records pipeline failures separately from scientific `INSUFFICIENT_DATA` outcomes so a green GitHub Actions badge cannot be mistaken for a valid astrometric result.

## Run 1 — Actions run 33618995907

The workflow itself completed, but 23 candidates with potential multi-epoch NIRCam coverage failed at MAST product retrieval because `astroquery==0.4.11` does not accept the `batch_size` keyword in `Observations.get_product_list`. The 22 candidates rejected earlier by epoch selection were unaffected. The aggregated `45/45 INSUFFICIENT_DATA` table from this run must not be interpreted scientifically.

## Run 2 — Actions run 33622532095

The `batch_size` issue was fixed, but the same 23 candidates then exposed a second compatibility bug: `Observations.get_product_list` internally joins observation IDs as strings, while the pipeline passed `numpy.int64` values. The error was `sequence item 0: expected str instance, numpy.int64 found`. Again, the 22 candidates with no independent NIRCam epoch were unaffected, but the aggregate `45/45 INSUFFICIENT_DATA` result is not a valid final PM classification.

## Safeguard added after Run 2

The pipeline now converts MAST observation IDs explicitly to strings and the aggregation workflow performs candidate-level QC. Any remaining `pm_status == ERROR` rows are committed for debugging but make the workflow fail visibly instead of producing a misleading green run. The same aggregation step also runs the 60 mas centroid-method/frame sensitivity audit.
