# Data provenance

`candidates_86_input.csv` is a text export of the original 86-row candidate catalog `step2_final_capotauro_candidates.xlsx`. Only catalog quantities needed by this audit are retained here: source ID, RA, Dec, and F444W photometry/flux-ratio fields.

`candidates_45.csv` is derived mechanically from those 86 rows by excluding:

1. candidate `282040`, which was audited separately as the reference mover; and
2. the 40 IDs in `validated_artifact_ids.txt`, whose F444W artifact behavior was independently checked exposure-by-exposure and then refined with full JWST gWCS positioning.

No external proper-motion or stationary/moving classification is included in `candidates_45.csv`.
