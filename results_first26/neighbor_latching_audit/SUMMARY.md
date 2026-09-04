# First26 historical neighbour-latching audit

Historical reference: `5426f33b0f4e92b36e4278bb0943035dd902a518`
Association limit tested: > 4.0 pix

## Verdict counts
- `NO_LARGE_OFFSET_LATCHING_EVIDENCE`: 24
- `HISTORICAL_LARGE_OFFSET_RECENTER_PRESENT_NOT_ASTROMETRIC`: 1
- `CONFIRMED_OLD_PM_CONTAMINATED_BY_NEIGHBOR_LATCHING`: 1

Confirmed old PM contaminated by large-offset neighbour latching: [643980]
Other candidates with historical >4-pix recentering to review: [575676]

A historical large-offset recenter is not itself a PM detection. The critical condition requires an old PM measurement plus astrometry-like >4-pix associations in both epochs. Current results remain the science result unless this audit identifies a separate reproducible issue.
