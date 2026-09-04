# Candidate 643980 forensic proper-motion audit

## Question
Why did candidate 643980 change from `CONSISTENT_WITH_ZERO / PM_MEASURED` in the earlier first26 run to `INSUFFICIENT_DATA` in the current pipeline?

## Result
The earlier PM measurement is not a valid measurement of candidate 643980. The old motion-aware recentering step searched out to 12 NIRCam LW pixels and latched onto unrelated detections roughly 10--11 pixels from the catalog-predicted target position. Those distant detections were then treated as the target in both epochs, producing a formally measurable same-filter F444W displacement.

The current 4-pixel association radius correctly prevents that failure mode. At the actual catalog-predicted position, the source is not securely detected/centroidable in both independent epochs, so the defensible astrometric result is `INSUFFICIENT_DATA`.

## Evidence from the earlier PM-producing run (commit 5426f33)
The reported same-filter F444W result was:

- baseline = 717.322 d = 1.9639 yr
- mu_alpha*cos(dec) = -24.56 +/- 16.32 mas/yr
- mu_delta = +19.48 +/- 18.56 mas/yr
- 2D significance = 1.83 sigma
- classification = `CONSISTENT_WITH_ZERO`

However, the exposure measurements used recentered detections far from the target prediction:

- early `jw01837003016_08201_00002_nrcblong_cal.fits`: recenter offset = 10.173 pix, clean S/N = 5.53
- late `jw05398275001_05201_00002_nrcblong_cal.fits`: recenter offset = 10.182 pix, clean S/N = 7.69
- late `jw05398275001_05201_00003_nrcblong_cal.fits`: recenter offset = 11.274 pix, clean S/N = 9.62

These offsets correspond to about 0.64--0.71 arcsec at the ~62.6 mas/pixel NIRCam LW scale, far outside a conservative source-association envelope for an autonomous nearest-source rescue. The old code explicitly allowed `SEARCH_RADIUS_PIX = 12.0`.

Several of the original catalog-position measurements in these exposures also carried problematic DQ information (including `DO_NOT_USE` / `JUMP_DET`), while the positive high-S/N measurement appeared only after shifting to the distant local detection. This is another indication that the PM solution was created by source mis-association rather than by the target itself.

## Current-pipeline check
The current recentering code uses `SEARCH_RADIUS_PIX = 4.0` (~0.25 arcsec) and states explicitly that larger offsets require cross-epoch trajectory linking rather than independent nearest-neighbour recentering.

For the best current F444W--F444W pair:

- covering exposures: early = 4, late = 3
- accepted astrometric detections: early = 0, late = 0

Nine ranked independent NIRCam pairs were tested and all failed the secure two-epoch centroid requirement.

## Final scientific interpretation
- The earlier numerical PM should be discarded; it was produced by a distant-neighbour association.
- There is no evidence here for significant proper motion of candidate 643980.
- There is also no reliable zero-PM measurement for the target itself.
- Final PM status: `INSUFFICIENT_DATA`.
- This regression is explained by a scientifically necessary association fix, not by random pipeline instability.

The SED/AIC galaxy-versus-brown-dwarf classification is independent of this astrometric conclusion.