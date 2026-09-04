# First26 final robustness audit

This audit perturbs S/N, forced-Gaussian width, axis ratio, positional uncertainty, reduced-chi2, centroid method, and searches every saved ranked pair attempt.
The <=4 pixel association envelope is held fixed to prevent the demonstrated 643980 neighbour-latching failure. DO_NOT_USE and SATURATED pixels are never relaxed.

## Verdict counts
- `ROBUST_INSUFFICIENT_ACROSS_TESTED_GATES`: 26

Candidates requiring fresh targeted rerun: []

No candidate acquires two-epoch astrometric availability under any tested plausible gate perturbation. The current INSUFFICIENT_DATA classifications are therefore robust to these threshold/centroid/pair-selection choices.
