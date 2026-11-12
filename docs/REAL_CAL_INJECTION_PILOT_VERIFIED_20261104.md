# Verified REAL CAL injection/recovery pilot, 2026-11-04

**Run:** https://github.com/physicsbahare/proper_motion_86/actions/runs/41813715740 (both candidate jobs passed).
**Actual uploaded trial artifacts:** 387132 = 13313298020; 1651050 = 13313607046.
This analysis inserted sources into original JWST Stage-2 CAL cutouts and reran `measure_exposure_recentered`, `register_exposures` and `infer_proper_motion`. It is not a unit-test simulation.

| Actual measured pilot quantity | 387132: F444W/F410M | 1651050: F444W/F356W |
| --- | ---: | ---: |
| Source sites, stationary + moving each | 20 | 20 |
| Stationary recovered / attempted | 20/20 | 15/20 |
| Moving recovered / attempted | 20/20 | 14/20 |
| Stationary null 95th-percentile PM amplitude (mas/yr) | 60.733 | 77.920 |
| Original candidate observed PM amplitude (mas/yr) | 23.025 | 57.616 |
| Stationary trials exceeding observed motion | 10/20 | 2/15 |
| Conditional empirical exceedance (+1 correction) | 0.52381 | 0.18750 |
| Stationary PM east/north mean (mas/yr) | -11.809, +9.518 | +8.405, -4.946 |
| Moving recovery bias east/north (mas/yr) | -8.918, +8.965 | +17.681, -6.124 |

**Interpretation:** At these actual flux levels, the *single-exposure cross-filter* null sometimes produces false motions as large as the candidate measurements. This weakens claims based on a single image pair but **cannot invalidate the paper's multi-exposure results**. These values are not cross-filter-only systematics, Gaussian sigma values, or publication-ready revised significance.

**Critical limitations:** only **one original CAL exposure per epoch** was used in this pilot, while the frozen independent refits use 4+2 for 387132 and 4+1 for 1651050. The PSF is a high-S/N compact local-source stamp independently selected in each filter, not a validated ultra-red ePSF; background/crowding, control motions and true source Poisson noise are imperfectly modeled. The empirical exceedance probabilities are conditional on successful detection (1651050 lost 5/20 stationary and 6/20 moving sources). Twenty samples are not useful for testing 5-sigma tails.

**Immediate follow-up started:** full-epoch real CAL injection workflow https://github.com/physicsbahare/proper_motion_86/actions/runs/41814514100 runs the original target-usable exposure set with consistent stationary/moving injection positions. It must be validated before reporting a publication-ready final conclusion. No manuscript or main-branch values should be changed from this pilot. Independent photometric number-density predictions remain unchanged.
