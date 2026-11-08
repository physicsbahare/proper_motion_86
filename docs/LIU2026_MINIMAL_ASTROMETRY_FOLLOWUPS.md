# Minimal follow-ups inspired by Liu et al. (2026): two Dave proper-motion candidates

**Status:** planning only, 2026-11-08. No new science tests have been run on this branch and no manuscript values were changed.  
**Base branch:** `review/cross-filter-systematics-20261103`.  
**Science scope:** *only* COSMOS2025 IDs **387132** (F444W/F410M) and **1651050** (F444W/F356W).  
**Reference:** Liu et al. (2026), *The cold substellar frontier in JWST Extragalactic surveys*, arXiv:2611.08346v1, especially Sections 3.1, 3.2.1–3.2.3. Their 50-F444W-exposure / five-epoch parallax solution **is not directly transplantable** to these two two-epoch, cross-filter candidates.

## 1. What is already DONE — do not duplicate

| Existing check | Actual status / source | What it does **not** establish |
|---|---|---|
| Detector-level original CAL and JWST-L2 parallax sensitivity | [Publication audit](https://github.com/physicsbahare/proper_motion_86/blob/publication-pm-parallax-audit-20261008/docs/PUBLICATION_ASTROMETRY_FINAL_FINDINGS.md), [PR #1](https://github.com/physicsbahare/proper_motion_86/pull/1) | Neither object's trigonometric parallax was directly measured. 1651050's intrinsic PM remains degenerate with assumed distance. |
| Empirical PSF (ePSFBuilder) versus Gaussian on original images | Same publication audit: 1651050 F444W ePSF–Gaussian north offsets -56.3 and +22.1 mas in different early images; F356W result about (-10.7,-16.8) mas | This **did not** validate a stable color-matched ePSF centroid; no reason to repeat the same generic ePSF fit unchanged. |
| Representative-exposure cross-filter held-out controls | [First cross-filter results](DAVE_CROSS_FILTER_FIRST_RESULTS_20261103.md), [covariance stability](PRIYANKA_COVARIANCE_STABILITY_VERIFIED_20261103.md) | Mixed stars/galaxies, very small N, no independent cross-filter PSF bias or Gaussian 5-sigma significance. |
| COSMOS2025/LePHARE control vetting | [Galaxy match record](PRIYANKA_STATIONARY_GALAXY_CONTROL_MATCH.md): 9 and 8 **photometric galaxy candidates**, not secure spectroscopic stationary galaxies | Photometric TYPE=0 alone is not independent astrometric stationarity or a color match to an ultra-red dwarf. |
| Full-epoch held-out galaxy tests | [Galaxy-only](PRIYANKA_FULL_EPOCH_GALAXY_REAL_RESULTS_20261106.md) and [original-reference-frame](PRIYANKA_ORIGINAL_FRAME_GALAXY_RESULTS_20261106.md) | Full six-CAL held-out null was unavailable for 387132 with galaxy-only registration; original-reference registration yielded 6 and 7 full-coverage held-out galaxies. Small-N null is not a 5-sigma calibration. |
| Real-CAL stationary and moving source injection, representative pair | [Verified 20-site pilot](REAL_CAL_INJECTION_PILOT_VERIFIED_20261104.md), [run 41813715740](https://github.com/physicsbahare/proper_motion_86/actions/runs/41813715740) | Single-pair injection is *not* the original multi-exposure estimator. |
| **Full-epoch real-CAL stationary/moving injection** | **Already executed**, [run 41814514100](https://github.com/physicsbahare/proper_motion_86/actions/runs/41814514100). Inspect existing uploaded artifacts before rerunning. 32 independent sky sites per target. For **387132**, recovered 30/32 stationary and 31/32 moving, zero of 30 recovered nulls exceeded observed PM; conditional p=(0+1)/(30+1)=**0.0323**. For **1651050**, recovered 24/32 stationary and 24/32 moving; 3/24 recovered nulls exceeded the observed PM; conditional p=(3+1)/(24+1)=**0.16**. | These are small-sample **conditional diagnostic tail fractions**, NOT Gaussian sigmas, formal final false-positive probabilities, or a publication-ready revision. Injections used local compact-source PSF stamps **not validated at the ultra-red target SED**; common WCS error and injection noise are incompletely modeled. Successful workflow means the trial ran, not that scientific validity is demonstrated. |

Original frozen PM vector (east,north): **387132 (-22.900,+2.400)** mas/yr; **1651050 (-55.230,+16.408)** mas/yr. Full-epoch independent original-frame refits give (-19.983,+2.056) and (-61.741,+17.150) mas/yr and label **both AMBIGUOUS** under the refit's internal criteria. These are not independent observations of the frozen results. Do not overwrite published outputs on `main`.

**Critical target-specific limitations**
- **387132:** four usable early F444W and two usable late F410M CAL exposures. Conditional JWST-parallax correction is small at the *assumed* 19.1 pc photometric distance. Relative motion is more robust; full-epoch galaxy null N is nevertheless small.
- **1651050:** four usable early F444W but **only one usable late F356W CAL**; F356W S/N about 5, nearby DQ/background structure, centroid-method sensitivity, and measured displacement can be substantially altered by **assumed** parallax/distance. Previous MIRI/HST searches **did not validate a third centroid**. Do not claim independently confirmed intrinsic PM, measured parallax, or a new 5-sigma detection.

## 2. Small, genuinely incremental future work — in priority order

### P0. Reconcile and preserve the **existing** full-epoch results (documentation, not a new test)
- [ ] Inspect the **existing** `run 41814514100` `trials.csv`, `result.json`, `original_exposure_quality.csv`, `site_selection_diagnostics.json` artifacts; save a short **one-page verified comparison** with counts, selection failures, null and moving-recovery bias, and exposure names for both targets. Logs already verify 32-site runs (see above). **Do not call a rerun 'new evidence'**.
- [ ] Reconcile frozen paper estimator versus current independent refit: exact CAL filenames, accepted/rejected exposures, centroid mode, association, DQ masking, registration and epoch weights; explain any mismatch without changing either stored measurement. Reuse `scripts/run_cross_filter_null.py`, `src/pm86/refit_reconciliation.py`, `src/pm86/galaxy_full_epoch.py` and the published diagnostics.
- [ ] **1651050 only:** inspect the single retained late `jw05893016008_05101_00001_nrcblong_cal.fits` F356W detection, along with the *rejection reasons* for its other covering exposures, SCI/ERR/DQ, nearby neighbors, and 2D Gaussian vs center-of-mass vs available ePSF centroids. Reuse existing PR #1 pixel panels first.

**Stop/gate:** If the sole F356W source association is not credible or the frozen estimator cannot be reproduced, stop astrometric 'confirmation' claims and explicitly label 1651050 inconclusive. Do not expand to a survey-wide reprocessing.

### P1. Minimal, **target-centered** PSF / cross-filter sanity test (Liu §3.2.1 adapted)
- [ ] Only after P0: on the *same* accepted CALs, compare the **current Gaussian centroid**, the **already attempted local empirical ePSF**, and **one color/SED-weighted STPSF/WebbPSF option** appropriate to each target/filter, if reliable inputs are available. For each, inspect source and nearby control residual stamps and record east/north centroid **difference in mas** for *each individual exposure*, plus sensitivity to one reasonable fitting window. Do not build a broad 50-image hybrid-PSF pipeline.
- [ ] Briefly visually vet the most influential existing galaxy controls, especially 1651050 catalog IDs **726373** (large residual) and **726167** (larger catalog-match offset), for blends, isolation, DQ and morphology. Reuse the existing 6 and 7 *fully covered* held-out-galaxy tests from the original reference frame. Do not claim galaxy colors alone calibrate ultra-red target chromatic shifts.
- [ ] Record whether any systematic direction/amplitude can explain the candidate's apparent displacement; report **method-to-method spread**, not an invented Gaussian PSF error floor.

**Stop/gate:** If color-aware PSF fitting cannot be independently validated or F356W remains too weak, report that explicitly. This is an honest finite result, not a reason for an open-ended new model grid.

### P2. Only **conditional** targeted injection extension (Liu §3.2.2–3.2.3 adapted)
- [ ] First inspect results of the existing full-epoch 32-site run. **Do not** rerun the representative 20-site pilot or increase trials merely to chase a nominal sigma.
- [ ] If P1 yields a validated, materially different PSF and a stable real-target centroid, modify the existing `scripts/real_cal_injection_pilot.py --full-epoch` to bracket a **stationary** and a **known moving** target at each filter's actual flux/SNR, using the exact accepted CALs. Include realistic pixel-phase/PSF choices; preserve failed recoveries and compute null tails *conditional on recovery* and completeness separately. Perform a modest reproducible repeat only if scientifically justified, e.g. 32–64 independent sky sites first, subject to real overlap and finite compute resources.
- [ ] Avoid pseudo-replication: multiple dithers of one fake sky location are **one trial**, not independent controls. Prefer matched method comparisons over an unvalidated >500-trial or >100,000-bootstrap exercise.
- [ ] If a sufficiently validated source-level null exists, report its empirical tail resolution and check **additional**, non-double-counted uncertainty beyond the current registration term. Do **not** simply add full held-out galaxy covariance to nominal PM errors or call a 2D quadrature score 'sigma'.

**Stop/gate:** If 1651050's stationary trials still recover false motion near the observed displacement or moving recovery is biased/incomplete, keep astrometry *inconclusive*; don't force a detection classification.

### P3. Conclude in a **short, decision-ready result**
- [ ] Write one compact comparison table (candidate, frozen motion, full-epoch refit, number of usable exposure frames, vetted galaxy-null N and scatter, injection detection/recovery fractions, centroid/PSF sensitivity, parallax caveat).
- [ ] Distinguish *apparent two-epoch displacement* from *intrinsic proper motion* and from *measured trigonometric parallax*. Give recommended wording for each target, stating exact evidence and any unresolved limitations. Do not change the paper's frozen tables or the photometric number-density estimate without explicit later scientific review.

## 3. Explicitly OUT OF SCOPE (avoid project divergence)

- **No five-parameter joint PM+parallax fit for 1651050 from only two independent epochs.** This is underconstrained; a parallax-only conditional fit is not a geometric distance.
- **No 850-reference-source, five-epoch Liu-style global connected-field solution or 20-arcsec block bootstrap** imported blindly into fields with only 6–9 vetted controls. Small-N source-level LOO/jackknife is already available; do not pretend a block bootstrap creates independent data.
- No new full COSMOS-Web/COSMOS2025 crossmatch, re-ranking all 125 F444W-only objects, PM-footprint model, MIRI pseudo-detection, atmosphere-grid refit, or re-running the same ePSF/injections with identical assumptions.
- No '5σ' conclusion from 24–30 recovered null trials, no independent full-epoch chromatic systematic deduced from ordinary galaxy scatter, no extrapolation beyond the color range of the controls.
- No blanket merger of PR #1 or #2 into `main`, and no publication-data modification as part of this *plan* branch.

## 4. Existing entry points and references

- Primary science-code base: branch `review/cross-filter-systematics-20261103` ([PR #2](https://github.com/physicsbahare/proper_motion_86/pull/2)); this planning branch starts **from that branch, not from `main`**.
- Independent JWST-parallax/ePSF audits: branch `publication-pm-parallax-audit-20261008` ([PR #1](https://github.com/physicsbahare/proper_motion_86/pull/1)).
- Code already available: `scripts/run_galaxy_null_original_frame.py`, `scripts/run_full_epoch_galaxy_null.py`, `scripts/evaluate_galaxy_only_controls.py`, `scripts/real_cal_injection_pilot.py`, `src/pm86/cross_filter_injection.py`, `src/pm86/galaxy_full_epoch.py`.
- Existing real-CAL artifacts: [full-epoch injection run 41814514100](https://github.com/physicsbahare/proper_motion_86/actions/runs/41814514100), [full-epoch galaxy-only run 42182049157](https://github.com/physicsbahare/proper_motion_86/actions/runs/42182049157), [original-reference galaxy-null run 42183463607](https://github.com/physicsbahare/proper_motion_86/actions/runs/42183463607).
- Relevant Liu methods to adapt *selectively*: PSF residual/centroid sensitivity (§3.2.1), coherent track-based injection/recovery and uncertainty propagation (§3.2.2), cross-band consistency and explicit null injections (§3.2.3). The science aim is **robustness of the existing two-epoch displacement**, not a replica of the JADES high-S/N parallax result.

**Definition of done:** A reproducible and concise record showing what was already tested; whether one small additional PSF/injection check actually changes the evidential status of either target; and an appropriately cautious, specific scientific conclusion for both. If neither new check is feasible, a rigorously documented unresolved limitation is a valid endpoint.
