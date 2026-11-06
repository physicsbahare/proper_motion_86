#!/usr/bin/env python3
"""Run galaxy-only full-epoch control nulls on real public JWST CAL detector data.

Candidate PMs from the paper are strictly frozen. No automatic reclassification.
"""
from __future__ import annotations
import argparse
import json
import sys
import traceback
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"first26"),str(ROOT/"scripts")]
from pm86.archive import get_pair_products,load_covering_cutouts,query_candidate_inventory
from pm86.galaxy_full_epoch import (choose_anchor,match_vetted,registered_solution,
                                     summarize_full_epoch_null)
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install as install_bounded_remote
from run_cross_filter_null import (EXPECTED_PAIRS,choose_frozen_pair,target_is_usable)


def plain(v):
    if isinstance(v,dict):return {str(k):plain(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [plain(x) for x in v]
    if isinstance(v,np.ndarray):return plain(v.tolist())
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,np.floating):return plain(float(v))
    if isinstance(v,np.bool_):return bool(v)
    if isinstance(v,float) and not np.isfinite(v):return None
    return v


def save(path,record):
    path.write_text(json.dumps(plain(record),indent=2,sort_keys=True)+"\n")


def run(cid,out):
    out.mkdir(parents=True,exist_ok=True)
    vetted=pd.read_csv(ROOT/"data/field_controls_priyanka/cosmos2025_photometric_galaxy_controls_17.csv")
    vetted=vetted[vetted.candidate_id.eq(cid)].copy()
    assert len(vetted) in (8,9) and vetted.catalog_id.is_unique
    catalog=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    frozen=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    row=catalog.loc[cid]
    orig=frozen.loc[cid]
    ra,dec=float(row.ra_deg),float(row.dec_deg)
    pair=choose_frozen_pair(
        query_candidate_inventory(cid,ra,dec),
        EXPECTED_PAIRS[cid],float(orig.baseline_days))
    early_products,late_products=get_pair_products(pair)
    early,ea=load_covering_cutouts(cid,ra,dec,"early",EXPECTED_PAIRS[cid][0],early_products)
    late,la=load_covering_cutouts(cid,ra,dec,"late",EXPECTED_PAIRS[cid][1],late_products)
    pd.DataFrame(ea+la).to_csv(out/"original_CAL_coverage.csv",index=False)
    if not early or not late:
        raise RuntimeError("No selected CAL exposures in one epoch")
    original_rows=[];accepted=[];maps={}
    for exp in early+late:
        m,c=measure_exposure_recentered(exp)
        original_rows.append(m)
        # Same astrometrically usable exposures as original independent
        # full-epoch candidate refit. Do not rescue failed later F356W frames.
        if target_is_usable(m):
            accepted.append(m)
            maps[exp.filename]=match_vetted(c,vetted)
    pd.DataFrame(original_rows).to_csv(out/"original_target_measurements_all_exposures.csv",index=False)
    measurements=pd.DataFrame(accepted)
    if measurements.empty or set(measurements.epoch_label)!={"early","late"}:
        raise RuntimeError("No original target-usable exposure sets for both epochs")
    # Document exact all target-usable CAL filenames and absence of extra late frames.
    anchor=choose_anchor(measurements,maps)
    n_total=int(len(measurements))
    all_control_rows=[]
    for filename,controls in maps.items():
        if not controls.empty:all_control_rows.append(controls)
    if all_control_rows:pd.concat(all_control_rows,ignore_index=True).to_csv(out/"identified_galaxies_per_exposure.csv",index=False)
    coverage=[{"galaxy_id":int(g.catalog_id),
               "catalog_ra_deg":float(g.ra_deg),"catalog_dec_deg":float(g.dec_deg),
               "n_target_usable_exposures":n_total,
               "n_detected_exposures":sum(int(g.catalog_id) in c.galaxy_id.values for c in maps.values()),
               "n_detected_early":sum(int(g.catalog_id) in maps[str(r.filename)].galaxy_id.values
                                      for _,r in measurements[measurements.epoch_label.eq("early")].iterrows()),
               "n_detected_late":sum(int(g.catalog_id) in maps[str(r.filename)].galaxy_id.values
                                     for _,r in measurements[measurements.epoch_label.eq("late")].iterrows())}
              for _,g in vetted.iterrows()]
    pd.DataFrame(coverage).to_csv(out/"galaxy_epoch_coverage.csv",index=False)
    reg,q,_,base=registered_solution(measurements,maps,anchor)
    if not reg.empty:reg.to_csv(out/"galaxy_only_registered_target_positions.csv",index=False)
    if not q.empty:q.to_csv(out/"galaxy_only_registration_quality.csv",index=False)
    runs=[];null_rows=[];target_rechecks=[]
    for galaxy in vetted.itertuples(index=False):
        gid=int(galaxy.catalog_id)
        r,q,source,check=registered_solution(measurements,maps,anchor,heldout=gid)
        check["galaxy_id"]=gid
        check["candidate_id"]=cid
        runs.append(check)
        if not source.empty:
            source["full_original_exposure_coverage"]=check.get("galaxy_appears_in_all_original_target_usable_exposures",False)
            null_rows.append(source)
        target_rechecks.append({"galaxy_id_omitted":gid,
            "n_exposures":check.get("n_registered",0),
            **check.get("target_refit_pm",{}),
            **({"galaxy_loo_delta_east_mas":check["heldout_galaxy_null"].get("delta_east_mas"),
                "galaxy_loo_delta_north_mas":check["heldout_galaxy_null"].get("delta_north_mas"),
                "galaxy_loo_status":check["heldout_galaxy_null"].get("status")}
               if "heldout_galaxy_null" in check else {})})
    if null_rows:pd.concat(null_rows,ignore_index=True).to_csv(out/"full_epoch_galaxy_positions_per_exposure.csv",index=False)
    pd.DataFrame(target_rechecks).to_csv(out/"leave_one_galaxy_out_target_refits.csv",index=False)
    save(out/"galaxy_full_epoch_loo_detail.json",runs)
    aggregate=summarize_full_epoch_null(runs)
    status={
        "candidate_id":cid,"scientific_status":"PRELIMINARY_ONLY_NEVER_ADOPT_PAPER_SIGNIFICANCE",
        "original_paper_pm_east_masyr":float(orig.mu_alpha_cosdec_masyr),
        "original_paper_pm_north_masyr":float(orig.mu_delta_masyr),
        "published_nominal_2d_quadrature_score":float(orig.significance_2d),
        "frozen_filters":EXPECTED_PAIRS[cid],"frozen_baseline_days":float(orig.baseline_days),
        "n_vetted_photometric_galaxy_candidates":len(vetted),
        "n_target_usable_exposures":n_total,
        "n_early_target_usable":int((measurements.epoch_label=="early").sum()),
        "n_late_target_usable":int((measurements.epoch_label=="late").sum()),
        "original_target_usable_CAL_files":measurements[["filename","filter","epoch_label","mjd","astrometric_snr"]].to_dict(orient="records"),
        "reference_early_frame":anchor,
        "all_galaxy_refit":base,"heldout_galaxy_null":aggregate,
        "n_target_loo_results":len(runs),
        "limitations":[
            "Photometric LePHARE TYPE=0 does not independently confirm stationary galaxies.",
            "No independent spectroscopic redshifts required or provided.",
            "Only original candidate-centroidable CAL exposures enter the full-epoch test.",
            "A single late-epoch CAL for 1651050 cannot calibrate within-epoch repeatability.",
            "Galaxy null combines intrinsic source centroid noise and registration, not isolated cross-filter PSF.",
            "Galaxy shapes/colors differ from the ultra-red point source, so results do not validate chromatic ePSF shifts.",
            "Control-galaxy nulls use original candidate exposure-error weights, not measured extended-galaxy centroid errors.",
            "Dither measurements are correlated and never independent null trials.",
            "Independent refit is not an exact replication of the frozen paper estimator.",
            "No new Gaussian sigma; no automatic PM reclassification; no number-density change."
        ]}
    save(out/"galaxy_full_epoch_summary.json",status)
    print(json.dumps(plain({"candidate_id":cid,"n_target_exposures":n_total,
                "early_late":[status["n_early_target_usable"],status["n_late_target_usable"]],
                "all_galaxy_refit":base.get("target_refit_pm"),
                "full_epoch_null":aggregate}),indent=2),flush=True)
    if aggregate["status"]!="PRELIMINARY_FULL_EPOCH_GALAXY_DIAGNOSTIC":
        raise RuntimeError("Full-epoch galaxy null insufficient; inspect artifact and do not claim calibration")
    return status


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidate",required=True,type=int,choices=EXPECTED_PAIRS.keys())
    ap.add_argument("--out",default="results_galaxy_full_epoch")
    args=ap.parse_args()
    install_bounded_remote()
    dest=Path(args.out)/str(args.candidate)
    try:run(args.candidate,dest)
    except Exception as e:
        dest.mkdir(parents=True,exist_ok=True)
        (dest/"fatal_error.txt").write_text(traceback.format_exc())
        print(f"FULL_EPOCH_GALAXY_AUDIT_FAILED:{type(e).__name__}:{e}",flush=True)
        raise


if __name__=="__main__":main()
