#!/usr/bin/env python3
"""Cross-filter full-epoch null using original mixed-field registration,
but evaluating only held-out COSMOS2025 photometric galaxy controls.
This is a DIFFERENT, separately labeled audit from the galaxy-only fit.
"""
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
from pm86.galaxy_full_epoch import (
    match_vetted,original_registration_galaxy_null,summarize_full_epoch_null)
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install as install_remote
from run_cross_filter_null import EXPECTED_PAIRS,choose_frozen_pair,target_is_usable

def safe(x):
    if isinstance(x,dict):return {str(k):safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [safe(a) for a in x]
    if isinstance(x,np.integer):return int(x)
    if isinstance(x,np.bool_):return bool(x)
    if isinstance(x,np.floating):return safe(float(x))
    if isinstance(x,float) and not np.isfinite(x):return None
    return x

def dump(path,obj):path.write_text(json.dumps(safe(obj),indent=2,sort_keys=True)+"\n")

def run(cid,out):
    out.mkdir(parents=True,exist_ok=True)
    vetted=pd.read_csv(ROOT/"data/field_controls_priyanka/cosmos2025_photometric_galaxy_controls_17.csv")
    vetted=vetted[vetted.candidate_id.eq(cid)].copy()
    catalog=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    orig=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id").loc[cid]
    ra,dec=float(catalog.loc[cid].ra_deg),float(catalog.loc[cid].dec_deg)
    inv=query_candidate_inventory(cid,ra,dec)
    pair=choose_frozen_pair(inv,EXPECTED_PAIRS[cid],float(orig.baseline_days))
    ep,lp=get_pair_products(pair)
    early,ea=load_covering_cutouts(cid,ra,dec,"early",EXPECTED_PAIRS[cid][0],ep)
    late,la=load_covering_cutouts(cid,ra,dec,"late",EXPECTED_PAIRS[cid][1],lp)
    pd.DataFrame(ea+la).to_csv(out/"CAL_coverage.csv",index=False)
    records=[];raw={};galaxy={}
    for exp in early+late:
        m,c=measure_exposure_recentered(exp)
        records.append(m)
        if target_is_usable(m):
            raw[exp.filename]=c
            galaxy[exp.filename]=match_vetted(c,vetted)
    pd.DataFrame(records).to_csv(out/"original_exposure_quality.csv",index=False)
    rows=pd.DataFrame([r for r in records if target_is_usable(r)])
    if rows.empty or set(rows.epoch_label)!={"early","late"}:
        raise RuntimeError("Original candidate lacks two usable detector epochs")
    if not raw or any(c.empty for c in raw.values()):
        raise RuntimeError("Original target-usable exposure missing all raw controls")
    rawcat=pd.concat(raw.values(),ignore_index=True)
    rawcat.to_csv(out/"original_field_controls_each_exposure.csv",index=False)
    galcat=pd.concat([c for c in galaxy.values() if not c.empty],ignore_index=True)
    galcat.to_csv(out/"galaxy_controls_each_exposure.csv",index=False)
    coverage=[{"galaxy_id":int(z.catalog_id),
               "n_detected":sum(int(z.catalog_id) in c.galaxy_id.values for c in galaxy.values()),
               "n_target_usable":len(rows)}
               for _,z in vetted.iterrows()]
    pd.DataFrame(coverage).to_csv(out/"galaxy_coverage.csv",index=False)
    base,q,_,b=original_registration_galaxy_null(rows,raw,galaxy)
    if not base.empty:base.to_csv(out/"original_frame_target_positions.csv",index=False)
    if not q.empty:q.to_csv(out/"original_frame_registration_quality.csv",index=False)
    results=[];positions=[];refits=[]
    for g in vetted.itertuples(index=False):
        gid=int(g.catalog_id)
        reg,qual,source,summary=original_registration_galaxy_null(rows,raw,galaxy,gid)
        summary.update({"candidate_id":cid,"galaxy_id":gid})
        results.append(summary)
        if not source.empty:
            positions.append(source)
        refits.append({"galaxy_id_omitted":gid,
                       "n_registered":summary.get("n_registered_exposures"),
                       **summary.get("target_refit_pm",{}),
                       "galaxy_null_status":summary.get("heldout_galaxy_null",{}).get("status"),
                       "galaxy_null_delta_east_mas":summary.get("heldout_galaxy_null",{}).get("delta_east_mas"),
                       "galaxy_null_delta_north_mas":summary.get("heldout_galaxy_null",{}).get("delta_north_mas")})
    if positions:pd.concat(positions,ignore_index=True).to_csv(out/"heldout_galaxy_positions_each_exposure.csv",index=False)
    pd.DataFrame(refits).to_csv(out/"galaxy_loo_original_frame_target_refits.csv",index=False)
    dump(out/"galaxy_loo_original_frame_all_details.json",results)
    aggregate=summarize_full_epoch_null(results)
    summary={"candidate_id":cid,"status":"ORIGINAL_MIXED_FIELD_REGISTRATION_VETTED_GALAXY_NULL",
             "source_selection":"photometric LePHARE type0 candidate, not independent stationary confirmation",
             "registration":"original mixed field controls, each galaxy excluded from every fit used to measure it",
             "n_target_usable_exposures":len(rows),
             "n_early":int((rows.epoch_label=="early").sum()),
             "n_late":int((rows.epoch_label=="late").sum()),
             "all_field_original_frame_target_refit":b.get("target_refit_pm"),
             "heldout_galaxy_full_epoch_null":aggregate,
             "frozen_mu_east_north_masyr":[float(orig.mu_alpha_cosdec_masyr),float(orig.mu_delta_masyr)],
             "NOT_a_Gaussian_significance":True,
             "NOT_a_galaxy_only_registration":True,
             "NOT_stationarity_confirmed":True,
             "NOT_an_additional_PSf_covariance":True,
             "additional_limitations":[
                 "Mixed registration includes potentially moving foreground field stars.",
                 "Photometric galaxy candidates are not individually confirmed stationary.",
                 "Original source selection is approximate and frozen manuscript estimator has not been replicated exactly.",
                 "Galaxy centroids and WCS registration scatter are mixed; no distinct chromatic PSF systematic has been isolated.",
                 "1651050 has a single usable late F356W exposure and no internal repeatability check.",
                 "All multi-exposure dither measurements share galaxy identities and correlations; source sample N, not number of CALs, sets null resolution."]}
    dump(out/"original_frame_full_epoch_galaxy_null_summary.json",summary)
    print(json.dumps(safe({"candidate_id":cid,"n_early_late":[summary["n_early"],summary["n_late"]],
            "original_registration_target":b.get("target_refit_pm"),"galaxy_null":aggregate}),indent=2),flush=True)
    if aggregate["status"]!="PRELIMINARY_FULL_EPOCH_GALAXY_DIAGNOSTIC":
        raise RuntimeError("No adequately complete galaxy null in original reference frame")
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--candidate",required=True,type=int,choices=list(EXPECTED_PAIRS))
    p.add_argument("--out",default="results_galaxy_original_frame_full_epoch")
    a=p.parse_args()
    install_remote()
    target=Path(a.out)/str(a.candidate)
    try:run(a.candidate,target)
    except Exception:
        target.mkdir(parents=True,exist_ok=True)
        (target/"fatal_error.txt").write_text(traceback.format_exc())
        raise

if __name__=="__main__":main()
