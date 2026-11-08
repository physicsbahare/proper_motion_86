#!/usr/bin/env python3
"""1651050 targeted original-estimator F444W frame leave-out robustness.

Restores original recentered target measurements on all real NIRCam CAL frames
and repeats exact pm86 local-control affine registration and PM estimation.
Drop each of the four originally accepted F444W early detector exposures in
turn; also jointly exclude the two early images with a >0 DQ hardbad *3x3*
core at the accepted target coordinate (as measured by the prior fixed-model
STPSF check). The sole originally accepted late F356W is always preserved.

Frame-removal systematic sensitivity only. No new PSF-based astrometry,
uncertainty significance or paper revisions.
"""
from __future__ import annotations
import argparse,json,traceback,sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"first26"),str(ROOT/"scripts")]
from pm86.archive import query_candidate_inventory,get_pair_products,load_covering_cutouts
from pm86.astrometry import register_exposures,infer_proper_motion
from run_cross_filter_null import choose_frozen_pair,EXPECTED_PAIRS,target_is_usable
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install
from audit_1651050_f356w_joint import safe_out
import pm86.measurement as measure

CID=1651050
def inspect_core(ex,row):
    xc,yc=int(round(float(row["x_2dg"]))),int(round(float(row["y_2dg"])))
    _,astrom_bad=measure.masks(
        np.asarray(ex.sci,float),np.asarray(ex.err,float),
        np.asarray(ex.dq,np.uint32))
    return int(np.count_nonzero(astrom_bad[yc-1:yc+2,xc-1:xc+2]))

def perform(keep,measurements,controls):
    sub=pd.DataFrame([measurements[e.filename] for e in keep])
    c=pd.concat([controls[e.filename] for e in keep],ignore_index=True)
    reg,q,matched=register_exposures(sub,c)
    if reg.empty or set(reg.epoch_label)!={"early","late"}:
        raise RuntimeError("registered target has missing epoch")
    r,pairs=infer_proper_motion(reg,q,"CROSS_FILTER_JWST",matched)
    keys=["pm_status","classification","baseline_days","mu_alpha_cosdec_masyr",
          "mu_delta_masyr","sigma_mu_alpha_cosdec_masyr","sigma_mu_delta_masyr",
          "significance_2d","centroid_method_disagreement","centroid_method_delta_masyr",
          "cross_filter_color_systematic_flag","pm_covariance_status","pm_ra_dec_rho"]
    result={k:safe_out(r.get(k)) for k in keys}
    result.update(n_registered=int(len(reg)),n_early=int((reg.epoch_label=="early").sum()),
                  n_late=int((reg.epoch_label=="late").sum()),
                  n_used_controls=int(q["n_used"].fillna(0).sum()) if len(q) and "n_used" in q else None)
    return result,reg,q,pairs

def run(out):
    catalog=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    frozen=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    ra,dec=float(catalog.loc[CID].ra_deg),float(catalog.loc[CID].dec_deg)
    pair=choose_frozen_pair(query_candidate_inventory(CID,ra,dec),
                            EXPECTED_PAIRS[CID],float(frozen.loc[CID].baseline_days))
    ep,lp=get_pair_products(pair)
    early,ea=load_covering_cutouts(CID,ra,dec,"early","F444W",ep)
    late,la=load_covering_cutouts(CID,ra,dec,"late","F356W",lp)
    pd.DataFrame(ea+la).to_csv(out/"coverage.csv",index=False)
    orig,controls={},{}
    for ex in early+late:
        row,c=measure_exposure_recentered(ex)
        orig[ex.filename]=row
        controls[ex.filename]=c
    accepted=[ex for ex in early+late if target_is_usable(orig[ex.filename])]
    eaccepted=[ex for ex in accepted if ex.epoch_label=="early"]
    laccepted=[ex for ex in accepted if ex.epoch_label=="late"]
    if len(eaccepted)!=4 or len(laccepted)!=1:
        raise RuntimeError(f"Unexpected accepted original CALs: {len(eaccepted)} early {len(laccepted)} late")
    quality=[]
    for ex in accepted:
        q=orig[ex.filename]
        quality.append(dict(filename=ex.filename,filter=ex.filter_name,
                            epoch=ex.epoch_label,originally_accepted=True,
                            astrometry_source=q.get("astrometry_source"),
                            target_snr=q.get("astrometric_snr"),
                            dq_bad_core_3x3=inspect_core(ex,q),
                            n_controls=int(len(controls[ex.filename]))))
    pd.DataFrame(quality).to_csv(out/"original_accepted_target_core_quality.csv",index=False)
    bad_files=[ex.filename for ex in eaccepted if inspect_core(ex,orig[ex.filename])>0]
    scenarios=[("all_original_accepted",[])]
    scenarios += [("drop_"+e.filename.replace("_nrcblong_cal.fits",""),[e.filename]) for e in eaccepted]
    scenarios += [("drop_all_flagged_F444W_cores",bad_files)]
    results=[]
    for name,drop in scenarios:
        keep=[ex for ex in accepted if ex.filename not in drop]
        info=dict(scenario=name,excluded_f444w_frames=drop,
                  n_original_accepted_early_before=4,
                  n_early_after=len([ex for ex in keep if ex.epoch_label=="early"]))
        try:
            fit,reg,q,pairs=perform(keep,orig,controls)
            info.update(fit)
            reg.to_csv(out/(name+"_registered.csv"),index=False)
            q.to_csv(out/(name+"_registration_quality.csv"),index=False)
            pairs.to_csv(out/(name+"_pairs.csv"),index=False)
        except Exception as exc:
            info["status"]="FAILED:"+type(exc).__name__+":"+str(exc)
        results.append(info)
        print("EXPOSURE_LEAVEOUT",json.dumps(safe_out(info)),flush=True)
    pd.DataFrame([{**x,"excluded_f444w_frames":";".join(x["excluded_f444w_frames"])} for x in results]).to_csv(
        out/"original_estimator_f444w_leaveout_summary.csv",index=False)
    (out/"summary.json").write_text(json.dumps(safe_out(dict(
       candidate_id=CID,
       original_frozen_pm_masyr=[float(frozen.loc[CID].mu_alpha_cosdec_masyr),
                                float(frozen.loc[CID].mu_delta_masyr)],
       scenarios=results,
       science_status="SENSITIVITY_ONLY_NOT_NEW_PROPER_MOTION",
       warnings=[
        "Exact original first26 recentered target estimator, but may differ in frozen paper processing settings; do not overwrite paper.",
        "References change with omitted exposures; cross-filter systematic and parallax remain unknown.",
        "Removing poor-DQ-core F444W original accepted exposures can reveal sensitivity, but it is not a proof that kept exposures are unbiased.",
        "Only two observational epochs: parallax and intrinsic proper motion are not separable.",
        "Formal Gaussian 2D score/classification remains uncalibrated with current tiny held-out galaxy controls."
       ])),indent=2))
    print("SUMMARY_WRITTEN",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("results_1651050_f444w_frame_dropout"))
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    install()
    try:run(a.out)
    except Exception:
        (a.out/"failure.txt").write_text(traceback.format_exc())
        raise

if __name__=="__main__":main()
