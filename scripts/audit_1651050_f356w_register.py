#!/usr/bin/env python3
"""1651050 F356W within-filter local-control registration sanity check.

Does NOT update the published proper motion. Reuses the target-specific shared
centroid diagnostic and tests whether the three rejected CAL frames still
constrain a common F356W location *after registering each detector's local
field sources to the already accepted F356W reference detector*.
All target fluxes are held out of control-source registration. Rejected target
detections are never promoted by numerical fit success alone.

Input: original public MAST NIRCam Stage-2 CAL + the pm86 control matcher.
"""
from __future__ import annotations
import json
import sys
import traceback
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"first26"),str(ROOT/"scripts")]
import pm86.measurement as meas
from pm86.astrometry import _match_controls, _fit_affine
from pm86.archive import query_candidate_inventory,get_pair_products,load_covering_cutouts
from run_cross_filter_null import choose_frozen_pair,EXPECTED_PAIRS,target_is_usable
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install
from audit_1651050_f356w_joint import make_frame,fit_shared,projected_matrix,safe_out,WIDTHS


def registered_projected_frame(ex,cat,ref_sky,coeff):
    # Original pm86 registration convention: [1,raw_east,raw_north] @ coeff
    # yields reference frame tangent-plane (east,north) in arcseconds.
    B=np.asarray(coeff,dtype=float)[1:,:]
    if not np.isfinite(B).all() or abs(np.linalg.det(B))<0.2:
        raise ValueError("singular or implausible affine WCS registration")
    inv=np.linalg.inv(B)
    origin=SkyCoord(float(cat[0])*u.deg,float(cat[1])*u.deg)
    en=meas.sky_to_tangent(float(ref_sky[0]),float(ref_sky[1]),
                           float(cat[0]),float(cat[1]))
    target_ref=np.array([float(en[0]),float(en[1])])
    target_raw=(target_ref-np.asarray(coeff[0,:]))@inv
    sky=origin.spherical_offsets_by(target_raw[0]*u.arcsec,target_raw[1]*u.arcsec)
    center,de,dn=projected_matrix(ex,float(sky.ra.deg),float(sky.dec.deg))
    Jraw=np.column_stack([de,dn])
    Jref=Jraw@inv.T
    frame=make_frame(ex.sci,ex.err,ex.dq,*center,Jref[0,0],Jref[0,1],
                     Jref[1,0],Jref[1,1],filename=ex.filename)
    return frame,target_ref,target_raw


def run(out):
    catdf=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    olddf=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    cid=1651050
    ra,dec=map(float,[catdf.loc[cid].ra_deg,catdf.loc[cid].dec_deg])
    pair=choose_frozen_pair(query_candidate_inventory(cid,ra,dec),
                            EXPECTED_PAIRS[cid],float(olddf.loc[cid].baseline_days))
    _,lp=get_pair_products(pair)
    exposures,coverage=load_covering_cutouts(cid,ra,dec,"late","F356W",lp)
    pd.DataFrame(coverage).to_csv(out/"mast_coverage.csv",index=False)
    records,controls={},{}
    for e in exposures:
        r,c=measure_exposure_recentered(e)
        records[e.filename]=r
        controls[e.filename]=c
    accepted=[x for x in exposures if target_is_usable(records[x.filename])]
    if len(accepted)!=1:
        raise RuntimeError(f"Expected one original accepted late F356W, got {len(accepted)}")
    ref=accepted[0]
    r=records[ref.filename]
    ref_sky=meas.local_pixel_to_sky(ref,r["x_2dg"],r["y_2dg"])
    cref=controls[ref.filename]
    if len(cref)<6:raise RuntimeError("Reference exposure has fewer than 6 field registration controls")
    frames_raw=[]
    frames_reg=[]
    diagnostics=[]
    for ex in exposures:
        try:
            p,de,dn=projected_matrix(ex,*ref_sky)
            frames_raw.append(make_frame(ex.sci,ex.err,ex.dq,*p,
                         de[0],dn[0],de[1],dn[1],filename=ex.filename))
            if ex.filename==ref.filename:
                coeff=np.array([[0.,0.],[1.,0.],[0.,1.]])
                matched=len(cref);good=np.ones(matched,bool)
                scatter=[0.,0.]
            else:
                current=controls[ex.filename]
                cidx,ridx=_match_controls(current,cref)
                p=current.iloc[cidx][["east_arcsec_raw_wcs",
                                      "north_arcsec_raw_wcs"]].to_numpy(float)
                q=cref.iloc[ridx][["east_arcsec_raw_wcs",
                                   "north_arcsec_raw_wcs"]].to_numpy(float)
                coeff,good,residual,sx,sy=_fit_affine(p,q)
                matched=len(cidx);scatter=[float(sx)*1000,float(sy)*1000]
            frame,target_ref,target_raw=registered_projected_frame(
                ex,(ra,dec),ref_sky,coeff)
            frames_reg.append(frame)
            diagnostics.append(dict(filename=ex.filename,
                                    originally_accepted=bool(target_is_usable(records[ex.filename])),
                                    n_field_controls=len(controls[ex.filename]),
                                    n_reference_controls=len(cref),
                                    n_matched=int(matched),n_used=int(np.sum(good)),
                                    registration_scatter_east_mas=float(scatter[0]),
                                    registration_scatter_north_mas=float(scatter[1]),
                                    affine_coeff=np.asarray(coeff).tolist(),
                                    inverse_projected_sky_east_arcsec=float(target_raw[0]),
                                    inverse_projected_sky_north_arcsec=float(target_raw[1]),
                                    source_position_ref_east_arcsec=float(target_ref[0]),
                                    source_position_ref_north_arcsec=float(target_ref[1])))
        except Exception as exc:
            diagnostics.append(dict(filename=ex.filename,
                                    originally_accepted=bool(target_is_usable(records[ex.filename])),
                                    status="REGISTRATION_UNAVAILABLE",
                                    error=type(exc).__name__+":"+str(exc)))
    (out/"per_frame_registration_quality.json").write_text(
        json.dumps(safe_out(diagnostics),indent=2))
    if len(frames_reg)<4 or len(frames_raw)<4:
        raise RuntimeError("Local registration did not cover all four F356W CAL images; check file")

    results=[]
    for model,frames in [("raw_wcs",frames_raw),("same_filter_field_registered",frames_reg)]:
        for w in WIDTHS:
            for leave in [None,ref.filename]:
                ff=[f for f in frames if f["filename"]!=leave]
                res=fit_shared(ff,float(w))
                res.update(model=model,left_out_frame=leave,
                           original_accepted_frame=ref.filename,
                           source_original_quality_unchanged=True)
                results.append(res)
                print("REGISTERED_DIAGNOSTIC",json.dumps(safe_out({
                    "model":model,"width":w,"leave":leave,
                    "offset_east_mas":res["offset_east_mas"],
                    "offset_north_mas":res["offset_north_mas"],
                    "total_chi2":res["total_chi2"]})),flush=True)
    (out/"registered_vs_raw_complete.json").write_text(json.dumps(safe_out(results),indent=2))
    pd.DataFrame([{k:v for k,v in res.items() if k!="per_frame"} for res in results]).to_csv(
        out/"registered_vs_raw_summary.csv",index=False)
    rows=[]
    for x in results:
        for p in x["per_frame"]:
            rows.append({**p,"model":x["model"],"width_pixels":x["width_pixels"],
                         "left_out_frame":x["left_out_frame"]})
    pd.DataFrame(rows).to_csv(out/"registered_vs_raw_per_frame_amplitudes.csv",index=False)

    summary={"candidate_id":cid,"n_f356w_covering":len(exposures),
             "n_original_accepted":len(accepted),
             "n_same_filter_registered":len(frames_reg),
             "accepted_original_frame":ref.filename,
             "registrations":diagnostics,
             "fitted_centroids":[{k:v for k,v in x.items() if k!="per_frame"} for x in results],
             "science_status":"UNVALIDATED_FIXED_GAUSSIAN_SHARED_LOCATION_DIAGNOSTIC_ONLY",
             "limitations":[
                 "This refines same-band local field registration but not ultra-red F356W PSF validity.",
                 "Raw and registered centroid offsets in this audit are relative to the originally accepted F356W sky position.",
                 "No independent trigonometric parallax, absolute astrometric uncertainties or F444W/F356W color-PMF calibration.",
                 "SNR values conditional on fixed/refitted position cannot be interpreted as independent detections.",
                 "Model still uses analytic Gaussian and constant sky background; neighboring source wings may bias the result.",
                 "Original rejection rules remain in effect, including DO_NOT_USE pixel flags.",
                 "No injected-source recovery or article PM is changed by this script."
             ]}
    (out/"summary.json").write_text(json.dumps(safe_out(summary),indent=2))
    print("AUDIT_SUMMARY",json.dumps(safe_out(summary)),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,default=Path("results_1651050_f356w_registered"))
    args=parser.parse_args()
    install()
    args.out.mkdir(parents=True,exist_ok=True)
    try:run(args.out)
    except Exception:
        (args.out/"failure.txt").write_text(traceback.format_exc())
        raise

if __name__=="__main__":main()
