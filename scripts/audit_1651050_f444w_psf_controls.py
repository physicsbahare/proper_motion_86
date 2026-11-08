#!/usr/bin/env python3
"""Bounded real-CAL 1651050 early-F444W / F356W PSF-versus-Gaussian audit.

Model-check independent local compact field-source controls in each detector,
then compare target detector centroids using STPSF and Gaussian on identical
unmodified SCI/ERR/DQ pixels. Only four *originally accepted* early F444W CALs
and the one accepted late F356W CAL. No promotion of rejected frames.

NOTE: field-source compactness is NOT spectroscopy-based verification that a
control is a star; STPSF oversampled interpolation is NOT an externally
calibrated mas-level PSF. All motion changes reported below are conditional,
raw-WCS method sensitivity, NOT new registered/PSF-corrected proper motion.
"""
from __future__ import annotations
import argparse
import json
import sys
import traceback
from pathlib import Path
from scipy.optimize import least_squares
import numpy as np
import pandas as pd
import astropy.units as u
from synphot import SourceSpectrum
from synphot.models import BlackBodyNorm1D

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"first26"),str(ROOT/"scripts")]
import pm86.measurement as meas
from pm86.archive import query_candidate_inventory,get_pair_products,load_covering_cutouts
from pm86.astrometry import _match_controls,_fit_affine
from run_cross_filter_null import choose_frozen_pair,EXPECTED_PAIRS,target_is_usable
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install
from audit_1651050_f356w_joint import safe_out,make_frame,projected_matrix
from audit_1651050_stpsf_neighbor import native_template

CID=1651050
FIT_HALF=6
MAX_CONTROLS_PER_CAL=3
MAX_TARGET_SHIFT_PIX=1.4
MAX_CONTROL_SHIFT_PIX=1.0
MAX_CAL_PER_FILTER=4
CONTROL_SNR_MIN=12.0
TEST_TEMPS=(500,900)              # hypothetical target source; NOT measured SED
CONTROL_TEMP=5700                 # hypothetical ordinary-star spectrum


def synthetic_self_test():
    # Check identifiability / fit residual output on a known smooth image.
    yy,xx=np.indices((41,41))
    center=(20.13,20.31)
    sci=3+55*np.exp(-.5*((xx-center[0])**2+(yy-center[1])**2)/1.1**2)
    err=np.ones_like(sci)*.6
    dq=np.zeros_like(sci,dtype=np.uint32)
    ex=type("E",(),dict(sci=sci,err=err,dq=dq))()
    f=stamp(ex,20,20,half=FIT_HALF)
    p=fit_gaussian(f,max_shift=.9)
    assert p["success"] and p["source_amplitude"]>0 and np.hypot(
        p["x_pix"]-center[0],p["y_pix"]-center[1])<0.15,p
    print("GAUSSIAN_STAMP_SELF_TEST_SUCCESS",json.dumps(safe_out(p)))


def stamp(ex,x0,y0,half=FIT_HALF):
    image=np.asarray(ex.sci,float)
    err=np.asarray(ex.err,float)
    dq=np.asarray(ex.dq,np.uint32)
    _,bad=meas.masks(image,err,dq)
    bad |= ~(np.isfinite(image)&np.isfinite(err)&(err>0))
    xc,yc=int(round(x0)),int(round(y0))
    cut=np.s_[yc-half:yc+half+1,xc-half:xc+half+1]
    if (xc<half or yc<half or xc+half>=image.shape[1] or
        yc+half>=image.shape[0]):
        raise ValueError("stamp outside detector local cutout")
    ok=~bad[cut]
    if ok.sum()<45:raise ValueError("too few usable pixels")
    yy,xx=np.mgrid[yc-half:yc+half+1,xc-half:xc+half+1]
    return dict(x0=float(x0),y0=float(y0),xx=xx[ok].astype(float),
                yy=yy[ok].astype(float),sci=image[cut][ok],
                err=err[cut][ok],n_valid=int(ok.sum()),
                n_bad=int(ok.size-ok.sum()),core_bad=int(
                    bad[max(0,yc-1):yc+2,max(0,xc-1):xc+2].sum()))


def fit_2d_source(f,psf=None,max_shift=MAX_TARGET_SHIFT_PIX):
    """Profile signed amplitude + planar background; optimize 2D source center."""
    def profile(x,y):
        if psf is None:
            model=np.exp(-.5*((f["xx"]-x)**2+(f["yy"]-y)**2)/1.05**2)
        else:
            model=native_template(psf,f["xx"],f["yy"],x,y)
        X=np.column_stack([model,np.ones(len(model)),
                           f["xx"]-f["x0"],f["yy"]-f["y0"]])
        w=1/f["err"]
        coeff, *_=np.linalg.lstsq(X*w[:,None],f["sci"]*w,rcond=None)
        r=(f["sci"]-X@coeff)*w
        return r,coeff
    def res(v):return profile(float(v[0]),float(v[1]))[0]
    x0,y0=f["x0"],f["y0"]
    result=least_squares(res,[x0,y0],
                         bounds=([x0-max_shift,y0-max_shift],
                                 [x0+max_shift,y0+max_shift]),
                         max_nfev=90,ftol=1e-7,xtol=1e-7,gtol=1e-7)
    r,c=profile(*result.x)
    return dict(success=bool(result.success),x_pix=float(result.x[0]),
                y_pix=float(result.x[1]),source_amplitude=float(c[0]),
                chi2=float(np.dot(r,r)),chi2_per_valid_pixel=float(np.mean(r**2)),
                n_valid=f["n_valid"],n_bad=f["n_bad"],core_bad=f["core_bad"],
                on_position_bound=bool(max(abs(result.x-[x0,y0]))>.95*max_shift))


def fit_gaussian(f,max_shift=MAX_TARGET_SHIFT_PIX):
    return fit_2d_source(f,psf=None,max_shift=max_shift)


def physical_psf(ex,x_local,y_local,temperature):
    import stpsf
    inst=stpsf.NIRCam()
    inst.filter=str(ex.filter_name)
    inst.detector=str(ex.detector).upper()
    inst.detector_position=(int(round(ex.x0+x_local)),int(round(ex.y0+y_local)))
    source=SourceSpectrum(BlackBodyNorm1D,temperature=temperature*u.K)
    hdul=inst.calc_psf(source=source,nlambda=7,fov_pixels=33,oversample=4)
    if "OVERDIST" not in hdul:
        raise RuntimeError("Missing STPSF detector-distorted oversampled PSF")
    image=np.asarray(hdul["OVERDIST"].data,float)
    if image.shape!=(132,132) or not np.isfinite(image).all() or image.sum()<=0:
        raise ValueError("Incorrect physical PSF sampling")
    image=np.maximum(image,0)
    image/=image.sum()
    hdul.close()
    return image,dict(stpsf_version=stpsf.__version__,filter=inst.filter,
                      detector=inst.detector,detector_position=list(inst.detector_position),
                      temperature_hypothesis_K=int(temperature))


def safe_fit(ex,x,y,temp,seed_type):
    f=stamp(ex,x,y)
    gaussian=fit_gaussian(f,max_shift=(MAX_TARGET_SHIFT_PIX if seed_type=="target" else MAX_CONTROL_SHIFT_PIX))
    img,meta=physical_psf(ex,x,y,temp)
    phys=fit_2d_source(f,psf=img,max_shift=(MAX_TARGET_SHIFT_PIX if seed_type=="target" else MAX_CONTROL_SHIFT_PIX))
    if np.isfinite([gaussian["x_pix"],gaussian["y_pix"],phys["x_pix"],phys["y_pix"]]).all():
        ra0,dec0=meas.local_pixel_to_sky(ex,gaussian["x_pix"],gaussian["y_pix"])
        ra1,dec1=meas.local_pixel_to_sky(ex,phys["x_pix"],phys["y_pix"])
        e,n=meas.sky_to_tangent(ra1,dec1,ra0,dec0)
        de,dn=float(e)*1000,float(n)*1000
    else:de=dn=np.nan
    return dict(gaussian=gaussian,stpsf=phys,
                delta_stpsf_minus_gaussian_east_mas=de,
                delta_stpsf_minus_gaussian_north_mas=dn,
                meta=meta)


def select_controls(ex,controls,tx,ty,nmax=MAX_CONTROLS_PER_CAL):
    """Strong isolated, unflagged compact *candidate* points, not certified stars."""
    d=controls.sort_values("snr",ascending=False)
    kept=[]
    rejected=[]
    xy=d[["x_local","y_local"]].to_numpy(float)
    sci=np.asarray(ex.sci,float)
    err=np.asarray(ex.err,float)
    _,bad=meas.masks(sci,err,np.asarray(ex.dq,np.uint32))
    ny,nx=sci.shape
    for _,row in d.iterrows():
        x,y=float(row.x_local),float(row.y_local)
        snr=float(row.snr)
        reason=None
        if snr<CONTROL_SNR_MIN:break
        if min(x,y)<FIT_HALF+2 or x>=nx-FIT_HALF-2 or y>=ny-FIT_HALF-2:
            reason="EDGE"
        elif np.hypot(x-tx,y-ty)<16:
            reason="TARGET_NEAR"
        elif any(np.hypot(x-q[0],y-q[1])<11 for q in kept):
            reason="NEAR_OTHER_SELECTED_CONTROL"
        else:
            nearest=[np.hypot(x-a,y-b) for a,b in xy if np.hypot(x-a,y-b)>.5]
            if len(nearest) and min(nearest)<10:
                reason="CROWDED_FIELD"
            else:
                yc,xc=int(round(y)),int(round(x))
                if np.any(bad[yc-2:yc+3,xc-2:xc+3]):
                    reason="BAD_CORE_DQ"
                elif not np.isfinite(err[yc-2:yc+3,xc-2:xc+3]).all():
                    reason="BAD_CORE_ERR"
        if reason:
            rejected.append(dict(x=x,y=y,snr=snr,reason=reason))
            continue
        kept.append((x,y,snr))
        if len(kept)>=nmax:break
    return kept,rejected


def run(out):
    cat=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    frozen=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    ra,dec=float(cat.loc[CID].ra_deg),float(cat.loc[CID].dec_deg)
    pair=choose_frozen_pair(query_candidate_inventory(CID,ra,dec),
                            EXPECTED_PAIRS[CID],float(frozen.loc[CID].baseline_days))
    ep,lp=get_pair_products(pair)
    early,ea=load_covering_cutouts(CID,ra,dec,"early","F444W",ep)
    late,la=load_covering_cutouts(CID,ra,dec,"late","F356W",lp)
    pd.DataFrame(ea+la).to_csv(out/"cal_coverage.csv",index=False)
    all_ex=early+late
    quality={}
    controls={}
    exmap={}
    for ex in all_ex:
        original,cs=measure_exposure_recentered(ex)
        quality[ex.filename]=original
        controls[ex.filename]=cs
        exmap[ex.filename]=ex
    usable=[ex for ex in all_ex if target_is_usable(quality[ex.filename])]
    if len(usable)!=5 or sum(ex.epoch_label=="late" for ex in usable)!=1:
        raise RuntimeError("Need exactly 4 originally target-accepted F444W + 1 F356W; got "+
                           str([(ex.filename,ex.epoch_label) for ex in usable]))
    audit=[]
    control_rows=[]
    target_rows=[]
    psf_meta=[]
    rejection_rows=[]
    for ex in usable:
        old=quality[ex.filename]
        tx,ty=float(old["x_2dg"]),float(old["y_2dg"])
        if not np.isfinite([tx,ty]).all():raise RuntimeError("Invalid original source centroid")
        chosen,rejected=select_controls(ex,controls[ex.filename],tx,ty)
        for item in rejected:
            rejection_rows.append(dict(filename=ex.filename,filter=ex.filter_name,**item))
        audit.append(dict(filename=ex.filename,filter=ex.filter_name,
                          original_astrometry_source=old.get("astrometry_source"),
                          original_astrometric_snr=old.get("astrometric_snr"),
                          originally_accepted=True,
                          n_detected_field_controls=int(len(controls[ex.filename])),
                          n_isolated_control_candidates=len(chosen)))
        for temp in TEST_TEMPS:
            try:
                t=safe_fit(ex,tx,ty,temp,"target")
                target_rows.append(dict(filename=ex.filename,filter=ex.filter_name,
                    epoch=ex.epoch_label,target_temperature_assumption_K=temp,
                    orig_x=tx,orig_y=ty,
                    gauss_x=t["gaussian"]["x_pix"],gauss_y=t["gaussian"]["y_pix"],
                    stpsf_x=t["stpsf"]["x_pix"],stpsf_y=t["stpsf"]["y_pix"],
                    gaussian_amp=t["gaussian"]["source_amplitude"],
                    stpsf_amp=t["stpsf"]["source_amplitude"],
                    gaussian_chi2=t["gaussian"]["chi2"],
                    stpsf_chi2=t["stpsf"]["chi2"],
                    delta_east_mas=t["delta_stpsf_minus_gaussian_east_mas"],
                    delta_north_mas=t["delta_stpsf_minus_gaussian_north_mas"],
                    stpsf_success=t["stpsf"]["success"],
                    stpsf_position_bound=t["stpsf"]["on_position_bound"],
                    stpsf_core_bad=t["stpsf"]["core_bad"]))
                psf_meta.append(dict(filename=ex.filename,source_type="target",**t["meta"]))
            except Exception as exc:
                target_rows.append(dict(filename=ex.filename,filter=ex.filter_name,
                                        epoch=ex.epoch_label,target_temperature_assumption_K=temp,
                                        status="FAILED:"+type(exc).__name__+":"+str(exc)))
            print("TARGET_PSF",ex.filter_name,ex.filename,temp,flush=True)
        for ix,(cx,cy,snr) in enumerate(chosen):
            try:
                t=safe_fit(ex,cx,cy,CONTROL_TEMP,"control")
                control_rows.append(dict(filename=ex.filename,filter=ex.filter_name,
                                         field_control_index=ix,
                                         candidate_compact_control_snr=snr,
                                         x_detect=cx,y_detect=cy,
                                         gaussian_x=t["gaussian"]["x_pix"],
                                         gaussian_y=t["gaussian"]["y_pix"],
                                         stpsf_x=t["stpsf"]["x_pix"],
                                         stpsf_y=t["stpsf"]["y_pix"],
                                         gaussian_chi2=t["gaussian"]["chi2"],
                                         stpsf_chi2=t["stpsf"]["chi2"],
                                         delta_east_mas=t["delta_stpsf_minus_gaussian_east_mas"],
                                         delta_north_mas=t["delta_stpsf_minus_gaussian_north_mas"],
                                         gaussian_amp=t["gaussian"]["source_amplitude"],
                                         stpsf_amp=t["stpsf"]["source_amplitude"],
                                         stpsf_success=t["stpsf"]["success"],
                                         stpsf_position_bound=t["stpsf"]["on_position_bound"],
                                         stpsf_core_bad=t["stpsf"]["core_bad"]))
                psf_meta.append(dict(filename=ex.filename,source_type="control",**t["meta"]))
            except Exception as exc:
                control_rows.append(dict(filename=ex.filename,filter=ex.filter_name,
                                         field_control_index=ix,
                                         status="FAILED:"+type(exc).__name__+":"+str(exc)))
            print("FIELD_CONTROL_PSF",ex.filter_name,ex.filename,ix,flush=True)
        # Save output even if one physical PSF fails in a subsequent exposure.
        pd.DataFrame(target_rows).to_csv(out/"target_f444w_f356w_stpsf_vs_gauss.csv",index=False)
        pd.DataFrame(control_rows).to_csv(out/"independent_field_compact_psf_controls.csv",index=False)
        (out/"physical_psf_meta.json").write_text(json.dumps(safe_out(psf_meta),indent=2))
    pd.DataFrame(audit).to_csv(out/"exposure_original_quality_and_control_counts.csv",index=False)
    pd.DataFrame(rejection_rows).to_csv(out/"rejected_field_control_candidates.csv",index=False)

    t=pd.DataFrame(target_rows)
    c=pd.DataFrame(control_rows)
    q=(t[t.get("stpsf_success",pd.Series(dtype=bool)).eq(True) &
         (pd.to_numeric(t.get("stpsf_amp",pd.Series(dtype=float)),errors="coerce")>0)].copy()
       if len(t) else pd.DataFrame())
    cc=(c[c.get("stpsf_success",pd.Series(dtype=bool)).eq(True) &
          (pd.to_numeric(c.get("stpsf_amp",pd.Series(dtype=float)),errors="coerce")>0)].copy()
        if len(c) else pd.DataFrame())
    per_filter={}
    for filt in ["F444W","F356W"]:
        ctrl=cc[cc["filter"].eq(filt)] if len(cc) else pd.DataFrame()
        tar=q[q["filter"].eq(filt)] if len(q) else pd.DataFrame()
        per_filter[filt]=dict(n_original_accepted_exposures=int(sum(ex.filter_name==filt for ex in usable)),
                              n_control_candidates=int(sum(len(select_controls(ex,controls[ex.filename],
                                  float(quality[ex.filename]["x_2dg"]),float(quality[ex.filename]["y_2dg"]))[0])
                                  for ex in usable if ex.filter_name==filt)),
                              n_usable_psf_control_fits=int(len(ctrl)),
                              control_model_shift_east_median_mas=float(ctrl["delta_east_mas"].median()) if len(ctrl) else None,
                              control_model_shift_north_median_mas=float(ctrl["delta_north_mas"].median()) if len(ctrl) else None,
                              control_model_shift_2d_median_mas=float(np.median(np.hypot(ctrl["delta_east_mas"],ctrl["delta_north_mas"]))) if len(ctrl) else None,
                              n_target_psf_model_fits=int(len(tar)))
    # Only a descriptive raw-WCS estimator sensitivity. NO new F444W/F356W
    # calibrated PM or uncertainty; the nominal single-late difference is known.
    diffs={}
    for temp in TEST_TEMPS:
        a=q[(q["filter"].eq("F444W")) & (q["target_temperature_assumption_K"].eq(temp))]
        b=q[(q["filter"].eq("F356W")) & (q["target_temperature_assumption_K"].eq(temp))]
        if len(a)==4 and len(b)==1:
            de=float(b["delta_east_mas"].iloc[0]-a["delta_east_mas"].median())
            dn=float(b["delta_north_mas"].iloc[0]-a["delta_north_mas"].median())
            diffs[str(temp)]=dict(delta_late_minus_median_early_east_mas=de,
                                  delta_late_minus_median_early_north_mas=dn,
                                  conditional_annualized_difference_east_masyr=de/(float(frozen.loc[CID].baseline_days)/365.25),
                                  conditional_annualized_difference_north_masyr=dn/(float(frozen.loc[CID].baseline_days)/365.25),
                                  NOTE="raw WCS PSF-vs-Gaussian centroid sensitivity; NOT a re-registered motion, correction, or Gaussian sigma")
    summary=dict(candidate_id=CID,status="METHOD_SENSITIVITY_AND_UNVETTED_COMPACT_CONTROL_ONLY",
                 original_accepted_exposures=5,
                 per_filter=per_filter,
                 shifts_between_two_filters=diffs,
                 tested_hypothetical_target_blackbodies_K=list(TEST_TEMPS),
                 fixed_hypothetical_control_star_blackbody_K=CONTROL_TEMP,
                 caution=[
                  "Selected controls are isolated bright compact detections, NOT spectroscopically verified stars or proven stationary; per-exposure few-source sample.",
                  "Same CALs and noisy pixels for Gaussian vs STPSF comparison, no ground truth; a model closer in chi2 need not be astrometrically unbiased.",
                  "STPSF oversampled interpolation and pixel phase were not independently validated at mas precision; color SEDs are blackbody hypotheses, not ultracool spectra.",
                  "No spectral or spatial PSF calibration from early vs late filters; do not apply raw method differences as PM corrections.",
                  "Current STPSF target fitting uses local stamps with planar background; neighbor wings may enter target image and controls in unmodeled ways.",
                  "This is not a validated new proper motion or new sigma, cannot distinguish JWST parallax with two epochs.",
                  "No rejected original target exposure or journal PM value is changed."
                 ])
    (out/"summary.json").write_text(json.dumps(safe_out(summary),indent=2))
    print("F444W_F356W_CONTROL_AUDIT_SUMMARY",json.dumps(safe_out(summary)),flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("results_1651050_f444w_controls"))
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        synthetic_self_test();return
    a.out.mkdir(parents=True,exist_ok=True)
    install()
    try:run(a.out)
    except Exception:
        (a.out/"failure.txt").write_text(traceback.format_exc())
        raise


if __name__=="__main__":
    main()
