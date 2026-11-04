"""Real JWST CAL cross-filter injection/recovery pilot.

This tests a representative original exposure pair, not the full epoch fit.
Requires public MAST exposures. PSFs are empirical compact field-source
stamps; they are not validated extreme-color ePSFs. Every trial re-runs the
existing target centroid/recenter, control detection, registration, and PM fit.
"""
import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import shift

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"first26"),str(ROOT/"scripts")]
from pm86.archive import query_candidate_inventory,get_pair_products,load_covering_cutouts
from pm86.cross_filter_injection import inject_at_sky
from pm86.injection_trial_summary import summarize_trials
from pm86.astrometry import register_exposures,infer_proper_motion
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install
from run_cross_filter_null import choose_frozen_pair,EXPECTED_PAIRS


def psf_from_controls(ex, controls):
    sci=np.asarray(ex.sci,float)
    dq=np.asarray(ex.dq)
    candidates=[]
    for _,row in controls.sort_values("snr",ascending=False).head(50).iterrows():
        x,y=float(row.x_local),float(row.y_local)
        xx,yy=int(round(x)),int(round(y))
        if min(xx,yy)<20 or xx+20>=sci.shape[1] or yy+20>=sci.shape[0]:
            continue
        p=sci[yy-7:yy+8,xx-7:xx+8]
        Y,X=np.indices(p.shape)
        r=np.hypot(X-7,Y-7)
        # JWST DQ has nonfatal flags. Reject defective PSF *cores*, but
        # exclude flagged wings from the empirical template instead of
        # requiring an impossible all-zero DQ square.
        mask=(dq[yy-7:yy+8,xx-7:xx+8] & 23)!=0
        mask=mask | ~np.isfinite(p)
        if np.any(mask[r<=2.5]) or np.count_nonzero(~mask)<160:
            continue
        sky=p[(r>5)&(~mask)]
        if len(sky)<25:
            continue
        p=np.maximum(np.where(mask,0,p-float(np.median(sky))),0)
        if p.sum()<=0:
            continue
        compact=float(p[r<2.5].sum()/p.sum())
        if compact<0.13 or p[r>5].sum()/p.sum()>0.6:
            continue
        p=np.maximum(shift(p,(yy-y,xx-x),order=3,mode="constant"),0)
        if p.sum()>0:
            candidates.append((compact,p/p.sum(),float(row.snr)))
    if not candidates:
        raise RuntimeError("No clean compact empirical PSF control in this CAL exposure")
    compact,p,snr=max(candidates,key=lambda k:k[0])
    return p,{"compactness":compact,"snr":snr,
              "warning":"unvetted field-source PSF, not an extreme-color ePSF"}


def clean_site(exposures,control_map,ra,dec):
    for e in exposures:
        try:
            xx,yy=e.gwcs.invert(ra,dec)
            x=float(xx)-e.x0
            y=float(yy)-e.y0
        except Exception:
            return False
        if not np.isfinite([x,y]).all():
            return False
        if min(x,y)<12 or x>e.sci.shape[1]-13 or y>e.sci.shape[0]-13:
            return False
        if np.hypot(x-(e.x_full-e.x0),y-(e.y_full-e.y0))<18:
            return False
        x,y=int(round(x)),int(round(y))
        cut=np.s_[y-7:y+8,x-7:x+8]
        if np.any(e.dq[cut]!=0) or not np.isfinite(e.sci[cut]).all() or not np.isfinite(e.err[cut]).all():
            return False
        c=control_map[e.filename]
        if len(c) and np.any(np.hypot(c.x_local-x,c.y_local-y)<16):
            return False
    return True


def measure_fake(exposures,psfs,flux,base,mu,moving,first_mjd):
    rows=[]
    controls=[]
    for e in exposures:
        dt=(e.mjd-first_mjd)/365.25
        pos=base.spherical_offsets_by(mu[0]*dt*u.mas,mu[1]*dt*u.mas) if moving else base
        try:
            fake,_=inject_at_sky(e,float(pos.ra.deg),float(pos.dec.deg),
                                 psfs[e.filename],flux[e.filename],strict=True)
            # Same fixed catalog sky coordinate at both epochs; the motion
            # must be found by centroiding, not provided to the fitter.
            x,y=e.gwcs.invert(float(base.ra.deg),float(base.dec.deg))
            fake=replace(fake,ra_deg=float(base.ra.deg),dec_deg=float(base.dec.deg),
                         x_full=float(x),y_full=float(y))
            row,c=measure_exposure_recentered(fake)
            rows.append(row)
            if len(c):
                controls.append(c)
        except Exception:
            return [np.nan,np.nan],"EXPOSURE_MEASUREMENT_FAILED"
    try:
        registered,q,matched=register_exposures(pd.DataFrame(rows),pd.concat(controls))
        if registered.empty or set(registered.epoch_label)!={"early","late"}:
            return [np.nan,np.nan],"NOT_BOTH_EPOCHS_REGISTERED"
        r,_=infer_proper_motion(registered,q,"CROSS_FILTER_JWST",matched)
        if r.get("pm_status")!="PM_MEASURED":
            return [np.nan,np.nan],r.get("pm_status","PM_FAILED")
        return [float(r["mu_alpha_cosdec_masyr"]),float(r["mu_delta_masyr"])],"RECOVERED"
    except Exception as err:
        return [np.nan,np.nan],f"PM_FIT_FAILED:{type(err).__name__}"


def run(cid,n,out):
    rng=np.random.default_rng(20261104+cid)
    catalog=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    original=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    target=catalog.loc[cid]
    frozen=original.loc[cid]
    pair=choose_frozen_pair(query_candidate_inventory(cid,float(target.ra_deg),float(target.dec_deg)),
                            EXPECTED_PAIRS[cid],float(frozen.baseline_days))
    early_p,late_p=get_pair_products(pair)
    early,a=load_covering_cutouts(cid,float(target.ra_deg),float(target.dec_deg),"early",EXPECTED_PAIRS[cid][0],early_p)
    late,b=load_covering_cutouts(cid,float(target.ra_deg),float(target.dec_deg),"late",EXPECTED_PAIRS[cid][1],late_p)
    pd.DataFrame(a+b).to_csv(out/"coverage.csv",index=False)
    if not early or not late:
        raise RuntimeError("No two-epoch original CAL images")
    # Representative original exposures, not an independent epoch-average.
    exposures=[early[0],late[0]]
    pmap,fmap,cmap,meta={}, {}, {}, {}
    for e in exposures:
        obs,controls=measure_exposure_recentered(e)
        p,info=psf_from_controls(e,controls)
        ap=float(obs["clean_flux"])
        if not np.isfinite(ap) or ap<=0:
            raise RuntimeError(f"Cannot preserve true faint-filter flux: {e.filename}")
        Y,X=np.indices(p.shape)
        fraction=float(p[np.hypot(X-7,Y-7)<=3].sum())
        if fraction<0.5:
            raise RuntimeError("Bad empirical PSF aperture correction")
        pmap[e.filename],fmap[e.filename],cmap[e.filename]=p,ap/fraction,controls
        meta[e.filename]={"filter":e.filter_name,"SNR":obs.get("astrometric_snr"),
                           "target_aperture_flux":ap,"injected_total_flux":ap/fraction,
                           "PSF":info}
    ref=early[0]
    centers=[]
    for _ in range(max(1000,120*n)):
        x=rng.uniform(13,ref.sci.shape[1]-14)
        y=rng.uniform(13,ref.sci.shape[0]-14)
        ra,dec=ref.gwcs(x+ref.x0,y+ref.y0)
        ra,dec=float(ra),float(dec)
        if not np.isfinite([ra,dec]).all():
            continue
        if not clean_site(exposures,cmap,ra,dec):
            continue
        sky=SkyCoord(ra*u.deg,dec*u.deg)
        if any(sky.separation(s).arcsec<0.7 for s in centers):
            continue
        centers.append(sky)
        if len(centers)>=n:
            break
    if len(centers)<5:
        raise RuntimeError("No adequate clean sites in original CAL overlap")
    true=np.array([float(frozen.mu_alpha_cosdec_masyr),float(frozen.mu_delta_masyr)])
    trials=[]
    for i,sky in enumerate(centers):
        for mode in ("stationary","moving"):
            pm,status=measure_fake(exposures,pmap,fmap,sky,true,mode=="moving",ref.mjd)
            trials.append({"trial":i,"mode":mode,"ra":float(sky.ra.deg),
                           "dec":float(sky.dec.deg),"mu_east":pm[0],
                           "mu_north":pm[1],"status":status})
        if (i+1)%10==0:
            print(f"{cid} actual CAL injection trial pairs: {i+1}/{len(centers)}",flush=True)
    df=pd.DataFrame(trials)
    df.to_csv(out/"trials.csv",index=False)
    null=df[df["mode"]=="stationary"][["mu_east","mu_north"]].to_numpy()
    moving=df[df["mode"]=="moving"][["mu_east","mu_north"]].to_numpy()
    result={"candidate_id":cid,"n_requested":n,"n_clean_sites":len(centers),
            "used_exposures":[e.filename for e in exposures],
            "filter_specific_flux_and_psf":meta,
            "status":"INJECTION_REAL_CAL_PILOT",
            "scientific_limitations":["Single original exposure pair, not exact multi-exposure paper fit",
                                    "Empirical PSF template not matched to ultracool target SED",
                                    "The blank-site injections do not independently calibrate common WCS errors",
                                    "500 trials cannot establish a Gaussian 5-sigma tail"]}
    if np.isfinite(null).all(axis=1).sum()>=2 and np.isfinite(moving).all(axis=1).sum()>=2:
        result["summary"]=summarize_trials(null,moving,true,true)
    else:
        result["status"]="INSUFFICIENT_RECOVERY"
        result["n_null_recovered"]=int(np.isfinite(null).all(axis=1).sum())
        result["n_moving_recovered"]=int(np.isfinite(moving).all(axis=1).sum())
    (out/"result.json").write_text(json.dumps(result,indent=2,default=float))
    print(json.dumps({k:v for k,v in result.items() if k!="filter_specific_flux_and_psf"},indent=2))
    if result["status"]!="INJECTION_REAL_CAL_PILOT":
        raise RuntimeError(result["status"])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--candidate",type=int,choices=[387132,1651050],required=True)
    p.add_argument("--trials",type=int,default=20)
    p.add_argument("--out",default="results_real_injection")
    a=p.parse_args()
    install()
    dest=Path(a.out)/str(a.candidate)
    dest.mkdir(parents=True,exist_ok=True)
    try:
        run(a.candidate,a.trials,dest)
    except Exception as exc:
        import traceback
        (dest/"fatal.txt").write_text(traceback.format_exc())
        print("REAL INJECTION FAILED:",repr(exc),flush=True)
        raise
