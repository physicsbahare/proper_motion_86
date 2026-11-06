"""Full-epoch two-filter astrometric control experiment for Dave PM candidates.

Use exactly the original candidate-centroidable CAL exposures. Re-register each
exposure to an early reference using vetted COSMOS2025 TYPE=0 candidates. For
each source, EXCLUDE IT from every fit before estimating its full-epoch PM null.
These are photometric galaxies, not independent motion-confirmed galaxies.
No adopted PM significance or cross-filter PSF systematic is computed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from .astrometry import _fit_affine, _position_errors, _weighted_epoch, infer_proper_motion

XY=("east_arcsec_raw_wcs","north_arcsec_raw_wcs")
IDENTITY=np.array([[0.,0.],[1.,0.],[0.,1.]])

def match_vetted(controls: pd.DataFrame, vetted: pd.DataFrame, radius_arcsec=.25):
    """Unique nearest-on-sky matches; controls must be independently detected."""
    cols=list(controls.columns)+["galaxy_id","catalog_separation_arcsec"]
    if controls.empty:
        return pd.DataFrame(columns=cols)
    required={"ra_deg","dec_deg","snr",*XY}
    if not required.issubset(controls.columns):
        raise ValueError("Original per-exposure control columns missing")
    a=controls.reset_index(drop=True).copy()
    for key in ["ra_deg","dec_deg","snr"]:
        a[key]=pd.to_numeric(a[key],errors="coerce")
    a=a[np.isfinite(a[["ra_deg","dec_deg","snr"]]).all(axis=1)&(a.snr>=10)].copy()
    if a.empty:return pd.DataFrame(columns=cols)
    ca=SkyCoord(a.ra_deg.values*u.deg,a.dec_deg.values*u.deg)
    cb=SkyCoord(vetted.ra_deg.to_numpy(float)*u.deg,
                vetted.dec_deg.to_numpy(float)*u.deg)
    idx,dist,_=ca.match_to_catalog_sky(cb)
    a["galaxy_id"]=vetted.iloc[idx].catalog_id.to_numpy(int)
    a["catalog_separation_arcsec"]=dist.arcsec
    a=a[a.catalog_separation_arcsec<=radius_arcsec].copy()
    return a.sort_values("catalog_separation_arcsec").drop_duplicates("galaxy_id").reset_index(drop=True)


def choose_anchor(measurements:pd.DataFrame, control_map:dict):
    early=measurements[measurements.epoch_label.eq("early")].copy()
    if early.empty:raise ValueError("No usable early epoch")
    early["_n"]=early.filename.map(lambda f:len(control_map.get(str(f),[])))
    return str(early.sort_values(["_n","mjd"],ascending=[False,True]).iloc[0].filename)


def calculate_affine(ref_controls,other_controls,exclude_id=None,min_fit=6):
    """Full 6-parameter affine fit using only shared, independently vetted galaxies."""
    ref=ref_controls.set_index("galaxy_id",drop=False)
    cur=other_controls.set_index("galaxy_id",drop=False)
    common=sorted(set(ref.index).intersection(cur.index)-({int(exclude_id)} if exclude_id is not None else set()))
    if len(common)<min_fit:
        return None,{"status":"INSUFFICIENT_COMMON_GALAXIES","n_common":len(common),"n_fit":0}
    src=cur.loc[common,list(XY)].to_numpy(float)
    dst=ref.loc[common,list(XY)].to_numpy(float)
    if np.linalg.matrix_rank(np.c_[np.ones(len(src)),src])<3:
        return None,{"status":"DEGENERATE_GEOMETRY","n_common":len(common),"n_fit":0}
    try:
        coef,good,resid,sx,sy=_fit_affine(src,dst,min_matches=min_fit)
    except (RuntimeError,np.linalg.LinAlgError,ValueError) as exc:
        return None,{"status":"AFFINE_FAILED","n_common":len(common),"reason":str(exc)}
    return coef,{"status":"OK","n_common":len(common),"n_fit":int(good.sum()),
                 "scatter_east_mas":float(sx*1000),"scatter_north_mas":float(sy*1000),
                 "galaxy_ids_fit":[common[i] for i in np.flatnonzero(good)]}


def transform(coef,east,north):
    return np.r_[1.,east,north]@coef


def combined_epoch_measurement(rows:pd.DataFrame,registered_target:pd.DataFrame):
    """Use exactly target inverse-variance epoch weights; not galaxy errors."""
    if rows.empty or registered_target.empty:return {"status":"NO_POSITIONS"}
    frames=registered_target.reset_index(drop=True)
    e,n=_position_errors(frames)
    wt=frames[["filename","mjd","epoch_label"]].copy()
    wt["target_weight_err_e"]=e
    wt["target_weight_err_n"]=n
    r=rows.merge(wt,on=["filename","mjd","epoch_label"],how="inner")
    early=r[r.epoch_label.eq("early")]
    late=r[r.epoch_label.eq("late")]
    if early.empty or late.empty:
        return {"status":"MISSING_EPOCH","n_early":len(early),"n_late":len(late)}
    ee,_=_weighted_epoch(early.east_mas,early.target_weight_err_e)
    en,_=_weighted_epoch(early.north_mas,early.target_weight_err_n)
    le,_=_weighted_epoch(late.east_mas,late.target_weight_err_e)
    ln,_=_weighted_epoch(late.north_mas,late.target_weight_err_n)
    dt=float((late.mjd.mean()-early.mjd.mean())/365.25)
    if not np.isfinite([ee,en,le,ln,dt]).all() or dt<=0:
        return {"status":"INVALID_COMBINATION"}
    de,dn=float(le-ee),float(ln-en)
    return {"status":"FULL_EPOCH_GALAXY_NULL_MEASURED",
            "n_early":int(len(early)),"n_late":int(len(late)),
            "baseline_days":dt*365.25,"delta_east_mas":de,"delta_north_mas":dn,
            "mu_east_masyr":de/dt,"mu_north_masyr":dn/dt,
            "weights":"same_target_error_weight_by_exposure_NOT_galaxy_centroid_uncertainties"}


def registered_solution(measurements, control_map, anchor, heldout=None, min_fit=6):
    """Per-exposure transform followed by original full-epoch target PM estimator."""
    anchor_controls=control_map[anchor]
    if (heldout is not None and int(heldout) not in anchor_controls.galaxy_id.values):
        return pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),{"status":"HELDOUT_NOT_IN_ANCHOR"}
    reg=[];quality=[];src_rows=[];skipped=[]
    for _,row in measurements.sort_values("mjd").iterrows():
        name=str(row.filename)
        controls=control_map.get(name,pd.DataFrame())
        if name==anchor:
            coef=IDENTITY
            fit={"status":"OK","n_common":len(anchor_controls),"n_fit":len(anchor_controls),
                 "scatter_east_mas":0.,"scatter_north_mas":0.}
        else:
            coef,fit=calculate_affine(anchor_controls,controls,heldout,min_fit)
        if coef is None:
            skipped.append({"filename":name,"epoch":str(row.epoch_label),**fit})
            continue
        transformed=transform(coef,float(row.east_2dg_arcsec_raw_wcs),
                                   float(row.north_2dg_arcsec_raw_wcs))
        out=row.to_dict()
        out["east_2dg_mas"]=float(transformed[0]*1000)
        out["north_2dg_mas"]=float(transformed[1]*1000)
        if np.isfinite(float(row.east_com_arcsec_raw_wcs)) and np.isfinite(float(row.north_com_arcsec_raw_wcs)):
            co=transform(coef,float(row.east_com_arcsec_raw_wcs),float(row.north_com_arcsec_raw_wcs))
            out["east_com_mas"]=float(co[0]*1000)
            out["north_com_mas"]=float(co[1]*1000)
        else:
            out["east_com_mas"]=np.nan;out["north_com_mas"]=np.nan
        out["reference_file"]=anchor
        out["registration_scatter_east_mas"]=fit["scatter_east_mas"]
        out["registration_scatter_north_mas"]=fit["scatter_north_mas"]
        reg.append(out)
        quality.append({"filename":name,"filter":row["filter"],"mjd":row.mjd,
                        "reference_file":anchor,"method":"GALAXY_ONLY_AFFINE_LOO" if heldout else "GALAXY_ONLY_AFFINE",
                        "n_match":fit["n_common"],"n_used":fit["n_fit"],
                        "scatter_east_mas":fit["scatter_east_mas"],
                        "scatter_north_mas":fit["scatter_north_mas"],
                        "color_corr_east_r":np.nan,"color_corr_east_p":np.nan,
                        "color_corr_north_r":np.nan,"color_corr_north_p":np.nan})
        if heldout is not None and not controls.empty:
            selected=controls[controls.galaxy_id.eq(int(heldout))]
            if len(selected)==1:
                sc=selected.iloc[0]
                pos=transform(coef,float(sc[XY[0]]),float(sc[XY[1]]))*1000
                src_rows.append({"galaxy_id":int(heldout),"filename":name,
                                 "mjd":float(row.mjd),"epoch_label":str(row.epoch_label),
                                 "east_mas":float(pos[0]),"north_mas":float(pos[1]),
                                 "snr":float(sc.snr),
                                 "n_other_galaxies":int(fit["n_fit"])})
    registered=pd.DataFrame(reg)
    q=pd.DataFrame(quality)
    sources=pd.DataFrame(src_rows)
    summary={"status":"REGISTERED" if not registered.empty else "NO_REGISTERED_FRAMES",
             "anchor_file":anchor,"n_input":len(measurements),
             "n_registered":len(registered),"n_early_registered":int((registered.epoch_label=="early").sum()) if not registered.empty else 0,
             "n_late_registered":int((registered.epoch_label=="late").sum()) if not registered.empty else 0,
             "skipped_exposures":skipped}
    if not registered.empty and set(registered.epoch_label)=={"early","late"}:
        try:
            pm,_=infer_proper_motion(registered,q,"CROSS_FILTER_JWST")
            summary["target_refit_pm"]={
                k:pm.get(k) for k in ("pm_status","classification","mu_alpha_cosdec_masyr",
                                      "mu_delta_masyr","sigma_mu_alpha_cosdec_masyr",
                                      "sigma_mu_delta_masyr","baseline_days","significance_2d",
                                      "centroid_method_disagreement")
            }
            summary["target_refit_is_independent_of_frozen_paper"]=False
            summary["target_significance_not_calibrated"]=True
        except Exception as exc:
            summary["target_refit_error"]=str(exc)
    if heldout is not None:
        null=combined_epoch_measurement(sources,registered)
        summary["heldout_galaxy_null"]=null
        summary["galaxy_appears_in_all_registered_target_exposures"]=bool(
            not registered.empty and len(sources)==len(registered))
        summary["galaxy_appears_in_all_original_target_usable_exposures"]=bool(
            len(sources)==len(measurements) and len(registered)==len(measurements))
    return registered,q,sources,summary


def summarize_full_epoch_null(results: list[dict]):
    """Assess SOURCE-level tails only, not correlated dithers or a 5-sigma test."""
    full=[x for x in results if x.get("galaxy_appears_in_all_original_target_usable_exposures")
          and x.get("heldout_galaxy_null",{}).get("status")=="FULL_EPOCH_GALAXY_NULL_MEASURED"]
    partial=[x for x in results if x.get("heldout_galaxy_null",{}).get("status")=="FULL_EPOCH_GALAXY_NULL_MEASURED"]
    info={"n_heldout_attempts":len(results),"n_full_exposure_galaxies":len(full),
          "n_at_least_two_epoch_galaxies":len(partial),"n_confirmed_stationary":0,
          "null_calibrated":False}
    if len(full)>=3:
        xy=np.array([[x["heldout_galaxy_null"]["delta_east_mas"],
                      x["heldout_galaxy_null"]["delta_north_mas"]] for x in full],float)
        bias=np.median(xy,axis=0)
        cov=np.cov(xy,rowvar=False,ddof=1)
        info.update({"status":"PRELIMINARY_FULL_EPOCH_GALAXY_DIAGNOSTIC",
                     "median_east_north_displacement_mas":bias.tolist(),
                     "scatter_east_north_displacement_mas":np.sqrt(np.diag(cov)).tolist(),
                     "covariance_mas2":cov.tolist(),
                     "minimum_tail_resolution":1/(len(xy)+1)})
    else:
        info["status"]="INSUFFICIENT_FULL_EXPOSURE_GALAXY_NULLS"
    info["warning"]="Only photometric galaxies; small source N. No validated stationarity, chromatic ePSF matching, independent WCS bias separation, frozen-PM reproduction, or calibrated sigma."
    return info
