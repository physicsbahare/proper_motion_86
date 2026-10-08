#!/usr/bin/env python3
"""Empirical effective PSF (ePSF) centroid audit from ORIGINAL JWST CAL exposures.

A local empirical PSF is built for each original JWST NIRCam Stage-2 exposure
using isolated high-S/N field controls IN THE SAME detector and filter.
No PSF is imported from another field; no data/motion tables are overwritten.

This is a diagnostic comparison, not an independent third epoch nor
a jointly constrained proper-motion/parallax fit.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter
from astropy.nddata import NDData
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from photutils.psf import EPSFBuilder, extract_stars

from pm86.archive import load_covering_cutouts
from pm86.measurement import measure_exposure, masks, local_pixel_to_sky, sky_to_tangent
from defensible.exhaustive_pm import measure_exposure_dq_aware

CANDIDATE_ID=1651050
RA=150.307303724545
DEC=2.27642841317716
BAD_DQ=1|2|4|16
CAL_FILENAME_SUFFIXES=(
    "jw01727065001_04101_00002_nrcblong_cal.fits",
    "jw01727065001_04101_00003_nrcblong_cal.fits",
    "jw05893016008_05101_00001_nrcblong_cal.fits",
)


def isolated_reference_sources(controls,sci,err,dq,tx,ty,min_snr=18.,min_n=4):
    """Explicit star-selection flags, includes isolation and DQ apertures.

    A high-S/N control is not automatically a PSF star; reject nearby
    comparators, edge cutouts, bad DQ cores, extended or double-peaked cores.
    """
    data=np.asarray(sci,dtype=float)
    rows=[]
    for _,r in controls.iterrows():
        x=float(r["x_local"]);y=float(r["y_local"]);snr=float(r["snr"])
        ix=int(round(x));iy=int(round(y))
        if min(ix,iy)<13 or ix>=data.shape[1]-13 or iy>=data.shape[0]-13:
            continue
        if np.hypot(x-tx,y-ty)<18 or snr<min_snr:
            continue
        neighbors=controls.loc[controls.index!=r.name]
        if len(neighbors):
            sep=np.hypot(
                neighbors["x_local"].to_numpy(float)-x,
                neighbors["y_local"].to_numpy(float)-y)
            if np.nanmin(sep)<19:
                continue
        sy=slice(iy-6,iy+7);sx=slice(ix-6,ix+7)
        sub=data[sy,sx];es=err[sy,sx];d=dq[sy,sx]
        if not (np.isfinite(sub).all() and np.isfinite(es).all() and (es>0).all()
                and np.count_nonzero((np.asarray(d,dtype=np.uint64)&BAD_DQ)!=0)==0):
            continue
        yy,xx=np.indices(sub.shape)
        rring=np.hypot(xx-6,yy-6)
        bg=np.nanmedian(sub[rring>4.5])
        positive=np.clip(sub-bg,0,None)
        core=positive[rring<=3.5].sum()
        if core<=0:
            continue
        # A basic morphology/isolation screen; not a star/galaxy classifier.
        pix_radius=np.sqrt(np.sum(positive[rring<=3.5]*rring[rring<=3.5]**2)/core)
        if not (.4<pix_radius<2.6):continue
        rows.append({"x":x,"y":y,"snr":snr,
                     "core_rms_radius_pix":pix_radius,"flux":float(r["flux"])})
    df=pd.DataFrame(rows).sort_values("snr",ascending=False) if rows else pd.DataFrame()
    if len(df)<min_n:raise RuntimeError(
      f"Only {len(df)} isolated SNR>{min_snr} reference sources (need {min_n}); "
      "cannot construct a defensible empirical PSF. NO FALLBACK TO GAUSSIAN.")
    return df.head(12).reset_index(drop=True)


def build_local_epsf(sci,err,dq,ref,oversampling=4):
    bad=~np.isfinite(sci)|~np.isfinite(err)|(err<=0)|((np.asarray(dq,dtype=np.uint64)&BAD_DQ)!=0)
    _,sky,_=sigma_clipped_stats(sci,mask=bad,sigma=3)
    data=np.asarray(sci,dtype=float)-float(sky)
    nd=NDData(data,mask=bad)
    cat=Table()
    cat["x"]=ref.x.to_numpy(float);cat["y"]=ref.y.to_numpy(float)
    stars=extract_stars(nd,cat,size=17)
    if stars.n_stars<4:raise RuntimeError(
        f"Only {stars.n_stars} valid ePSF reference stars after extraction")
    builder=EPSFBuilder(oversampling=oversampling,fit_shape=(9,9),
                        recentering_boxsize=(7,7),
                        maxiters=4,fitter_maxiters=100,progress_bar=False)
    with warnings.catch_warnings(record=True) as cw:
        warnings.simplefilter("always")
        output=builder(stars)
    epsf=output.epsf
    if not np.isfinite(epsf.data).all() or not np.any(epsf.data>0):
        raise ValueError("ePSF contains invalid values")
    return epsf,data,{"n_input_stars":int(stars.n_stars),
                      "n_good_stars":int(getattr(output.fitted_stars,"n_good_stars",stars.n_stars)),
                      "converged":bool(getattr(output,"converged",False)),
                      "iterations":int(getattr(output,"iterations",0)),
                      "builder_warnings":[str(w.message) for w in cw[:6]]}


def model_pixels(epsf,xx,yy,cx,cy,flux,bg):
    m=epsf.copy()
    m.x_0.value=float(cx)
    m.y_0.value=float(cy)
    m.flux.value=float(flux)
    return np.asarray(m(xx,yy),dtype=float)+float(bg)


def psf_target_fit(epsf,data,err,dq,tx,ty,seed=None,half=5):
    """Fitted ePSF flux and centroid with bounded source location."""
    Y,X=np.indices(data.shape)
    keep=(np.abs(X-tx)<=half)&(np.abs(Y-ty)<=half)&np.isfinite(data)&np.isfinite(err)&(err>0)&((dq.astype(np.uint64)&BAD_DQ)==0)
    x=X[keep].astype(float);y=Y[keep].astype(float)
    obs=data[keep].astype(float);sigma=err[keep].astype(float)
    if len(x)<40:raise RuntimeError("Too few good pixels for 11x11 PSF fit")
    # Local image scaling: a PSF template has unit integrated flux.
    ix=int(round(tx));iy=int(round(ty))
    peak=float(np.nanmax(data[iy-2:iy+3,ix-2:ix+3]))
    integral=max(peak*5,1e-4)
    x0,y0=seed if seed is not None and np.isfinite(seed).all() else (tx,ty)
    x0=float(np.clip(x0,tx-2.49,tx+2.49));y0=float(np.clip(y0,ty-2.49,ty+2.49))
    guess=np.array([x0,y0,integral,0.0],float)
    low=[tx-2.5,ty-2.5,0.0,-np.inf]
    high=[tx+2.5,ty+2.5,np.inf,np.inf]
    def residual(p):
        return (model_pixels(epsf,x,y,*p)-obs)/sigma
    opt=least_squares(residual,guess,bounds=(low,high),max_nfev=500,loss="linear")
    cx,cy,flux,bg=opt.x
    boundary=any(abs(opt.x[i]-[low[0],low[1]][i])<.02 or abs(opt.x[i]-[high[0],high[1]][i])<.02 for i in (0,1))
    redchi=float(np.sum(opt.fun**2)/max(1,len(x)-len(opt.x)))
    # only a *formal* uncertainty: assumes fixed ePSF, uncorrelated pixels,
    # no PSF mismatch, no registration uncertainty.
    try:
        cov=np.linalg.pinv(opt.jac.T@opt.jac)
        errx,erry=np.sqrt(np.maximum(np.diag(cov)[:2],0))
    except Exception:
        errx=erry=np.nan
    return {
        "x_epsf":float(cx),"y_epsf":float(cy),
        "flux_epsf":float(flux),"bg_epsf":float(bg),
        "formal_sigma_x_pix":float(errx),"formal_sigma_y_pix":float(erry),
        "reduced_chi2_epsf":redchi,
        "optimizer_success":bool(opt.success),
        "fit_at_centroid_bound":bool(boundary),
        "n_pixels":int(len(x)),
        "fit_halfsize":int(half),
    }


def target_diagnostics(exp,epsf,data,initial,fit,outdir):
    tx=exp.x_full-exp.x0;ty=exp.y_full-exp.y0
    # direct comparison to original adopted 2DG using embedded gWCS
    xg=float(initial["x_2dg"]);yg=float(initial["y_2dg"])
    ra_epsf,dec_epsf=local_pixel_to_sky(exp,fit["x_epsf"],fit["y_epsf"])
    ra_gauss,dec_gauss=local_pixel_to_sky(exp,xg,yg)
    e,n=sky_to_tangent(ra_epsf,dec_epsf,ra_gauss,dec_gauss)
    delta_e=float(e*1000);delta_n=float(n*1000)
    fig,ax=plt.subplots(1,4,figsize=(15,3.9))
    ix=int(round(tx));iy=int(round(ty));half=11
    y0,y1=iy-half,iy+half+1;x0,x1=ix-half,ix+half+1
    zz=data[y0:y1,x0:x1]
    yy,xx=np.indices(zz.shape);xx=xx+x0;yy=yy+y0
    ps=model_pixels(epsf,xx,yy,fit["x_epsf"],fit["y_epsf"],
                    fit["flux_epsf"],fit["bg_epsf"])
    err=np.asarray(exp.err,dtype=float)[y0:y1,x0:x1]
    obs_norm=np.nanpercentile(zz,99)
    for a,v,title in zip(ax[:3],[zz,ps,(zz-ps)/err],
                        ["Original CAL SCI - background","Fitted local ePSF + background",
                         "Residual / ERR"]):
        a.imshow(v,origin="lower",cmap="RdBu_r" if "Residual" in title else "gray",
                 vmin=-5 if "Residual" in title else np.nanpercentile(zz,3),
                 vmax=5 if "Residual" in title else obs_norm)
        a.set_title(title,fontsize=9)
        if "Residual" not in title:
            a.scatter([tx-x0],[ty-y0],marker="x",c="lime",s=70,label="Catalog")
            a.scatter([xg-x0],[yg-y0],marker="o",facecolors="none",edgecolors="cyan",s=75,label="2DG")
            a.scatter([fit["x_epsf"]-x0],[fit["y_epsf"]-y0],marker="+",c="orange",s=85,label="ePSF")
            a.legend(fontsize=7)
    ax[3].imshow(epsf.data,origin="lower",cmap="viridis")
    ax[3].set_title("Same-exposure empirical ePSF",fontsize=9)
    fig.suptitle(
      f"{exp.filter_name} {exp.filename} | {len(data)} px original CAL\n"
      f"ePSF-Gaussian delta: east {delta_e:+.1f}, north {delta_n:+.1f} mas | "
      f"ePSF reduced chi2 {fit['reduced_chi2_epsf']:.2f}",fontsize=10)
    fig.tight_layout(rect=[0,0,1,.88])
    p=outdir/f"{exp.filename.replace('_cal.fits','')}_epsf_review.png"
    fig.savefig(p,dpi=160,bbox_inches="tight");plt.close(fig)
    return {"delta_east_epsf_minus_gauss_mas":delta_e,
            "delta_north_epsf_minus_gauss_mas":delta_n,
            "delta_pix_epsf_minus_gauss":float(np.hypot(fit["x_epsf"]-xg,fit["y_epsf"]-yg)),
            "qa_panel":str(p),"gaussian_x":xg,"gaussian_y":yg,
            "catalog_x_local":float(tx),"catalog_y_local":float(ty)}


def fit_exposure(rec,root):
    filename=str(rec.filename)
    filt=str(rec["filter"])
    products=pd.DataFrame([{"productFilename":filename,
                           "dataURI":f"mast:JWST/product/{filename}",
                           "productType":"SCIENCE"}])
    cutouts, audit=load_covering_cutouts(
      CANDIDATE_ID,RA,DEC,"epsf",filt,products)
    if len(cutouts)!=1:raise RuntimeError(f"CAL exposure inaccessible: {audit}")
    exp=cutouts[0]
    initial,controls=measure_exposure_dq_aware(exp)
    if not bool(initial.get("astrometrically_usable",False)):
        raise RuntimeError("Previously accepted target now fails DQ/centroid quality; abort ePSF")
    tx=exp.x_full-exp.x0;ty=exp.y_full-exp.y0
    ref=isolated_reference_sources(controls,exp.sci,exp.err,exp.dq,tx,ty)
    ref.to_csv(root/"reference_stars"/f"{filename}_stars.csv",index=False)
    epsf,data,binfo=build_local_epsf(exp.sci,exp.err,exp.dq,ref)
    gauss=np.array([float(initial["x_2dg"]),float(initial["y_2dg"])])
    fit=psf_target_fit(epsf,data,np.asarray(exp.err,dtype=float),
                       np.asarray(exp.dq,dtype=np.uint32),tx,ty,
                       seed=gauss,half=5)
    if fit["fit_at_centroid_bound"] or not fit["optimizer_success"]:
        raise RuntimeError("Empirical PSF fit failed or hit astrometric boundary")
    alt=psf_target_fit(epsf,data,np.asarray(exp.err,dtype=float),
                       np.asarray(exp.dq,dtype=np.uint32),tx,ty,
                       seed=gauss,half=4)
    diff=float(np.hypot(alt["x_epsf"]-fit["x_epsf"],alt["y_epsf"]-fit["y_epsf"]))
    meta=target_diagnostics(exp,epsf,data,initial,fit,root/"panels")
    np.savez_compressed(root/"models"/f"{filename}.npz",
        epsf_data=np.asarray(epsf.data,dtype=np.float32),
        epsf_oversampling=int(4),
        source_filename=filename,
        reference_stars_used=ref[["x","y"]].to_numpy(float),
        measured_epsf_xy=np.array([fit["x_epsf"],fit["y_epsf"]]))
    return {"filename":filename,"filter":filt,"mjd":float(exp.mjd),
        "n_reference_stars":len(ref),"epsf_builder_converged":binfo["converged"],
        "epsf_builder_iterations":binfo["iterations"],
        "epsf_builder_warnings":" | ".join(binfo["builder_warnings"][:4]),
        "half4_vs_half5_centroid_shift_pix":diff,
        "raw_gaussian_snr":float(initial["astrometric_snr"]),
        "dq_or_r4":int(initial["dq_or_r4"]),
        **fit,**meta}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",default="data/1651050_astrometric_epoch_manifest.csv")
    ap.add_argument("--out",default="publication_check/epsf_1651050")
    args=ap.parse_args()
    root=Path(args.out)
    for name in ("panels","models","reference_stars"):
        (root/name).mkdir(parents=True,exist_ok=True)
    manifest=pd.read_csv(args.input)
    accepted=manifest[manifest.astrometrically_usable.astype(bool)]
    expected=set(CAL_FILENAME_SUFFIXES)
    if set(accepted.filename)!=expected:
        raise RuntimeError(f"Input manifest mismatch: {set(accepted.filename)^expected}")
    results=[];errors=[]
    for _,rec in accepted.iterrows():
        try:
            print(f"\n>>> Building {rec['filter']} local empirical PSF: {rec.filename}",flush=True)
            data=fit_exposure(rec,root)
            results.append(data)
            print(json.dumps(data,indent=2,default=float),flush=True)
        except Exception as exc:
            msg=f"{type(exc).__name__}: {exc}"
            print("ERROR "+msg,flush=True)
            errors.append({"filename":str(rec.filename),"error":msg})
    if errors:
        pd.DataFrame(errors).to_csv(root/"exposure_errors.csv",index=False)
    if results:
        pd.DataFrame(results).to_csv(root/"epsf_centroid_comparison.csv",index=False)
    guide=[
        "# Empirical effective PSF comparison: candidate 1651050",
        "",
        "Uses original Stage-2 JWST CAL detector pixels and photutils EPSFBuilder. "
        "Each exposure has a separately built local empirical PSF from isolated field "
        "sources with clean DQ. F356W and F444W NEVER share a PSF model.",
        "",
        "The fitted target position is compared directly to the original 2DG "
        "centroid, and local WCS converts their difference to east and north "
        "milliarcseconds. This isolates the PSF-fitting method shift but does "
        "not independently register reference sources across epochs.",
        "",
        "The report includes a secondary fitting-window check using exactly the "
        "same empirical PSF, and counts the stars used to construct that PSF.",
        "",
        "The photutils effective PSF is an empirical image model; it is not a "
        "JWST WebbPSF optical simulation. Faint objects of extremely unusual color "
        "can retain chromatic PSF centroid biases relative to ordinary field stars.",
        "",
        "**Do not call the candidate confirmed solely on the basis of this report.** "
        "Cross-filter astrometry, unmodeled morphology/color effects and "
        "two-epoch parallax degeneracy still matter.",
    ]
    (root/"READ_ME_FIRST.md").write_text("\n".join(guide)+"\n")
    if len(results)!=3:
        raise RuntimeError(
            f"Only {len(results)}/3 accepted exposure ePSF centroid results. "
            "See exposure_errors.csv; insufficient for cross-filter QA.")


if __name__=="__main__":main()
