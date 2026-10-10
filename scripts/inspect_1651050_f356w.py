#!/usr/bin/env python3
"""Independent visual review of original JWST/NIRCam CAL pixels for 1651050.

For each exact CAL filename in the archived September-2026 astrometry:
- use its embedded distortion-aware JWST gWCS, not a mosaic/WCS guess;
- read only a local 401x401 SCI/ERR/DQ stamp by HTTP range from public MAST;
- rerun the repository's unmodified DQ-aware target measurement;
- show image, error-weighted S/N, bad pixels, model residual, context, profiles;
- export raw local SCI/ERR/DQ .npz and a reviewer decision sheet.
This is QA, not a new proper-motion or parallax fit. No science table is changed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
import numpy as np
import pandas as pd
from scipy.ndimage import maximum_filter
from astropy.stats import sigma_clipped_stats

from pm86.archive import load_covering_cutouts
from pm86.measurement import masks
from defensible.exhaustive_pm import measure_exposure_dq_aware

ID=1651050
RA=150.307303724545
DEC=2.27642841317716
BAD_BITS=1|2|4|16


def window2d(a, x, y, half):
    a=np.asarray(a)
    x0=max(0,int(round(x))-half)
    x1=min(a.shape[1],int(round(x))+half+1)
    y0=max(0,int(round(y))-half)
    y1=min(a.shape[0],int(round(y))+half+1)
    return a[y0:y1,x0:x1],x0,y0


def centroid_distance_pixels(row):
    q=np.array([row.get("x_2dg",np.nan),row.get("y_2dg",np.nan),
                row.get("x_com",np.nan),row.get("y_com",np.nan)],dtype=float)
    if not np.isfinite(q).all():
        return np.nan
    return float(np.hypot(q[0]-q[2],q[1]-q[3]))


def local_peaks(snr,x,y,radius=13,threshold=3.5):
    """Find candidate neighboring peaks for visual review, not automated ID matching."""
    a=np.asarray(snr,dtype=float)
    v=np.isfinite(a)&(a>=threshold)&(a==maximum_filter(a,size=5))
    yy,xx=np.where(v)
    d=np.hypot(xx-x,yy-y)
    keep=(d>2.5)&(d<radius)
    return sorted([(float(xx[i]),float(yy[i]),float(a[yy[i],xx[i]]))
                   for i in np.where(keep)[0]],key=lambda z:-z[2])[:8]


def best_offset(row,x,y):
    """Exact original coordinates; never silently substitute sky centroids."""
    vals=[]
    for label,kx,ky in (("Gauss","x_2dg","y_2dg"),("COM","x_com","y_com")):
        px=float(row.get(kx,np.nan));py=float(row.get(ky,np.nan))
        if np.isfinite(px) and np.isfinite(py):
            vals.append((label,px-x,py-y))
    return vals


def make_image(exp,row,outroot):
    sci=np.asarray(exp.sci,dtype=float)
    err=np.asarray(exp.err,dtype=float)
    dq=np.asarray(exp.dq,dtype=np.uint32)
    badphot,badastrom=masks(sci,err,dq)
    sky=float(sigma_clipped_stats(sci,mask=badphot)[1])
    data=sci-sky
    tx=exp.x_full-exp.x0
    ty=exp.y_full-exp.y0
    snr=np.divide(data,err,out=np.full_like(data,np.nan),where=err>0)
    bad=((dq.astype(np.uint64)&np.uint64(BAD_BITS))!=0)
    cent_sep=centroid_distance_pixels(row)
    peak_candidates=local_peaks(snr,tx,ty)
    h=20
    z,x0,y0=window2d(data,tx,ty,h)
    sz,_,_=window2d(snr,tx,ty,h)
    dz,_,_=window2d(dq,tx,ty,h)
    bz,_,_=window2d(bad,tx,ty,h)
    wide,wx,wy=window2d(data,tx,ty,34)
    xx,yy=np.meshgrid(np.arange(x0,x0+z.shape[1]),np.arange(y0,y0+z.shape[0]))
    gauss=np.full_like(z,np.nan)
    if bool(row.get("forced_fit_accepted",False)):
        px=float(row.get("x_2dg",np.nan));py=float(row.get("y_2dg",np.nan))
        sx=float(row.get("forced_sigma_x_pix",np.nan))
        sy=float(row.get("forced_sigma_y_pix",np.nan))
        amp=float(row.get("forced_amp",np.nan))
        if np.isfinite([px,py,sx,sy,amp]).all() and sx>0 and sy>0:
            gauss=amp*np.exp(-0.5*((xx-px)/sx)**2-0.5*((yy-py)/sy)**2)
    residual=np.divide(z-gauss,err[y0:y0+z.shape[0],x0:x0+z.shape[1]],
                       out=np.full_like(z,np.nan),where=err[y0:y0+z.shape[0],x0:x0+z.shape[1]]>0)
    # Store original pixels with an unambiguous full-detector coordinate mapping.
    near_sci,na,nb=window2d(sci,tx,ty,34)
    near_err,_,_=window2d(err,tx,ty,34)
    near_dq,_,_=window2d(dq,tx,ty,34)
    stem=exp.filename.replace("_cal.fits","")
    np.savez_compressed(outroot/"pixel_stamps"/f"{stem}.npz",
        SCI=near_sci.astype(np.float32), ERR=near_err.astype(np.float32),
        DQ=near_dq.astype(np.uint32), pixel_window_x0=int(na),
        pixel_window_y0=int(nb), full_detector_window_x0=int(na+exp.x0),
        full_detector_window_y0=int(nb+exp.y0),
        catalog_x_local=float(tx),catalog_y_local=float(ty),
        gaussian_x_local=float(row.get("x_2dg",np.nan)),
        gaussian_y_local=float(row.get("y_2dg",np.nan)),
        com_x_local=float(row.get("x_com",np.nan)),
        com_y_local=float(row.get("y_com",np.nan)),
        ra_deg=float(RA),dec_deg=float(DEC),mjd=float(exp.mjd),
        filename=str(exp.filename),filter=str(exp.filter_name))
    fig,axes=plt.subplots(2,3,figsize=(15,10))
    vmax=np.nanpercentile(z[np.isfinite(z)],98) if np.isfinite(z).any() else 1.
    norm=AsinhNorm(linear_width=max(abs(vmax)*.08,1e-6),
                   vmin=np.nanpercentile(z,5),vmax=max(vmax,1e-5))
    axes[0,0].imshow(z,origin="lower",cmap="gray",norm=norm)
    axes[0,1].imshow(sz,origin="lower",cmap="RdBu_r",vmin=-6,vmax=8)
    axes[0,2].imshow(bz.astype(int),origin="lower",cmap="Reds",vmin=0,vmax=1)
    if np.isfinite(residual).any():
        axes[1,0].imshow(residual,origin="lower",cmap="RdBu_r",vmin=-5,vmax=5)
    else:
        axes[1,0].text(.5,.5,"No accepted forced Gaussian model",
                        ha="center",va="center",transform=axes[1,0].transAxes)
    wnorm=AsinhNorm(linear_width=max(np.nanpercentile(abs(wide),90)*.1,1e-6),
                     vmin=np.nanpercentile(wide,4),vmax=max(np.nanpercentile(wide,99.5),1e-5))
    axes[1,1].imshow(wide,origin="lower",cmap="gray",norm=wnorm)
    for ax,title in zip(axes.flat,
      ["Background-subtracted SCI (41 px)","(SCI - sky) / ERR","DQ astrometric mask (red = rejected)",
       "Forced Gaussian residual / ERR","Wider field (69 px)","Central 1D pixel profiles"]):
        ax.set_title(title,fontsize=11)
    for ax,origin in [(axes[0,0],(x0,y0)),(axes[0,1],(x0,y0)),
                      (axes[0,2],(x0,y0)),(axes[1,0],(x0,y0)),
                      (axes[1,1],(wx,wy))]:
        ax.scatter([tx-origin[0]],[ty-origin[1]],marker="x",s=90,
                   c="lime",linewidths=1.7,label="Catalog/WCS")
        positions=[
            ("Gaussian/forced","x_2dg","y_2dg","cyan","o"),
            ("Center of mass","x_com","y_com","yellow","+"),
        ]
        for name,kx,ky,color,marker in positions:
            x=float(row.get(kx,np.nan));y=float(row.get(ky,np.nan))
            if np.isfinite(x) and np.isfinite(y):
                ax.scatter([x-origin[0]],[y-origin[1]],marker=marker,
                           s=95,facecolors="none" if marker=="o" else color,
                           edgecolors=color if marker=="o" else None,
                           color=color if marker!="o" else None,
                           linewidths=1.6,label=name)
        if ax is axes[1,1]:
            for x,y,s in peak_candidates:
                ax.plot(x-origin[0],y-origin[1],marker="s",markersize=6,
                        markeredgecolor="magenta",markerfacecolor="none")
    axes[0,0].legend(loc="lower right",fontsize=7,framealpha=.7)
    iy=int(np.clip(round(ty-y0),0,z.shape[0]-1))
    ix=int(np.clip(round(tx-x0),0,z.shape[1]-1))
    axes[1,2].plot(np.arange(z.shape[1])+x0-tx,z[iy,:],label="Row through WCS")
    axes[1,2].plot(np.arange(z.shape[0])+y0-ty,z[:,ix],label="Column through WCS")
    if np.isfinite(gauss).any():
        axes[1,2].plot(np.arange(z.shape[1])+x0-tx,gauss[iy,:],"--",label="Fitted Gaussian row")
    axes[1,2].axvline(0,c="gray",linestyle=":")
    axes[1,2].set_xlabel("Offset from WCS position (pixels)")
    axes[1,2].set_ylabel("Background-subtracted SCI units")
    axes[1,2].legend(fontsize=8)
    subtitle=(
        f"{exp.filter_name} | {stem} | MJD {exp.mjd:.5f} | "
        f"usable {bool(row.get('astrometrically_usable',False))} | "
        f"S/N {float(row.get('astrometric_snr',np.nan)):.2f} | "
        f"Gaussian-COM {cent_sep:.2f} px"
    )
    fig.suptitle(f"COSMOS-Web target {ID}: Original calibrated pixels\n{subtitle}",fontsize=12)
    fig.text(.03,.015,"'WCS' is fixed catalog sky position, not a proper-motion prediction. "
             "Magenta squares flag local bright peaks, not confirmed neighbors. "
             "Visual approval requires independent confirmation.",fontsize=9)
    fig.tight_layout(rect=[0,.035,1,.93])
    image_path=outroot/"panels"/f"{stem}_review.png"
    fig.savefig(image_path,dpi=150,bbox_inches="tight")
    plt.close(fig)
    center_dq=int(dq[int(round(ty)),int(round(tx))])
    qa={
        "candidate_id":ID,"filter":exp.filter_name,"filename":exp.filename,
        "mjd":float(exp.mjd),"astrometrically_usable":bool(row.get("astrometrically_usable",False)),
        "astrometry_source":str(row.get("astrometry_source","")),
        "raw_snr":float(row.get("raw_snr_err",np.nan)),
        "clean_snr":float(row.get("clean_snr_err",np.nan)),
        "astrometric_snr":float(row.get("astrometric_snr",np.nan)),
        "forced_fit_accepted":bool(row.get("forced_fit_accepted",False)),
        "forced_chi2_reduced":float(row.get("forced_reduced_chi2",np.nan)),
        "centroid_2dg_com_distance_pix":cent_sep,
        "pixel_scale_mas":float(row.get("pixel_scale_mas",np.nan)),
        "centroid_2dg_com_distance_mas":cent_sep*float(row.get("pixel_scale_mas",np.nan)),
        "gaussian_vs_catalog_pix":float(np.hypot(
            float(row.get("x_2dg",np.nan))-tx,float(row.get("y_2dg",np.nan))-ty)),
        "dq_at_catalog":center_dq,
        "dq_or_r4":int(row.get("dq_or_r4",0)),
        "dq_bad_at_gaussian_centroid":None,
        "n_nearby_snr_peaks_13pix":len(peak_candidates),
        "nearest_peak_proximity_pix":float(min((np.hypot(x-tx,y-ty)
                for x,y,_ in peak_candidates),default=np.nan)),
        "panel":str(image_path),
        "stamp_file":str(outroot/"pixel_stamps"/f"{stem}.npz"),
    }
    gx=float(row.get("x_2dg",np.nan));gy=float(row.get("y_2dg",np.nan))
    if np.isfinite(gx) and np.isfinite(gy):
        cx,cy=int(round(gx)),int(round(gy))
        if 0<=cx<dq.shape[1] and 0<=cy<dq.shape[0]:
            qa["dq_bad_at_gaussian_centroid"]=bool(int(dq[cy,cx])&BAD_BITS)
    qa["requires_visual_review"]=bool(
        not qa["astrometrically_usable"] or
        (np.isfinite(cent_sep) and cent_sep>0.5) or
        (qa["dq_or_r4"]&BAD_BITS)!=0 or
        len(peak_candidates)>0
    )
    return qa


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",default="data/1651050_targeted_review_products.csv")
    ap.add_argument("--output",default="publication_check/centroid_visual")
    args=ap.parse_args()
    out=Path(args.output)
    for sub in ("panels","pixel_stamps"):
        (out/sub).mkdir(parents=True,exist_ok=True)
    inp=pd.read_csv(args.input)
    if len(inp)!=6 or inp["filename"].duplicated().any() or sorted(inp["filter"].value_counts().to_dict().items())!=[("F356W",4),("F444W",2)]:
        raise ValueError("Unexpected source audit list; expected 4 F356W and 2 F444W unique CAL files")
    summaries=[];errors=[]
    for _,rec in inp.iterrows():
        filt=str(rec["filter"]);filename=str(rec["filename"])
        print(f"Downloading original CAL cutout: {filt} {filename}",flush=True)
        df=pd.DataFrame([{
            "productFilename":filename,
            "dataURI":f"mast:JWST/product/{filename}",
            "productType":"SCIENCE"}])
        try:
            cutouts,audits=load_covering_cutouts(ID,RA,DEC,"review",filt,df)
            if len(cutouts)!=1:
                raise RuntimeError(f"Original CAL not covering target: {audits}")
            exp=cutouts[0]
            row,controls=measure_exposure_dq_aware(exp)
            qa=make_image(exp,row,out)
            summaries.append(qa)
            print(json.dumps(qa,allow_nan=True,default=str),flush=True)
        except Exception as exc:
            message=f"{type(exc).__name__}: {exc}"
            errors.append({"filter":filt,"filename":filename,"error":message})
            print("EXPOSURE FAILED: "+message,flush=True)
    summary=pd.DataFrame(summaries)
    if len(summary):
        summary.to_csv(out/"centroid_quality_summary.csv",index=False)
        review=summary[["filter","filename","mjd","astrometrically_usable",
                        "centroid_2dg_com_distance_pix","dq_or_r4","panel"]].copy()
        for col in ("target_visible","gaussian_center_on_target",
                    "neighbor_or_blend_problem","bad_pixels_affect_target",
                    "accept_for_astrometry","reviewer_notes"):
            review[col]=""
        review.to_csv(out/"MANUAL_REVIEW_CHECKLIST.csv",index=False)
    if errors:
        pd.DataFrame(errors).to_csv(out/"download_errors.csv",index=False)
    report=[
        "# Candidate 1651050 F356W centroid and neighboring-source visual QA",
        "",
        "This archive uses individual original JWST NIRCam Stage-2 CAL "
        "SCI/ERR/DQ images; it is not a mosaic-only simulated drawing.",
        "",
        f"Generated {len(summaries)} of {len(inp)} exposure panels.",
        "",
        "## Mandatory visual questions",
        "1. Is there a compact source at the green catalog/WCS cross in EACH exposure?",
        "2. Does the cyan forced-Gaussian circle lie on that source, rather than a neighbor?",
        "3. Is the yellow center-of-mass sign displaced because of visible blending, "
        "an extended wing, detector structure, or background noise?",
        "4. Are red-marked DQ pixels sufficiently close to bias the flux centroid?",
        "5. Does the target appear consistently in the other F356W exposures even "
        "if their automatic centroid fits were rejected?",
        "6. Compare accepted F356W centroids with both accepted F444W panels. "
        "Different filters have different PSFs, and images are centered on the fixed "
        "catalog WCS position, not aligned to expected proper motion.",
        "",
        "A reviewer must explicitly sign off on the CSV checklist. "
        "A visually convincing point source does not alone confirm proper motion.",
        "",
        "### Existing known issue",
        "The September 17 archived accepted F356W exposure showed "
        "2.18 pixel (approximately 138 mas) separation between the two centroid "
        "methods, 3 other F356W covering CALs were unusable for astrometry, "
        "and a JUMP_DET pixel was flagged within 4 pix of the catalog position. "
        "Current panels independently recompute these measurements.",
        "",
        "No detected source in a different filter is proof of a mover, "
        "and non-detections never imply a stationary source.",
    ]
    (out/"READ_ME_FIRST.md").write_text("\n".join(report)+"\n")
    if len(summaries) < len(inp):
        raise RuntimeError(f"Only {len(summaries)}/{len(inp)} images retrieved; inspect download_errors.csv")


if __name__=="__main__":
    main()
