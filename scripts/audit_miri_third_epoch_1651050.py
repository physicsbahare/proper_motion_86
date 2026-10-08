#!/usr/bin/env python3
"""Targeted original MIRI CAL cutout/forced photometry audit for 1651050.

Uses five public MAST MIRI observation IDs from the archived Sep 2026
inventory. Query new archive search separately for newly released data.
These are tests of target DETECTABILITY only, not a proper-motion solution.
A robust third astrometric epoch would require an unambiguous centroid,
independent local field registration, and appropriate MIRI chromatic PSF errors.
"""
from __future__ import annotations
import argparse,json,re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from astroquery.mast import Observations
from astropy.stats import sigma_clipped_stats

from pm86.archive import load_covering_cutouts, normalise_filter

CID=1651050
RA=150.307303724545
DEC=2.27642841317716
OBS={
 "200632589":"F770W",  # 2024 F770W
 "266380940":"F1000W", # 2025 F1000W
 "266380942":"F2100W", # 2025 F2100W
 "368595154":"F1000W", # later MIRI observations
 "368595153":"F2100W",
}
FILTER_APERTURES={"F770W":3.0,"F1000W":3.5,"F2100W":5.0}
BAD_DQ=1|2|4|16


def measure_forced_aperture(sci,err,dq,x,y,filter_name):
    sci=np.asarray(sci,dtype=float);err=np.asarray(err,dtype=float)
    dq=np.asarray(dq,dtype=np.uint64)
    yy,xx=np.indices(sci.shape,dtype=float)
    rr=np.hypot(xx-x,yy-y)
    r=FILTER_APERTURES[filter_name]
    valid=np.isfinite(sci)&np.isfinite(err)&(err>0)&((dq&BAD_DQ)==0)
    ring=valid&(rr>=max(12,2.5*r))&(rr<=max(23,4*r))
    aperture=valid&(rr<=r)
    if ring.sum()<35 or aperture.sum()<10:
        return {"status":"INSUFFICIENT_LOCAL_VALID_PIXELS"}
    bg=float(np.nanmedian(sci[ring]))
    residual=sci-bg
    total=float(np.sum(residual[aperture]))
    per_pixel_unc=np.sqrt(np.sum(err[aperture]**2))
    robust_sky=float(1.4826*np.nanmedian(np.abs(residual[ring]-np.median(residual[ring]))))
    sky_unc=robust_sky*np.sqrt(aperture.sum())
    # conservative of the two estimates: empirical scatter may include structured background
    max_err=max(per_pixel_unc,sky_unc)
    pixel_peak=float(np.nanmax(residual[valid&(rr<4)]))
    return dict(status="FORCED_APERTURE_ONLY",
        aperture_radius_pix=r,aperture_pixels=int(aperture.sum()),
        local_background=float(bg),net_aperture_flux_native=float(total),
        pipeline_flux_unc_native=float(per_pixel_unc),
        empirical_sky_unc_native=float(sky_unc),
        snr_pipeline_only=float(total/per_pixel_unc) if per_pixel_unc>0 else np.nan,
        snr_conservative=float(total/max_err) if max_err>0 else np.nan,
        local_peak_above_background_native=pixel_peak,
        bad_dq_pixels_aperture=int(np.count_nonzero((rr<=r)&~valid)),
        note="No aperture correction, no unit conversion, no independent centroid or registration; candidate may be displaced from catalog WCS by PM+parallax.")


def plot_exposure(exp,diagnostic,output):
    sci=np.asarray(exp.sci,float);err=np.asarray(exp.err,float)
    dq=np.asarray(exp.dq,dtype=np.uint64)
    tx=float(exp.x_full-exp.x0);ty=float(exp.y_full-exp.y0)
    cut_half=30
    xc=int(round(tx));yc=int(round(ty))
    x0=max(xc-cut_half,0);x1=min(xc+cut_half+1,sci.shape[1])
    y0=max(yc-cut_half,0);y1=min(yc+cut_half+1,sci.shape[0])
    sci=sci[y0:y1,x0:x1];err=err[y0:y1,x0:x1]
    dq=dq[y0:y1,x0:x1]
    good=np.isfinite(sci)&np.isfinite(err)&(err>0)&((dq&BAD_DQ)==0)
    bg=float(diagnostic.get("local_background",np.nanmedian(sci[good])))
    arr=sci-bg
    snr=np.divide(arr,err,out=np.full_like(arr,np.nan),where=good)
    fig,ax=plt.subplots(1,3,figsize=(14,5))
    finite=arr[np.isfinite(arr)]
    lo=float(np.nanpercentile(finite,2));hi=float(np.nanpercentile(finite,99.5))
    if hi<=lo: hi=lo+1e-6
    ax[0].imshow(arr,origin="lower",cmap="gray",norm=AsinhNorm(linear_width=max((hi-lo)/15,1e-6),vmin=lo,vmax=hi))
    ax[1].imshow(snr,origin="lower",cmap="RdBu_r",vmin=-5,vmax=7)
    ax[2].imshow(((dq&BAD_DQ)!=0).astype(float),origin="lower",cmap="Reds",vmin=0,vmax=1)
    for a,t in zip(ax,["MIRI CAL SCI minus local background","Per-pixel (SCI - background)/ERR","Masked DQ pixels"]):
        a.set_title(t,fontsize=11)
        x=tx-x0;y=ty-y0
        a.scatter([x],[y],marker="x",c="lime",s=85,lw=1.5)
        a.add_patch(plt.Circle((x,y),FILTER_APERTURES[exp.filter_name],fill=False,color="cyan",lw=1.2))
        a.set_xlim(max(0,x-22),min(sci.shape[1],x+23))
        a.set_ylim(max(0,y-22),min(sci.shape[0],y+23))
    fig.suptitle(f"Candidate {CID} | {exp.filter_name} | {exp.filename}\n"
      f"MJD {exp.mjd:.5f} | catalog/WCS green X | conservative forced S/N "
      f"{diagnostic.get('snr_conservative',float('nan')):.2f}\n"
      "NO astrometric centroid has been measured from this MIRI image",fontsize=11)
    fig.tight_layout(rect=[0,0,1,.88])
    stem=exp.filename.replace("_cal.fits","")
    out=output/"panels"/f"{stem}_MIRI_review.png"
    fig.savefig(out,dpi=160,bbox_inches="tight");plt.close(fig)
    np.savez_compressed(output/"pixel_stamps"/f"{stem}.npz",
       SCI=sci.astype(np.float32),ERR=err.astype(np.float32),DQ=dq.astype(np.uint32),
       stamp_x0_in_full_detector=int(exp.x0+x0),stamp_y0_in_full_detector=int(exp.y0+y0),
       catalog_x_in_stamp=float(tx-x0),catalog_y_in_stamp=float(ty-y0),
       filter=str(exp.filter_name),filename=str(exp.filename),mjd=float(exp.mjd),
       ra_deg=RA,dec_deg=DEC)
    return str(out)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",default="publication_check/miri_third_epoch")
    args=ap.parse_args()
    root=Path(args.output);root.mkdir(parents=True,exist_ok=True)
    for sub in ("panels","pixel_stamps"): (root/sub).mkdir(exist_ok=True)
    obsids=list(OBS)
    t=Observations.get_product_list(obsids).to_pandas()
    t=t[t.productFilename.astype(str).str.endswith("_mirimage_cal.fits")].copy()
    if "productType" in t:
        t=t[t.productType.astype(str).str.upper().eq("SCIENCE")]
    t=t.drop_duplicates("dataURI")
    t.to_csv(root/"MIRI_CANDIDATE_CAL_PRODUCTS.csv",index=False)
    print("MIRI CAL products discovered:",len(t),flush=True)
    if t.empty:
        raise RuntimeError("No calibrated MIRI science products returned for the archived five obs IDs")
    receipts=[]
    for _,row in t.iterrows():
        fname=str(row.productFilename);uri=str(row.dataURI)
        # Get the MIRI band from the actual product group rather than
        # applying the NIRCam three-digit filter regex to F1000W/F2100W.
        obsid=str(row.obsid)
        if obsid not in OBS:
            receipts.append(dict(filename=fname,obsid=obsid,status="OBS_GROUP_NOT_IDENTIFIED"))
            continue
        filt=OBS[obsid]
        print(f"MIRI WCS + original pixels: {obsid} {filt} {fname}",flush=True)
        try:
            exposures,cover=load_covering_cutouts(CID,RA,DEC,
                 "independent_MIRI_epoch",filt,pd.DataFrame([row]))
            if not exposures:
                st=cover[0].get("status","NO_COVERAGE") if cover else "NO_COVERAGE"
                receipts.append(dict(filename=fname,obsid=obsid,filter=filt,status=st))
                continue
            exp=exposures[0]
            # MIRI four-digit filter names are not represented by the repository's
            # NIRCam-centric normalization, so ensure our per-observation mapping.
            exp.filter_name=filt
            tx=exp.x_full-exp.x0;ty=exp.y_full-exp.y0
            measurement=measure_forced_aperture(exp.sci,exp.err,exp.dq,tx,ty,filt)
            panel=plot_exposure(exp,measurement,root) if measurement["status"]=="FORCED_APERTURE_ONLY" else ""
            receipts.append(dict(filename=fname,obsid=obsid,filter=filt,
                mjd=float(exp.mjd),detector=str(exp.detector),
                status="COVERS_TARGET",panel=panel,**measurement))
        except Exception as exc:
            receipts.append(dict(filename=fname,obsid=obsid,filter=filt,
                                 status="ANALYSIS_ERROR",reason=f"{type(exc).__name__}: {exc}"))
        pd.DataFrame(receipts).to_csv(root/"MIRI_ORIGINAL_CAL_COVERAGE_AND_FORCED_SNR.csv",index=False)
    d=pd.DataFrame(receipts)
    r=[
      "# Original MIRI epochs: candidate 1651050", "",
      "This is an explicitly *exploratory* forced-aperture and WCS coverage test.",
      "A covered CAL image or nominal aperture S/N is NOT an independently",
      "validated astrometric centroid. MIRI's wider PSF, structured backgrounds,",
      "different wavelengths and detector astrometric systematics need",
      "separate control-star registration and reliable matching.", "",
      f"Products considered {len(d)}; verified covering CAL products {(d.status=='COVERS_TARGET').sum()}.", "",
      "Any credible MIRI target detection should be followed with a",
      "PSF-model centroid, per-exposure field-registration null test,",
      "and cross-instrument distortion/color systematic budget, before adding",
      "it to a three-epoch position/proper-motion/parallax fit.", "",
      "All native MIRI CAL fluxes retain unconverted SCI units; conservative",
      "S/N checks use ERR and sky RMS but are not an official MIRI photometric",
      "calibration or aperture-corrected measurement. Non-detections are not",
      "evidence of zero proper motion.", "",
      "The source's 500-K atmosphere is intrinsically red; late MIRI",
      "images may give a qualitatively useful additional-epoch test.",
    ]
    (root/"READ_ME_FIRST.md").write_text("\n".join(r)+"\n")
    print("\n".join(r),flush=True)


if __name__=="__main__":main()
