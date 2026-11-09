#!/usr/bin/env python3
"""1651050 independent compact-field-control PSF audit, F444W and F356W.

A *bounded* validation of model-to-Gaussian centroid offsets in original CAL
detector images, using isolated high-SNR field control sources, excluding target.
No star/galaxy labels are inferred: compact morphology only. STPSF models are
spectral hypotheses, not empirically certified PSFs or fitted Teff.

Crucially: the PSF fit uses a correct 7x7 central stamp; target is excluded.
No paper PM or original acceptance labels are modified.
"""
from __future__ import annotations
import argparse,json,sys,traceback
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.ndimage import maximum_filter
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"first26"),str(ROOT/"scripts")]
import pm86.measurement as meas
from pm86.archive import query_candidate_inventory,get_pair_products,load_covering_cutouts
from run_cross_filter_null import choose_frozen_pair,EXPECTED_PAIRS,target_is_usable
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install
from audit_1651050_stpsf_neighbor import get_psf,native_template
from audit_1651050_f356w_joint import safe_out

CID=1651050
MAX_CONTROLS_PER_EXPOSURE=8
FIT_HALF=5
MIN_SNR=12.0
ISOLATION_PIX=10.0
PSF_TEMPS={"F444W":(900,5700),"F356W":(900,5700)}


def control_candidates(ex,controls,target_x,target_y):
    sci=np.asarray(ex.sci,float);err=np.asarray(ex.err,float)
    _,bad=meas.masks(sci,err,ex.dq)
    cc=controls.sort_values("snr",ascending=False)
    out=[]
    for _,row in cc.iterrows():
        x,y=float(row.x_local),float(row.y_local)
        if not np.isfinite([x,y,row.snr]).all() or row.snr<MIN_SNR:continue
        if np.hypot(x-target_x,y-target_y)<ISOLATION_PIX:continue
        if min(x,y,sci.shape[1]-1-x,sci.shape[0]-1-y)<FIT_HALF+3:continue
        if any(np.hypot(x-a,y-b)<ISOLATION_PIX for a,b in out):continue
        # Exclude sources with another detected control in 10-pixel radius.
        nearby=cc[(cc.x_local-x)**2+(cc.y_local-y)**2<ISOLATION_PIX**2]
        if len(nearby)>1:continue
        xi,yi=int(round(x)),int(round(y))
        patch_bad=bad[yi-2:yi+3,xi-2:xi+3]
        if patch_bad.shape!=(5,5) or patch_bad.any():continue
        out.append((x,y))
        if len(out)>=MAX_CONTROLS_PER_EXPOSURE:break
    return out


def fit_star(ex,x,y,psf):
    sci=np.asarray(ex.sci,float);err=np.asarray(ex.err,float)
    _,bad=meas.masks(sci,err,ex.dq)
    xi,yi=int(round(x)),int(round(y))
    ys=slice(yi-FIT_HALF,yi+FIT_HALF+1)
    xs=slice(xi-FIT_HALF,xi+FIT_HALF+1)
    z=sci[ys,xs];s=err[ys,xs];mask=bad[ys,xs]
    yy,xx=np.mgrid[ys,xs]
    good=(~mask)&np.isfinite(z)&np.isfinite(s)&(s>0)
    if good.sum()<75:return dict(ok=False,reason="too_few_unmasked_pixels")
    xx=xx[good].astype(float);yy=yy[good].astype(float)
    z=z[good];s=s[good]
    def linfit(cx,cy):
        template=native_template(psf,xx,yy,cx,cy)
        A=np.column_stack((template,np.ones(len(xx))))
        coef,*_=np.linalg.lstsq(A/s[:,None],z/s,rcond=None)
        residual=(z-A@coef)/s
        return residual,coef
    def residual(p):return linfit(p[0],p[1])[0]
    opt=least_squares(residual,[x,y],bounds=([x-1.4,y-1.4],[x+1.4,y+1.4]),
                      max_nfev=65,xtol=1e-6,ftol=1e-6)
    rr,coef=linfit(*opt.x)
    return dict(ok=bool(opt.success and coef[0]>0),
                x_stpsf=float(opt.x[0]),y_stpsf=float(opt.x[1]),
                x_gaussian=float(x),y_gaussian=float(y),
                dx_pix=float(opt.x[0]-x),dy_pix=float(opt.x[1]-y),
                chi2=float(np.dot(rr,rr)),n_valid=int(good.sum()),
                amplitude=float(coef[0]),at_bound=bool(np.max(np.abs(opt.x-[x,y]))>1.35))


def run(out):
    c=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    f=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    ra,dec=map(float,[c.loc[CID].ra_deg,c.loc[CID].dec_deg])
    pair=choose_frozen_pair(query_candidate_inventory(CID,ra,dec),EXPECTED_PAIRS[CID],
                            float(f.loc[CID].baseline_days))
    early,late=get_pair_products(pair)
    rows=[];coverage=[]
    for epoch,filter_name,products in [("early","F444W",early),("late","F356W",late)]:
        exposures,log=load_covering_cutouts(CID,ra,dec,epoch,filter_name,products)
        coverage.extend(log)
        print("FILTER",epoch,filter_name,"exposures",len(exposures),flush=True)
        if len(exposures)<1:raise RuntimeError("No CAL images for "+filter_name)
        for ex in exposures:
            original,controls=measure_exposure_recentered(ex)
            tx=ex.x_full-ex.x0;ty=ex.y_full-ex.y0
            candidates=control_candidates(ex,controls,tx,ty)
            print("FRAME",ex.filename,"eligible",len(candidates),"controls",len(controls),flush=True)
            if len(candidates)<3:
                rows.append(dict(epoch=epoch,filter=filter_name,filename=ex.filename,
                                 status="INSUFFICIENT_ISOLATED_COMPACT_CONTROLS",
                                 n_eligible=len(candidates)))
                continue
            # Use a fixed set of compact field controls across both spectra.
            for temp in PSF_TEMPS[filter_name]:
                psf,meta=get_psf(ex,temp)
                for i,(x,y) in enumerate(candidates):
                    try:
                        res=fit_star(ex,x,y,psf)
                        if res.get("ok"):
                            a=meas.local_pixel_to_sky(ex,x,y)
                            b=meas.local_pixel_to_sky(ex,res["x_stpsf"],res["y_stpsf"])
                            de,dn=meas.sky_to_tangent(b[0],b[1],a[0],a[1])
                            res.update(offset_east_mas=float(de*1000),
                                       offset_north_mas=float(dn*1000))
                        rows.append(dict(epoch=epoch,filter=filter_name,filename=ex.filename,
                                         temperature_hypothesis_K=temp,control_index=i,
                                         control_x=x,control_y=y,originally_accepted_target=bool(target_is_usable(original)),
                                         status="FIT" if res.get("ok") else "FAILED",
                                         **res))
                    except Exception as exc:
                        rows.append(dict(epoch=epoch,filter=filter_name,filename=ex.filename,
                                         temperature_hypothesis_K=temp,control_index=i,
                                         status="ERROR",error=type(exc).__name__+":"+str(exc)))
                pd.DataFrame(rows).to_csv(out/"all_control_fit_results.csv",index=False)
    pd.DataFrame(coverage).to_csv(out/"mast_coverage.csv",index=False)
    table=pd.DataFrame(rows)
    stats=[]
    for (epoch,flt,temp),g in table[table.status.eq("FIT")&(~table.at_bound.fillna(True))].groupby(
            ["epoch","filter","temperature_hypothesis_K"]):
        if len(g)<3:continue
        e=g.offset_east_mas.to_numpy(float);n=g.offset_north_mas.to_numpy(float)
        stats.append(dict(epoch=epoch,filter=flt,temp_hypothesis_K=int(temp),
                          n_controls=int(len(g)),n_frames=int(g.filename.nunique()),
                          median_east_mas=float(np.median(e)),
                          median_north_mas=float(np.median(n)),
                          scatter_east_mas=float(meas.robust_sigma(e)),
                          scatter_north_mas=float(meas.robust_sigma(n)),
                          p90_2d_offset_mas=float(np.percentile(np.hypot(e,n),90))))
    pd.DataFrame(stats).to_csv(out/"compact_control_psf_summary.csv",index=False)
    summary=dict(candidate_id=CID,n_fitted=int(table.status.eq("FIT").sum()),
                 stats=stats,
                 interpretation="Empirical compact-control STPSF minus 2D-Gaussian centroid scatter; NOT spectroscopically confirmed stars and NOT a cross-filter systematic by itself.",
                 stop_gate="If fewer than 3 isolated good controls in any filter or model scatter large, no new significance or PM claim; do not promote rejected F356W frames.",
                 no_change_to_paper=True)
    (out/"summary.json").write_text(json.dumps(safe_out(summary),indent=2))
    print("FINAL",json.dumps(safe_out(summary)),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,default=Path("results_1651050_psf_controls"))
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    install()
    try:run(a.out)
    except Exception:
        (a.out/"failure.txt").write_text(traceback.format_exc())
        raise
if __name__=="__main__":main()
