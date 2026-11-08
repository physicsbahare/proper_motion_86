#!/usr/bin/env python3
"""1651050-only F356W rejected-frame / joint-source-location feasibility audit.

PUBLIC JWST CAL, original masks and per-detector WCS. This tests a strictly
bounded, *fixed-width analytic Gaussian* with common sky location in four
F356W detector images. Each image has its own nonnegative source amplitude and
constant local background. Frames rejected by the original pipeline REMAIN
rejected: this is an independent feasibility/robustness diagnostic, NOT a
replacement paper PM, ePSF, chromatic PSF, detection sigma, or parallax.

Outputs preserve every rejected frame, pixel DQ and fixed-location flux test,
including leave-one-frame-out and accepted-frame-excluded solutions.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "first26"), str(ROOT / "scripts")]
import pm86.measurement as measure
from pm86.archive import query_candidate_inventory, get_pair_products, load_covering_cutouts
from pm86.config import CONFIG
from run_cross_filter_null import choose_frozen_pair, EXPECTED_PAIRS, target_is_usable
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install


CID = 1651050
WIDTHS = (0.8, 1.05, 1.3)
WINDOW_HALF = 7
MAX_SHIFT_MAS = 160.0


def finite_float(v):
    try:
        z = float(v)
        return z if np.isfinite(z) else None
    except Exception:
        return None


def make_frame(sci, err, dq, x0, y0, dx_de, dx_dn, dy_de, dy_dn,
               filename="synthetic", *, h=WINDOW_HALF):
    """Build a fixed-grid SCI/ERR/DQ stamp anchored on a specified sky location."""
    sci = np.asarray(sci, dtype=float)
    err = np.asarray(err, dtype=float)
    dq = np.asarray(dq, dtype=np.uint32)
    bad = (~np.isfinite(sci) | ~np.isfinite(err) | (err <= 0)
           | ((dq & np.uint32(CONFIG.astrom_bad_bits)) != 0))
    xc, yc = int(round(x0)), int(round(y0))
    xlo, xhi = xc - h, xc + h + 1
    ylo, yhi = yc - h, yc + h + 1
    if min(xlo, ylo) < 0 or xhi > sci.shape[1] or yhi > sci.shape[0]:
        raise ValueError("Detector cutout does not contain full fitting stamp")
    yy, xx = np.mgrid[ylo:yhi, xlo:xhi]
    cut = (slice(ylo, yhi), slice(xlo, xhi))
    ok = ~bad[cut]
    if np.count_nonzero(ok) < 30:
        raise ValueError("Fewer than 30 unmasked SCI/ERR pixels")
    return dict(filename=filename, xx=xx[ok].astype(float), yy=yy[ok].astype(float),
                sci=sci[cut][ok], err=err[cut][ok],
                x0=float(x0), y0=float(y0),
                dx_de=float(dx_de), dx_dn=float(dx_dn),
                dy_de=float(dy_de), dy_dn=float(dy_dn),
                n_valid=int(np.count_nonzero(ok)),
                n_masked=int(ok.size - np.count_nonzero(ok)),
                center_dq=int(dq[yc, xc]),
                dq_core_fraction=float(np.count_nonzero(bad[
                    max(0,yc-1):yc+2,max(0,xc-1):xc+2])/9.0),
                dq_stamp_fraction=float(1.0 - np.mean(ok)),
                image_stamp=sci[cut], err_stamp=err[cut], dq_stamp=dq[cut],
                stamp_xlo=xlo, stamp_ylo=ylo)


def proj(frame, east_mas, north_mas):
    x = frame["x0"] + frame["dx_de"]*east_mas + frame["dx_dn"]*north_mas
    y = frame["y0"] + frame["dy_de"]*east_mas + frame["dy_dn"]*north_mas
    return x, y


def best_linear_flux(frame, e, n, width):
    """Profile out image amplitude >=0 and constant local background."""
    x,y = proj(frame, e, n)
    psf = np.exp(-.5*((frame["xx"]-x)**2+(frame["yy"]-y)**2)/width**2)
    weights = 1/frame["err"]
    A = np.column_stack((psf*weights, weights))
    coeff = np.linalg.lstsq(A, frame["sci"]*weights, rcond=None)[0]
    amp = float(coeff[0])
    if amp < 0:
        amp = 0.0
        bkg = float(np.sum(frame["sci"]*weights**2)/np.sum(weights**2))
    else:
        bkg = float(coeff[1])
    resid = (frame["sci"]-amp*psf-bkg)*weights
    null_bkg = float(np.sum(frame["sci"]*weights**2)/np.sum(weights**2))
    null_chisq = float(np.sum(((frame["sci"]-null_bkg)*weights)**2))
    # This is a conditional amplitude-to-error ratio at fixed center/width.
    # Do NOT label it a per-frame independent source-detection significance.
    S = np.array([[np.sum((psf*weights)**2), np.sum(psf*weights**2)],
                  [np.sum(psf*weights**2), np.sum(weights**2)]])
    try:
        covariance = np.linalg.pinv(S)
        amplitude_err = math.sqrt(max(0.0, float(covariance[0,0])))
    except Exception:
        amplitude_err = np.nan
    return resid, dict(amplitude=amp, background=bkg, conditional_amp_snr=(
        amp/amplitude_err if np.isfinite(amplitude_err) and amplitude_err>0 else np.nan),
        chi2=float(np.sum(resid**2)), null_chi2=null_chisq,
        source_center_x=x,source_center_y=y)


def fit_shared(frames,width,shift_limit=MAX_SHIFT_MAS):
    if not frames:
        raise ValueError("No frames passed the common WCS fit cut")
    def residual(v):
        return np.concatenate([best_linear_flux(f,float(v[0]),float(v[1]),width)[0]
                               for f in frames])
    opt=least_squares(residual,x0=[0.,0.],bounds=([-shift_limit,-shift_limit],
                      [shift_limit,shift_limit]),max_nfev=140,
                      xtol=1e-8,ftol=1e-8,gtol=1e-8)
    e,n=map(float,opt.x)
    rows=[]
    for f in frames:
        _,stat=best_linear_flux(f,e,n,width)
        stat["filename"]=f["filename"]
        stat["n_valid"]=f["n_valid"]
        rows.append(stat)
    return dict(success=bool(opt.success),offset_east_mas=e,offset_north_mas=n,
                shift_magnitude_mas=float(np.hypot(e,n)),
                at_position_bound=bool(max(abs(e),abs(n))>0.98*shift_limit),
                width_pixels=width,total_chi2=float(np.sum(opt.fun**2)),
                nfev=int(opt.nfev),n_frames=len(frames),per_frame=rows)


def projected_matrix(ex,ref_ra,ref_dec):
    anchor=SkyCoord(ref_ra*u.deg,ref_dec*u.deg)
    def px(pos):
        q=ex.gwcs.invert(float(pos.ra.deg),float(pos.dec.deg))
        return np.asarray(q,dtype=float)-[ex.x0,ex.y0]
    center=px(anchor)
    eps=20.0
    ee=px(anchor.spherical_offsets_by(eps*u.mas,0*u.mas))
    nn=px(anchor.spherical_offsets_by(0*u.mas,eps*u.mas))
    return center, (ee-center)/eps, (nn-center)/eps


def frame_diagnostic_plot(frames,original_records,out):
    if not frames: return
    fig,axs=plt.subplots(len(frames),3,figsize=(11,3.1*len(frames)),squeeze=False)
    for i,f in enumerate(frames):
        s=f["image_stamp"]
        good=np.isfinite(s)
        vmin,vmax=np.percentile(s[good],[5,99]) if np.any(good) else (0,1)
        for j,(data,name) in enumerate([(s,"SCI"),(f["err_stamp"],"ERR"),
                                          (f["dq_stamp"],"DQ")]):
            ax=axs[i,j]
            ax.imshow(data,origin="lower",cmap="gray" if j<2 else "viridis",
                      vmin=vmin if j==0 else None,vmax=vmax if j==0 else None,
                      interpolation="nearest")
            cx=f["x0"]-f["stamp_xlo"]
            cy=f["y0"]-f["stamp_ylo"]
            ax.plot(cx,cy,marker="+",color="red",markersize=9)
            ax.set_title((f["filename"].split("_")[-2] if i==0 else "")+" "+name,fontsize=8)
            ax.set_xticks([]);ax.set_yticks([])
            if j==0:
                ax.set_ylabel(f["filename"].replace("_cal.fits","")[-26:],fontsize=7)
    fig.suptitle("1651050 F356W original detector stamps; cross = accepted-frame sky location",fontsize=11)
    fig.tight_layout()
    fig.savefig(out,dpi=145,bbox_inches="tight")
    plt.close(fig)


def safe_out(x):
    if isinstance(x,dict): return {str(k):safe_out(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [safe_out(v) for v in x]
    if isinstance(x,(np.integer,)):return int(x)
    if isinstance(x,(np.floating,)):return safe_out(float(x))
    if isinstance(x,float) and not np.isfinite(x):return None
    return x


def self_test():
    # Unit-level recovery of a common true detector shift with nonidentical
    # pixel phases, missing pixels and different frame amplitudes.
    rng=np.random.default_rng(1651050)
    fs=[]
    truth=np.array([23.,-16.])
    for i in range(4):
        shape=(41,41)
        x0=19.9+0.16*i
        y0=20.1-0.09*i
        m=np.array([[.013+.0001*i,.003],[-.002,.014]])
        xx,yy=np.indices(shape)[::-1]
        x=x0+np.dot(m[0],truth)
        y=y0+np.dot(m[1],truth)
        sci=10.+(9.+2*i)*np.exp(-0.5*((xx-x)**2+(yy-y)**2)/1.05**2)
        sci+=rng.normal(0,.45,size=shape)
        err=np.full(shape,.45)
        dq=np.zeros(shape,dtype=np.uint32)
        dq[18+i,17]=4
        fs.append(make_frame(sci,err,dq,x0,y0,*m[0],*m[1],filename=str(i)))
    f=fit_shared(fs,1.05)
    assert f["success"] and np.linalg.norm(
        np.asarray([f["offset_east_mas"],f["offset_north_mas"]])-truth
    )<15, f
    print("SELF_TEST_SUCCESS",json.dumps({"fitted_mas":[f["offset_east_mas"],
                                                        f["offset_north_mas"]]}))


def run(out):
    cid=CID
    catalog=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    frozen=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    ra,dec=map(float,[catalog.loc[cid].ra_deg,catalog.loc[cid].dec_deg])
    original=frozen.loc[cid]
    pair=choose_frozen_pair(query_candidate_inventory(cid,ra,dec),
                            EXPECTED_PAIRS[cid],float(original.baseline_days))
    _,late_products=get_pair_products(pair)
    late,audit=load_covering_cutouts(cid,ra,dec,"late","F356W",late_products)
    pd.DataFrame(audit).to_csv(out/"mast_f356w_coverage.csv",index=False)
    if not late: raise RuntimeError("No F356W CAL frames available")

    records,expositions={},{}
    for e in late:
        row,controls=measure_exposure_recentered(e)
        records[e.filename]=row
        expositions[e.filename]=e
    accepted=[x for x in late if target_is_usable(records[x.filename])]
    if len(accepted)!=1:
        raise RuntimeError(f"Expected exactly one baseline-accepted F356W frame; saw {len(accepted)}")
    ref=accepted[0]
    row=records[ref.filename]
    if not np.isfinite([row["x_2dg"],row["y_2dg"]]).all():
        raise RuntimeError("Reference accepted forced Gaussian missing centroid")
    ref_ra,ref_dec=measure.local_pixel_to_sky(ref,row["x_2dg"],row["y_2dg"])
    original_rows=[]
    for e in late:
        r=records[e.filename]
        d={k:safe_out(v) for k,v in r.items() if k in (
            "filename","filter","epoch_label","astrometric_snr","clean_snr_err",
            "raw_snr_err","astrometry_source","forced_fit_attempted",
            "forced_fit_accepted","forced_snr","forced_x","forced_y",
            "forced_x_err_pix","forced_y_err_pix","forced_axis_ratio",
            "forced_sigma_x_pix","forced_sigma_y_pix","forced_reduced_chi2",
            "forced_npix","dq_or_r4","dq_center","dq_jump_det","x_2dg","y_2dg",
            "x_com","y_com","aperture_clean_snr_err")}
        d.update(filename=e.filename,
                 accepted_original=bool(target_is_usable(r)),
                 same_sky_source_verified=False)
        original_rows.append(d)
    pd.DataFrame(original_rows).to_csv(out/"original_f356w_gate_and_forced_diagnostics.csv",index=False)

    frames=[]
    for e in late:
        try:
            center,de,dn=projected_matrix(e,ref_ra,ref_dec)
            f=make_frame(e.sci,e.err,e.dq,*center,
                         de[0],dn[0],de[1],dn[1],filename=e.filename)
            frames.append(f)
        except Exception as exc:
            print("FRAME_NOT_FITTABLE",e.filename,type(exc).__name__,str(exc),flush=True)
    frame_diagnostic_plot(frames,records,out/"four_f356w_original_stamps.png")
    if len(frames)!=len(late):
        raise RuntimeError("Not all four frames can be included; see output diagnostics")

    full=[]
    for width in WIDTHS:
        for leave in [None,*[f["filename"] for f in frames]]:
            selected=[f for f in frames if f["filename"]!=leave]
            if len(selected)==0:continue
            res=fit_shared(selected,width)
            res["left_out_frame"]=leave
            res["included_accepted_frame"]=any(
                f["filename"]==ref.filename for f in selected)
            full.append(res)
            print("JOINT_RESULT",json.dumps(safe_out({
                "width":width,"left_out":leave,
                "offset_east":res["offset_east_mas"],
                "offset_north":res["offset_north_mas"],
                "on_bound":res["at_position_bound"]})),flush=True)
    (out/"joint_diagnostic_full.json").write_text(json.dumps(safe_out(full),indent=2))
    overview=pd.DataFrame([{k:v for k,v in r.items() if k!="per_frame"} for r in full])
    overview.to_csv(out/"joint_loo_width_summary.csv",index=False)
    fluxes=[]
    for res in full:
        for f in res["per_frame"]:
            fluxes.append({**f,"width_pixels":res["width_pixels"],
                           "left_out_frame":res["left_out_frame"]})
    pd.DataFrame(fluxes).to_csv(out/"conditional_per_frame_amplitudes.csv",index=False)
    summary={
        "candidate_id":cid,
        "n_f356w_covered":len(late),
        "n_f356w_original_accepted":len(accepted),
        "accepted_frame":ref.filename,
        "initial_center_provenance":"original forced_gaussian sky position in accepted F356W CAL",
        "model":"shared F356W WCS sky centroid, fixed circular Gaussian width per sensitivity run; independent nonnegative amp and constant background each frame",
        "analyses":"full 4-frame vs leave-one-frame-out, including accepted-frame-excluded solution, width 0.8/1.05/1.3 detector pixel",
        "n_frames_in_joint":len(frames),
        "full_joint_by_width":[{k:v for k,v in r.items() if k!="per_frame"}
                               for r in full if r["left_out_frame"] is None],
        "without_accepted_by_width":[{k:v for k,v in r.items() if k!="per_frame"}
                                     for r in full if r["left_out_frame"]==ref.filename],
        "science_status":"FEASIBILITY_ONLY_NOT_RESCUED_NOT_INDEPENDENT_PM",
        "cautions":[
            "No color-dependent JWST physical PSF, PSF spatial variation or extended neighbor model.",
            "The seed is the sole originally accepted F356W centroid; an apparent joint detection could inherit that seed.",
            "Frames rejected by original pipeline have NOT been admitted into published astrometry.",
            "Nonnegative amplitudes at fixed fitted sky position are selection-conditioned; their SNRs are not independent detection significances.",
            "No joint centroid covariance, calibrated WCS systematics, calibrated injection false-positive rate, or parallax.",
            "If excluding accepted frame shifts the fit to a boundary or across PSF widths, the rejected data cannot independently confirm the sky location."
        ]
    }
    (out/"summary.json").write_text(json.dumps(safe_out(summary),indent=2))
    print("AUDIT_SUMMARY",json.dumps(safe_out(summary)),flush=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,default=Path("results_1651050_f356w_joint"))
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test();return
    args.out.mkdir(parents=True,exist_ok=True)
    install()
    try:run(args.out)
    except Exception:
        (args.out/"failure.txt").write_text(traceback.format_exc())
        raise


if __name__=="__main__":main()
