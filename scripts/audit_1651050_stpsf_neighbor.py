#!/usr/bin/env python3
"""1651050 bounded *physical* F356W PSF + bright-neighbor sensitivity pilot.

Stage-2 original MAST CAL data. STPSF/NIRCam on NRCB5 with 5700/900/500 K
hypothetical blackbody source spectra, synphot throughput weighting. The
5700 K case is a compact-star proxy, NOT the target SED; 500/900 K are ONLY
color-bracketing placeholders, NOT atmosphere models or fitted Teff.

A common F356W sky centroid is fit across the same four detector exposures
after same-filter local-control affine registration. Per exposure all source
amplitudes are UNCONSTRAINED signed weighted-LS parameters, plus background
constant or plane. Test a known bright local neighbor (position from accepted
F356W exposure) included/excluded and hold out original accepted detector.

IMPORTANT: This produces method-dependent centroid shifts and residual checks,
not a calibrated PSF correction, a new PM, new independent detections, or
a parallax measurement. Original rejection gates remain unchanged.
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, maximum_filter, map_coordinates
from scipy.optimize import least_squares
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"first26"),str(ROOT/"scripts")]
import pm86.measurement as meas
from pm86.archive import query_candidate_inventory,get_pair_products,load_covering_cutouts
from pm86.astrometry import _match_controls,_fit_affine
from run_cross_filter_null import choose_frozen_pair,EXPECTED_PAIRS,target_is_usable
from recenter_measurement import measure_exposure_recentered
from remote_io_patch import install
from audit_1651050_f356w_joint import make_frame, projected_matrix, safe_out
from audit_1651050_f356w_register import registered_projected_frame

CID=1651050
FOV=33
OVERSAMPLE=4
N_LAMBDA=7
TEMPERATURES=(500,900,5700)    # sensitivity scenarios; NOT Teff results
FITTING_HALF=15
MAX_OFFSET_MAS=85
MIN_REFS=6


def get_psf(ex,temperature):
    import stpsf
    from synphot import SourceSpectrum
    from synphot.models import BlackBodyNorm1D
    inst=stpsf.NIRCam()
    inst.filter="F356W"
    inst.detector=str(ex.detector).upper()   # NRCBLONG aliases NRCB5
    inst.detector_position=(int(round(ex.x_full)),int(round(ex.y_full)))
    source=SourceSpectrum(BlackBodyNorm1D,temperature=temperature*u.K)
    h=inst.calc_psf(source=source,nlambda=N_LAMBDA,
                    fov_pixels=FOV,oversample=OVERSAMPLE)
    if "OVERDIST" not in h:
        raise RuntimeError("STPSF OVERDIST oversampled detector PSF unavailable")
    data=np.asarray(h["OVERDIST"].data,float)
    if data.shape!=(FOV*OVERSAMPLE,FOV*OVERSAMPLE) or not np.isfinite(data).all():
        raise RuntimeError("Unexpected oversampled STPSF PSF shape or nonfinite pixels")
    if np.min(data)<-1e-10 or data.sum()<=0:
        raise RuntimeError("Unphysical STPSF oversampled PSF")
    data=np.maximum(data,0)
    data=data/data.sum()
    meta=dict(stpsf_version=stpsf.__version__,detector=inst.detector,
              detector_position=list(inst.detector_position),filter=inst.filter,
              nlambda=N_LAMBDA,source_temperature_kelvin=int(temperature),
              oversample=OVERSAMPLE,shape=list(data.shape),
              output_extension="OVERDIST",source_spectrum="synphot BlackBodyNorm1D")
    h.close()
    return data,meta


def native_template(image,xx,yy,x,y):
    """Interpolate oversampled DETECTOR-DISTORTED PSF at science pixel centers.

    This is a bounded pilot: interpolation/pixel integration needs independent
    pixel-phase calibration before mas-level systematics may be claimed.
    """
    size=OVERSAMPLE
    center=np.array(image.shape)/2.0-.5
    coords=np.array([center[0]+size*(yy-y),center[1]+size*(xx-x)])
    return np.asarray(map_coordinates(image,coords,order=3,
                                      mode="constant",cval=0.0),float)*size*size


def detect_neighbor(ref,accepted_x,accepted_y):
    """Select strongest isolated local brightness peak 7-14 detector pixels away."""
    sci=np.asarray(ref.sci,float)
    err=np.asarray(ref.err,float)
    dq=np.asarray(ref.dq,np.uint32)
    _,bad=meas.masks(sci,err,dq)
    yy,xx=np.indices(sci.shape)
    rad=np.hypot(xx-accepted_x,yy-accepted_y)
    ring=(rad>=15)&(rad<=20)&(~bad)&np.isfinite(sci)
    if ring.sum()<35:
        raise RuntimeError("Insufficient background pixels near target")
    background=float(np.nanmedian(sci[ring]))
    significance=np.zeros_like(sci)
    finite=(~bad)&np.isfinite(sci)&np.isfinite(err)&(err>0)
    significance[finite]=(sci[finite]-background)/err[finite]
    smooth=gaussian_filter(significance,1.0,mode="constant")
    candidates=(smooth==maximum_filter(smooth,size=5))
    candidates &= (rad>=7) & (rad<=14) & finite
    yy0,xx0=np.nonzero(candidates)
    peaks=sorted([(float(smooth[y,x]),float(x),float(y))
                  for y,x in zip(yy0,xx0)],reverse=True)
    if not peaks or peaks[0][0]<5.0:
        raise RuntimeError("No defensible >5-smooth-SNR local neighbor peak")
    significance_peak,x,y=peaks[0]
    return x,y,dict(peak_smoothed_snr=significance_peak,
                    neighbor_distance_pixels=float(np.hypot(x-accepted_x,y-accepted_y)),
                    local_background=background,
                    source_selection="strongest smoothed SNR local peak 7-14 px from accepted target")


def fit_fluxes(frame, psf, east, north, neighbor=None,
               background="plane"):
    x,y=frame["x0"]+frame["dx_de"]*east+frame["dx_dn"]*north,\
        frame["y0"]+frame["dy_de"]*east+frame["dy_dn"]*north
    cols=[native_template(psf,frame["xx"],frame["yy"],x,y)]
    if neighbor is not None:
        xn=frame["x0"]+frame["dx_de"]*neighbor[0]+frame["dx_dn"]*neighbor[1]
        yn=frame["y0"]+frame["dy_de"]*neighbor[0]+frame["dy_dn"]*neighbor[1]
        cols.append(native_template(psf,frame["xx"],frame["yy"],xn,yn))
    cols.append(np.ones(len(frame["xx"])))
    if background=="plane":
        cols.extend([frame["xx"]-frame["x0"],frame["yy"]-frame["y0"]])
    A=np.column_stack(cols)
    w=1/np.asarray(frame["err"],float)
    coeff, *_=np.linalg.lstsq(A*w[:,None],frame["sci"]*w,rcond=1e-8)
    resid=(frame["sci"]-A@coeff)*w
    stat=dict(filename=frame["filename"],target_amplitude=float(coeff[0]),
              neighbor_amplitude=float(coeff[1]) if neighbor is not None else None,
              background=float(coeff[2]) if neighbor is not None else float(coeff[1]),
              residual_chi2=float(np.dot(resid,resid)),n_pixels=len(resid),
              target_flux_sign=("positive" if coeff[0]>0 else "nonpositive"))
    return resid,stat


def model_fit(frames, psfs, neighbor, fit_neighbor, background,
              without_accepted, accepted_filename):
    subset=[f for f in frames if not(without_accepted and f["filename"]==accepted_filename)]
    def residual(v):
        pieces=[]
        for f in subset:
            rr,_=fit_fluxes(f,psfs[f["filename"]],v[0],v[1],
                            neighbor if fit_neighbor else None,background)
            pieces.append(rr)
        return np.concatenate(pieces)
    result=least_squares(residual,[0.,0.],
                         bounds=([-MAX_OFFSET_MAS]*2,[MAX_OFFSET_MAS]*2),
                         max_nfev=45,xtol=1e-5,ftol=1e-5,gtol=1e-5)
    per_frame=[fit_fluxes(f,psfs[f["filename"]],
                          float(result.x[0]),float(result.x[1]),
                          neighbor if fit_neighbor else None,background)[1]
               for f in subset]
    return dict(success=bool(result.success),
                centroid_offset_east_mas=float(result.x[0]),
                centroid_offset_north_mas=float(result.x[1]),
                centroid_offset_amplitude_mas=float(np.linalg.norm(result.x)),
                is_position_bound=bool(np.max(np.abs(result.x))>.97*MAX_OFFSET_MAS),
                chi2=float(np.dot(result.fun,result.fun)),
                n_frames=len(subset),
                excluded_accepted_frame=bool(without_accepted),
                neighbor_mode="two_source" if fit_neighbor else "target_only",
                background_mode=background,
                per_frame=per_frame,
                n_negative_target_amplitude=int(sum(x["target_amplitude"]<=0 for x in per_frame)))


def gaussian_selftest():
    shape=(31,31)
    yy,xx=np.indices(shape)
    psf=np.exp(-.5*((yy-15)**2+(xx-15)**2)/1.4**2)
    oversamp=np.kron(psf,np.ones((OVERSAMPLE,OVERSAMPLE)))
    oversamp=oversamp/oversamp.sum()
    f=dict(filename="test",xx=xx.ravel(),yy=yy.ravel(),sci=(
        9+40*np.exp(-.5*((yy-15.13)**2+(xx-14.87)**2)/1.4**2)).ravel(),
        err=np.ones(shape).ravel(),x0=15,y0=15,
        dx_de=.015,dx_dn=.001,dy_de=.001,dy_dn=.015)
    r=model_fit([f],{"test":oversamp},neighbor=None,fit_neighbor=False,
                background="constant",without_accepted=False,accepted_filename="none")
    if not r["success"] or r["centroid_offset_amplitude_mas"]>50:
        raise RuntimeError("Oversampled PSF sampler / joint solver selftest failure:"+str(r))
    print("SYNTHETIC_SELFTEST_SUCCESS",json.dumps({k:v for k,v in r.items() if k!="per_frame"}))


def run(out):
    catalog=pd.read_csv(ROOT/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    frozen=pd.read_csv(ROOT/"data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    ra,dec=map(float,[catalog.loc[CID].ra_deg,catalog.loc[CID].dec_deg])
    pair=choose_frozen_pair(query_candidate_inventory(CID,ra,dec),
                            EXPECTED_PAIRS[CID],float(frozen.loc[CID].baseline_days))
    _,lp=get_pair_products(pair)
    exposures,coverage=load_covering_cutouts(CID,ra,dec,"late","F356W",lp)
    pd.DataFrame(coverage).to_csv(out/"mast_coverage.csv",index=False)
    if len(exposures)!=4:
        raise RuntimeError("Four original F356W CAL frames required; got "+str(len(exposures)))
    original,controls={},{}
    for ex in exposures:
        r,c=measure_exposure_recentered(ex)
        original[ex.filename]=r
        controls[ex.filename]=c
    accepted=[ex for ex in exposures if target_is_usable(original[ex.filename])]
    if len(accepted)!=1:
        raise RuntimeError("Unexpected accepted original F356W count: "+str(len(accepted)))
    ref=accepted[0]
    row=original[ref.filename]
    refsky=meas.local_pixel_to_sky(ref,row["x_2dg"],row["y_2dg"])
    x_neigh,y_neigh,neighbor_meta=detect_neighbor(ref,row["x_2dg"],row["y_2dg"])
    ns=meas.local_pixel_to_sky(ref,x_neigh,y_neigh)
    ne,nn=meas.sky_to_tangent(ns[0],ns[1],refsky[0],refsky[1])
    neighbor=(float(ne)*1000,float(nn)*1000)
    neighbor_meta["east_north_from_target_mas"]=list(neighbor)
    neighbor_meta["accepted_frame"]=ref.filename
    (out/"neighbor_selection.json").write_text(json.dumps(safe_out(neighbor_meta),indent=2))
    print("NEIGHBOR",json.dumps(safe_out(neighbor_meta)),flush=True)

    reference_controls=controls[ref.filename]
    if len(reference_controls)<MIN_REFS:raise RuntimeError("Insufficient reference control sources")
    frames=[]
    reg_info=[]
    for ex in exposures:
        if ex.filename==ref.filename:
            coeff=np.array([[0.,0.],[1.,0.],[0.,1.]])
            n_used=len(reference_controls)
            scat=[0.,0.]
        else:
            cidx,ridx=_match_controls(controls[ex.filename],reference_controls)
            a=controls[ex.filename].iloc[cidx][["east_arcsec_raw_wcs","north_arcsec_raw_wcs"]].to_numpy(float)
            b=reference_controls.iloc[ridx][["east_arcsec_raw_wcs","north_arcsec_raw_wcs"]].to_numpy(float)
            coeff,good,res,sx,sy=_fit_affine(a,b)
            n_used=int(good.sum());scat=[float(sx)*1000,float(sy)*1000]
        frame,_,_=registered_projected_frame(ex,(ra,dec),refsky,coeff)
        # Redo stamp with generous spatial window to INCLUDE bright neighbor
        # near ~11 detector pixels from target.
        ff=make_frame(ex.sci,ex.err,ex.dq,frame["x0"],frame["y0"],
                      frame["dx_de"],frame["dx_dn"],frame["dy_de"],frame["dy_dn"],
                      filename=ex.filename,h=FITTING_HALF)
        frames.append(ff)
        reg_info.append(dict(filename=ex.filename,originally_accepted=bool(target_is_usable(original[ex.filename])),
                             n_field_controls=int(len(controls[ex.filename])),
                             n_matched_controls=int(n_used),
                             reference_registration_scatter_mas=scat,
                             masked_fraction=ff["dq_stamp_fraction"],
                             center_dq=int(ff["center_dq"])))
    (out/"registration_and_masks.json").write_text(json.dumps(safe_out(reg_info),indent=2))

    all_results=[]
    model_meta=[]
    for temperature in TEMPERATURES:
        psfs={}
        for ex in exposures:
            model,meta=get_psf(ex,int(temperature))
            psfs[ex.filename]=model
            meta["filename"]=ex.filename
            model_meta.append(meta)
        for background,with_neighbor in [("constant",False),("plane",False),
                                           ("constant",True),("plane",True)]:
            for without in [False,True]:
                r=model_fit(frames,psfs,neighbor,with_neighbor,
                            background,without,ref.filename)
                r["hypothetical_temperature_kelvin"]=int(temperature)
                all_results.append(r)
                print("STPSF_RESULT",json.dumps(safe_out({k:v for k,v in r.items() if k!="per_frame"})),flush=True)

        # checkpoint as soon as each spectral scenario is complete
        (out/"fit_results.json").write_text(json.dumps(safe_out(all_results),indent=2))
        (out/"model_provenance.json").write_text(json.dumps(safe_out(model_meta),indent=2))

    records=[]
    for r in all_results:
        records.append({k:v for k,v in r.items() if k!="per_frame"})
    pd.DataFrame(records).to_csv(out/"physical_psf_model_comparisons.csv",index=False)
    signs=[]
    for r in all_results:
        for pr in r["per_frame"]:
            signs.append({**pr,**{k:v for k,v in r.items() if k not in ("per_frame","chi2")}})
    pd.DataFrame(signs).to_csv(out/"source_amplitudes_by_frame.csv",index=False)
    summary=dict(candidate_id=CID,workflow="STPSF F356W same-filter controlled PSF and neighbor pilot",
                 input_images=4,n_accepted_original=1,
                 neighbor=neighbor_meta,
                 spectral_scenarios_kelvin=list(TEMPERATURES),
                 model_repetitions=len(all_results),
                 fit_runs=[{k:v for k,v in r.items() if k!="per_frame"} for r in all_results],
                 registration=reg_info,
                 scientific_status="SENSITIVITY_ONLY_NOT_A_NEW_PM_OR_VALIDATED_ASTROMETRY",
                 caveats=[
                  "Blackbody source assumptions only bracket within-filter PSF color; real ultra-cool dwarf has molecular absorption and spectrum is not supplied.",
                  "Interpolation of oversampled STPSF pixels into CAL detector pixels is not independently calibrated at milli-arcsecond scale.",
                  "Detector optical path difference time variability and PSF validation on truly identified stars not established.",
                  "The bright neighbor is selected from same accepted F356W frame and held fixed in the registered position; neighbor flux may vary by detector.",
                  "A planar sky is not a complete neighbor/deblending model and 31x31 stamp includes other weak sources.",
                  "Every reported centroid is tied to the originally accepted F356W sky position; rejected frames remain rejected.",
                  "No new injection calibration or 5-sigma claims, original frozen PM and parallax analyses untouched."
                 ])
    (out/"summary.json").write_text(json.dumps(safe_out(summary),indent=2))
    print("FINAL_AUDIT_SUMMARY",json.dumps(safe_out({k:v for k,v in summary.items() if k!="fit_runs"})),flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("results_1651050_stpsf_pilot"))
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        gaussian_selftest()
        return
    a.out.mkdir(parents=True,exist_ok=True)
    install()
    try:
        run(a.out)
    except Exception:
        (a.out/"failure.txt").write_text(traceback.format_exc())
        raise


if __name__=="__main__":
    main()
