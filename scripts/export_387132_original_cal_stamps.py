#!/usr/bin/env python3
"""Portable source-image package for COSMOS-Web candidate 387132 (legacy 282040).

Reads six individual original JWST/NIRCam Stage-2 CAL products by remote
byte range; saves sky-anchored local SCI ERR DQ stamps with provenance.
No new independent PM measurement or parallax fit is attempted.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from pm86.archive import load_covering_cutouts

CID, RA, DEC=387132,150.1513055210011,1.9720770310321503

def process(row,outdir):
    filename=str(row.productFilename); filt=str(row["filter"])
    products=pd.DataFrame([{"productFilename":filename,
                            "dataURI":"mast:JWST/product/"+filename,
                            "productType":"SCIENCE"}])
    c,rcpt=load_covering_cutouts(CID,RA,DEC,str(row.epoch),filt,products)
    if len(c)!=1:raise RuntimeError(str(rcpt))
    e=c[0]
    tx=float(e.x_full-e.x0);ty=float(e.y_full-e.y0)
    x=int(round(tx));y=int(round(ty));h=36
    x0=max(0,x-h);y0=max(0,y-h);x1=min(e.sci.shape[1],x+h+1);y1=min(e.sci.shape[0],y+h+1)
    if min(x1-x0,y1-y0)<55:raise RuntimeError("Source too near the cropped CAL boundary")
    img=np.asarray(e.sci[y0:y1,x0:x1],dtype=np.float32)
    err=np.asarray(e.err[y0:y1,x0:x1],dtype=np.float32)
    dq=np.asarray(e.dq[y0:y1,x0:x1],dtype=np.uint32)
    stem=filename[:-9] if filename.endswith("_cal.fits") else filename
    path=outdir/"pixel_stamps"/(stem+".npz")
    np.savez_compressed(path,SCI=img,ERR=err,DQ=dq,
        pixel_window_x0=int(x0),pixel_window_y0=int(y0),
        full_detector_window_x0=int(e.x0+x0),
        full_detector_window_y0=int(e.y0+y0),
        catalog_x_local=float(tx-x0),catalog_y_local=float(ty-y0),
        filter=str(filt),epoch=str(row.epoch),
        filename=filename,mjd=float(e.mjd),
        ra_deg=float(RA),dec_deg=float(DEC),
        detector=str(e.detector),
        note="Original calibrated detector pixel stamp, not a registered image.")
    return dict(filename=filename,filter=filt,epoch=str(row.epoch),
                mjd=float(e.mjd),stamp_name=path.name,
                sci_shape=str(img.shape),catalog_x_in_stamp=tx-x0,
                catalog_y_in_stamp=ty-y0,
                dq_center=int(dq[int(round(ty-y0)),int(round(tx-x0))]),
                pixel_origin_full_x=int(e.x0+x0),pixel_origin_full_y=int(e.y0+y0))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",default="publication_check/387132_original_cal_stamps")
    a=ap.parse_args()
    root=Path(a.output)
    (root/"pixel_stamps").mkdir(parents=True,exist_ok=True)
    names=[
    ("F444W","early","jw01727139001_04101_00001_nrcblong_cal.fits"),
    ("F444W","early","jw01727139001_04101_00002_nrcblong_cal.fits"),
    ("F444W","early","jw01727139001_04101_00003_nrcblong_cal.fits"),
    ("F444W","early","jw01727139001_04101_00004_nrcblong_cal.fits"),
    ("F410M","late","jw06434296001_06201_00001_nrcalong_cal.fits"),
    ("F410M","late","jw06434296001_06201_00002_nrcalong_cal.fits")]
    results=[];errors=[]
    for filt,epoch,filename in names:
        print(f"Retrieving {filt} {epoch} original CAL {filename}",flush=True)
        try:results.append(process(pd.Series(dict(filter=filt,epoch=epoch,
                                         productFilename=filename)),root))
        except Exception as ex:
            msg=f"{filename}: {type(ex).__name__}: {ex}"
            print(msg,flush=True)
            errors.append(msg)
    pd.DataFrame(results).to_csv(root/"ORIGINAL_EXPOSURE_MANIFEST.csv",index=False)
    (root/"README.md").write_text(
       "# 387132 original JWST detector stamps\n"
       "These six CAL cutouts use embedded JWST distortion-aware gWCS for target coverage. "
       "Each image is an individual detector stamp and its sky position is defined "
       "only by the fixed catalog coordinate. They have NOT been registered together. "
       "Published proper motion should be read from the independent PM audit, "
       "not inferred from differences between displayed stamp centers.\n"
       "For accurate PM, propagate registration errors and JWST orbital parallax.\n")
    if errors:(root/"ERRORS.txt").write_text("\n".join(errors))
    if len(results)!=6:raise RuntimeError(f"Only {len(results)}/6 stamps")
    print("All six original CAL stamps generated",flush=True)

if __name__=="__main__":main()
