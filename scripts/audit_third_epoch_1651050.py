#!/usr/bin/env python3
"""Up-to-date public MAST coverage audit for third astrometric epoch of 1651050.

Crucial distinction:
- MAST query returns *candidate observation-level* footprint intersections.
- Product-level original CAL WCS checking determines actual sky coverage.
- A covered detector DOES NOT imply the target is detected.
- 3+ nominal visits DOES NOT imply three independent astrometric positions.

Download original product metadata and WCS only; never automatically upgrade a
proper motion or count optical/blue non-detections as a new UCD position.
"""
from __future__ import annotations

import argparse, copy, json, re, traceback
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astroquery.mast import Observations
from stdatamodels import asdf_in_fits

from pm86.archive import mast_download_url, _open_remote_cal

CID=1651050
RA=150.307303724545
DEC=2.27642841317716
RELEVANT = {"F444W","F410M","F430M","F356W","F335M","F322W","F277W",
            "F770W","F1000W","F1280W","F1500W","F1800W","F2100W",
            "F160W"}
def normalize_instrument_filter(value):
    """JWST MIRI uses four-digit filters (e.g. F1000W, F2100W)."""
    match=re.search(r"F\d{3,4}[WMN]",str(value).upper())
    return match.group(0) if match else None


RE_FITS=re.compile(r"\.(fits|fit)$",re.IGNORECASE)


def product_type(filename, instrument):
    name=str(filename).lower()
    inst=str(instrument).upper()
    if not RE_FITS.search(name):
        return None
    if "JWST" in inst:
        # Stage-2 image calibrations only. Avoid association mosaics that combine epochs.
        return "JWST_CAL" if name.endswith("_cal.fits") else None
    if "HST" in inst:
        # Separate native exposures prefer flt/flc over drizzled combined visits.
        if name.endswith(("_flt.fits","_flc.fits")):
            return "HST_FLT_FLC"
        if name.endswith(("_drz.fits","_drc.fits")):
            return "HST_DRZ_DRC"
    return None


def list_inventory(radius_arcsec=3.0):
    position=SkyCoord(RA*u.deg,DEC*u.deg)
    tab=Observations.query_region(position,radius=radius_arcsec*u.arcsec)
    if not len(tab):
        raise RuntimeError("No MAST observations at sky coordinate")
    obs=tab.to_pandas()
    obs=obs[obs.obs_collection.astype(str).str.upper().isin(["JWST","HST"])].copy()
    if "dataproduct_type" in obs:
        obs=obs[obs.dataproduct_type.astype(str).str.lower().eq("image")]
    if "dataRights" in obs:
        rights=obs.dataRights.astype(str).str.upper()
        obs=obs[rights.isin(["PUBLIC","","NAN","NONE"])].copy()
    obs["filter_norm"]=obs.filters.map(normalize_instrument_filter)
    obs["start_mjd"]=pd.to_numeric(obs.t_min,errors="coerce")
    obs["end_mjd"]=pd.to_numeric(obs.t_max,errors="coerce")
    obs["red_relevant"]=obs.filter_norm.isin(RELEVANT)
    return obs.sort_values("start_mjd")


def all_products(obs):
    ids=obs["obsid"].dropna().astype("int64").astype(str).drop_duplicates().tolist()
    frames=[]
    # Group to avoid overwhelming MAST with one multi-hundred-observation request
    for start in range(0,len(ids),15):
        table=Observations.get_product_list(ids[start:start+15])
        if len(table):
            frames.append(table.to_pandas())
    if not frames:
        return pd.DataFrame()
    products=pd.concat(frames,ignore_index=True).drop_duplicates(subset="dataURI")
    filt=products.productFilename.astype(str).str.contains(RE_FITS,regex=True,na=False)
    products=products.loc[filt].copy()
    return products


def check_jwst_wcs(filename,data_uri):
    """Read JWST FITS headers and ASDF-in-FITS WCS, not the SCI pixel array."""
    with _open_remote_cal(str(data_uri)) as hdul:
        if "SCI" not in hdul:
            return {"status":"NO_SCI"}
        ny,nx=hdul["SCI"].shape
        with asdf_in_fits.open(hdul) as af:
            g=copy.deepcopy(af.tree["meta"]["wcs"])
            pix=g.invert(float(RA),float(DEC))
        px=float(np.asarray(pix[0]).squeeze())
        py=float(np.asarray(pix[1]).squeeze())
        in_detector=bool(np.isfinite([px,py]).all() and 0<=px<nx and 0<=py<ny)
        hdr=hdul[0].header
        return {"status":"COVERS_TARGET" if in_detector else "OUTSIDE_DETECTOR",
                "detector":str(hdr.get("DETECTOR","")),
                "mjd":float(hdr.get("EXPSTART",np.nan)),
                "x_full":px,"y_full":py,"image_nx":nx,"image_ny":ny}


def check_hst_wcs(filename,data_uri):
    with fits.open(mast_download_url(data_uri),lazy_load_hdus=True,memmap=False,
                   use_fsspec=True,fsspec_kwargs={"block_size":1024*1024,"cache_type":"readahead"}) as h:
        for i,ext in enumerate(h):
            if ext.name not in ("SCI","PRIMARY"):
                continue
            nx=int(ext.header.get("NAXIS1",0))
            ny=int(ext.header.get("NAXIS2",0))
            if nx<=0 or ny<=0:
                continue
            try:
                w=WCS(ext.header)
                x,y=w.all_world2pix(RA,DEC,0)
                if np.isfinite([x,y]).all() and 0<=x<nx and 0<=y<ny:
                    return {"status":"COVERS_TARGET","detector":str(h[0].header.get("DETECTOR","")),
                            "mjd":float(h[0].header.get("EXPSTART",np.nan)),
                            "x_full":float(x),"y_full":float(y),"image_nx":int(nx),"image_ny":int(ny),
                            "sci_extension_index":i}
            except Exception:
                pass
    return {"status":"OUTSIDE_DETECTOR_OR_BAD_WCS"}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",default="publication_check/third_epoch_inventory")
    ap.add_argument("--radius",type=float,default=3.0)
    args=ap.parse_args()
    root=Path(args.output);root.mkdir(parents=True,exist_ok=True)
    obs=list_inventory(args.radius)
    basic=["obs_collection","instrument_name","filters","filter_norm","start_mjd","end_mjd",
           "obs_id","obsid","dataRights","red_relevant"]
    obs[basic].to_csv(root/"PUBLIC_MAST_OBSERVATIONS.csv",index=False)
    p=all_products(obs)
    if p.empty:
        (root/"RESULTS.md").write_text("No FITS products found. No third epoch established.\n")
        raise RuntimeError("No public candidate FITS products")
    join=obs.drop_duplicates("obsid")
    p=p.merge(join[["obsid","obs_collection","instrument_name","filter_norm",
                    "start_mjd","end_mjd"]],on="obsid",how="left")
    p["file_type"]=[
        product_type(row.productFilename,
          str(row.obs_collection)+"/"+str(row.instrument_name)) for _,row in p.iterrows()
    ]
    p["priority"]=p.filter_norm.isin(RELEVANT)
    meta=["obsid","obs_collection","instrument_name","filter_norm","start_mjd",
          "end_mjd","productFilename","dataURI","file_type","priority"]
    p[meta].to_csv(root/"ALL_PUBLIC_FITS_METADATA.csv",index=False)
    cand=p[(p.priority)&p.file_type.notna()].copy().sort_values(["start_mjd","filter_norm","productFilename"])
    cand[meta].to_csv(root/"RELEVANT_PRODUCT_CANDIDATES.csv",index=False)
    # For HST prefer single exposure FLT/FLC to mosaics, but preserve mosaics as fallback.
    receipts=[]
    for _,rec in cand.iterrows():
        fname=str(rec.productFilename);uri=str(rec.dataURI)
        print("Checking product WCS",rec.filter_norm,fname,flush=True)
        try:
            if rec.file_type=="JWST_CAL":
                result=check_jwst_wcs(fname,uri)
            else:
                result=check_hst_wcs(fname,uri)
        except Exception as exc:
            result={"status":"WCS_ERROR",
                    "reason":str(exc)[:500]}
        receipts.append({**{k:rec[k] for k in meta},**result})
        pd.DataFrame(receipts).to_csv(root/"PER_EXPOSURE_WCS_RECEIPTS.csv",index=False)
    checks=pd.DataFrame(receipts)
    good=checks[checks.status.eq("COVERS_TARGET")].copy()
    good.to_csv(root/"ACTUAL_COVERING_EXPOSURES.csv",index=False)
    if len(good):
        good["time_for_group"]=pd.to_numeric(good["mjd"],errors="coerce").fillna(pd.to_numeric(good["start_mjd"],errors="coerce"))
        good=good.sort_values("time_for_group")
        # Instruments/filters may vary; separate visits >30 d, but detections
        # and astrometric errors are not yet determined.
        good["visit_id"]=(good["time_for_group"].diff().fillna(1000)>30).cumsum()
        visits=good.groupby("visit_id").agg(first_mjd=("time_for_group","min"),
            last_mjd=("time_for_group","max"),
            n_exposures=("productFilename","size"),
            instruments=("instrument_name",lambda s:";".join(sorted(set(s.astype(str))))),
            filters=("filter_norm",lambda s:";".join(sorted(set(s.astype(str)))))).reset_index()
    else:
        visits=pd.DataFrame(columns=["visit_id","first_mjd","last_mjd","n_exposures","instruments","filters"])
    visits.to_csv(root/"POTENTIAL_EPOCHS_NOT_DETECTIONS.csv",index=False)
    report=[
        "# Third-epoch audit: candidate 1651050",
        "",
        "This audit distinguishes: public observation inventory, individual CAL "
        "coverage verified using WCS, and an *actual measured astrometric centroid*. "
        "Only the first two are tested here. **No new PM or parallax is measured.**",
        "",
        f"Public MAST observation rows: {len(obs)}.",
        f"Priority candidate individual CAL/HST products: {len(cand)}.",
        f"Original files with verified sky-position coverage: {len(good)}.",
        "",
        "## Candidate visits separated by >30 days (still NOT detections)",
        "",
        visits.to_markdown(index=False) if len(visits) else "None confirmed by detector WCS.",
        "",
        "## Scientific limitations",
        "- MIRI filters have significantly larger diffraction-limited PSFs than NIRCam: "
        "a nominally covered source must be detected at sufficient S/N and locally "
        "registered against field controls before it counts as a third epoch.",
        "- Blue HST/WFC3/ACS observations are extremely unlikely to show a "
        "500 K UCD and a non-detection cannot constrain its PM.",
        "- MAST records with a broad t_min/t_max interval may include visits "
        "separated by days; exact CAL EXPSTART is used where available.",
        "- A third astrometric epoch requires an independent, securely matched "
        "centroid with realistic cross-instrument/cross-filter systematic errors.",
        "- If no usable third centroid is found, do not claim a fully independent "
        "PM/parallax solution for 1651050.",
    ]
    (root/"RESULTS.md").write_text("\n".join(report)+"\n")
    print("\n".join(report))


if __name__=="__main__":
    main()
