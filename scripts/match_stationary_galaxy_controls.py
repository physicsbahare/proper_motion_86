#!/usr/bin/env python3
"""Cross-match Dave's field controls with the user's local COSMOS2025 FITS.

Run after downloading the dave-cross-filter-astrometry-audit GitHub artifact.
This is a vetting aid: a photo-z galaxy is NOT a spectroscopically confirmed
stationary object, and low measured relative motion alone does not confirm
an extragalactic nature.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

DEFAULT_CATALOG = "/run/media/bahareh/Seagate Hub/COSMOSWeb_mastercatalog_v1.fits"
REPO = "physicsbahare/proper_motion_86"
ARTIFACT_RUN = "41716345445"

def name_norm(s):
    return "".join(c for c in s.lower() if c.isalnum())

def find_name(names, *aliases):
    exact = {name_norm(s): s for s in names}
    for a in aliases:
        if name_norm(a) in exact:
            return exact[name_norm(a)]
    return None

def col(table, *aliases):
    name = find_name(table.names, *aliases)
    return None if name is None else np.asarray(table[name])

def num(v):
    if v is None:
        return np.nan
    try:
        if isinstance(v, bytes): v=v.decode(errors="replace")
        r=float(v)
        return r if np.isfinite(r) else np.nan
    except (ValueError, TypeError):
        return np.nan

def clean(v):
    if v is None: return None
    if isinstance(v, (bytes, np.bytes_)): v=v.decode(errors="replace")
    if isinstance(v, (float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,np.ndarray):return str(v)
    return v

def find_table(hdul, kind):
    if kind=="positions":
        for h in hdul[1:]:
            if not isinstance(h, fits.BinTableHDU): continue
            names=h.columns.names
            if find_name(names,"ra","RAJ2000","alpha_j2000","ra_deg") and find_name(names,"dec","DEJ2000","delta_j2000","dec_deg"):
                return h
    if kind=="lephare":
        for h in hdul[1:]:
            if isinstance(h, fits.BinTableHDU) and ("LEPHARE" in h.name.upper() or (
              find_name(h.columns.names,"zfinal","zpdf_med") is not None and
              find_name(h.columns.names,"type") is not None)):
                return h
    return None

def read_controls(controls_dir):
    root=Path(controls_dir)
    original=pd.read_csv(Path(__file__).resolve().parents[1]/"data/candidates_step3_magnitude_cut_strict_125.csv").set_index("candidate_id")
    allrows=[]
    for cid in (387132,1651050):
        fp=root/str(cid)/"leave_one_out_residuals.csv"
        if not fp.exists():raise FileNotFoundError(f"{fp} not found; obtain Action artifact {ARTIFACT_RUN}")
        data=pd.read_csv(fp)
        source=original.loc[cid]
        center=SkyCoord(float(source.ra_deg)*u.deg,float(source.dec_deg)*u.deg)
        for i,row in data.iterrows():
            if not np.isfinite([row.early_east_arcsec,row.early_north_arcsec]).all():continue
            sky=center.spherical_offsets_by(float(row.early_east_arcsec)*u.arcsec,float(row.early_north_arcsec)*u.arcsec)
            item={k:clean(v) for k,v in row.to_dict().items()}
            item.update({"candidate_id":cid,"local_control_index":int(i),
                         "control_ra_deg":float(sky.ra.deg),"control_dec_deg":float(sky.dec.deg)})
            allrows.append(item)
    return pd.DataFrame(allrows)

def auto_download(path):
    try:
        subprocess.run(["gh","run","download",ARTIFACT_RUN,
                        "--repo",REPO,"--name","dave-cross-filter-astrometry-audit",
                        "--dir",str(path)], check=True)
    except Exception as exc:
        raise RuntimeError("Could not download the previous field-control artifact. "
                           "Install/login with GitHub CLI (gh auth login), or download run "
                           f"{ARTIFACT_RUN} manually and pass --controls-dir. {exc}") from exc

def vet(catalog, controls, max_sep_arcsec=0.25):
    with fits.open(catalog,memmap=True,lazy_load_hdus=True) as hdul:
        phot_hdu=find_table(hdul,"positions")
        if phot_hdu is None:
            raise ValueError("No FITS table with both RA and Dec columns; inspect catalog extensions")
        phot=phot_hdu.data
        ra=col(phot,"ra","RAJ2000","alpha_j2000","ra_deg")
        dec=col(phot,"dec","DEJ2000","delta_j2000","dec_deg")
        if ra is None or dec is None:
            raise ValueError("Position columns missing in COSMOS2025 photometry")
        ra=np.asarray(ra,dtype=float);dec=np.asarray(dec,dtype=float)
        good=np.isfinite(ra)&np.isfinite(dec)&(dec>=-90)&(dec<=90)
        idx=np.flatnonzero(good)
        cosmos=SkyCoord(ra[idx]*u.deg,dec[idx]*u.deg)
        observed=SkyCoord(controls.control_ra_deg.to_numpy(float)*u.deg,
                          controls.control_dec_deg.to_numpy(float)*u.deg)
        match_index,sep,_=observed.match_to_catalog_sky(cosmos)
        cosmos_indices=idx[match_index]
        leph_hdu=find_table(hdul,"lephare")
        lep=leph_hdu.data if leph_hdu is not None else None
        aligned=False
        if lep is not None and len(lep)==len(phot):
            phot_id=col(phot,"id","source_id")
            lep_id=col(lep,"id","source_id")
            if phot_id is not None and lep_id is not None:
                aligned=bool(np.array_equal(phot_id[:min(1000,len(phot))],lep_id[:min(1000,len(lep))]) and
                             np.array_equal(phot_id[-min(1000,len(phot)):],lep_id[-min(1000,len(lep)):]))
            else:
                # same row count alone is NOT proof of aligned row IDs
                aligned=False
        fields={
          "id":("id","source_id"),
          "flag_star":("flag_star",),
          "flag_blend":("flag_blend",),
          "zspec":("zspec","z_spec","z_spectro","spec_z","zspecbest","specz"),
          "reff":("re_f444w","r_eff_f444w","radius_50","flux_radius_f444w"),
          "class_star":("class_star","class_star_f444w"),
          "f444w_snr":("snr_f444w","s_n_f444w"),
        }
        rec=[]
        for i,(original_idx,s) in enumerate(zip(cosmos_indices,sep.arcsec)):
            row=controls.iloc[i].to_dict()
            row.update({"matched_phot_hdu":phot_hdu.name,"catalog_match_sep_arcsec":float(s),
                        "catalog_position_match":bool(s<=max_sep_arcsec),
                        "catalog_index":int(original_idx),
                        "lephare_row_alignment_verified":bool(aligned)})
            for key,aliases in fields.items():
                cc=find_name(phot.names,*aliases)
                row[key]=clean(phot[cc][original_idx]) if cc is not None else None
            if lep is not None and aligned:
                typecol=find_name(lep.names,"type","lephare_type")
                zcol=find_name(lep.names,"zfinal","zpdf_med","zphot")
                row["lephare_type"]=clean(lep[typecol][original_idx]) if typecol else None
                row["lephare_zphot"]=clean(lep[zcol][original_idx]) if zcol else None
                zspeclep=find_name(lep.names,"zspec","z_spec","specz","spec_z")
                if zspeclep is not None and not np.isfinite(num(row["zspec"])):
                    row["zspec"]=clean(lep[zspeclep][original_idx])
            else:
                row["lephare_type"]=None;row["lephare_zphot"]=None
            exists=bool(s<=max_sep_arcsec)
            galaxy=exists and num(row["lephare_type"])==0
            spec=exists and num(row["zspec"])>0.01
            not_bad=exists and (row["flag_blend"] is None or num(row["flag_blend"])==0)
            eligible=galaxy and not_bad
            row["galaxy_photometric_candidate"]=bool(eligible)
            row["spectroscopic_redshift_support"]=bool(spec)
            # Do not call spec-z a "confirmation" unless its measured reliability is checked.
            row["stationarity_status"]="NOT_INDEPENDENTLY_CONFIRMED"
            row["vetting_level"]=("REVIEW_SPECZ_GALAXY" if eligible and spec else
                                 "REVIEW_PHOTOMETRIC_GALAXY" if eligible else
                                 "UNKNOWN_OR_REJECTED")
            rec.append(row)
        output=pd.DataFrame(rec)
        info={
          "catalog":str(catalog),"phot_hdu":phot_hdu.name,"phot_n":len(phot),
          "lephare_hdu":leph_hdu.name if leph_hdu else None,
          "lephare_row_alignment_verified":aligned,
          "catalog_match_radius_arcsec":max_sep_arcsec,
          "matched":int(output.catalog_position_match.sum()),
          "possible_photometric_galaxies":int(output.galaxy_photometric_candidate.sum()),
          "specz_galaxy_candidates":int((output.galaxy_photometric_candidate&output.spectroscopic_redshift_support).sum()),
          "confirmed_stationary_galaxies":0,
          "warning":"Spectroscopic evidence is unverified and no independent stationarity test is yet applied. Never adopt the full control covariance as an isolated chromatic error."
        }
        return output,info

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--catalog",default=DEFAULT_CATALOG)
    ap.add_argument("--controls-dir",default="results_cross_filter_null")
    ap.add_argument("--download-controls",action="store_true")
    ap.add_argument("--out",default="results_stationary_galaxy_controls")
    ap.add_argument("--match-radius-arcsec",type=float,default=0.25)
    args=ap.parse_args()
    catalog=Path(args.catalog)
    if not catalog.exists():ap.error(f"Local FITS missing: {catalog} (this path must exist on the machine executing the script)")
    controls_dir=Path(args.controls_dir)
    if args.download_controls and not any(controls_dir.glob("*/leave_one_out_residuals.csv")):auto_download(controls_dir)
    sources=read_controls(controls_dir)
    result,summary=vet(catalog,sources,args.match_radius_arcsec)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    result.to_csv(out/"cosmos2025_matched_field_controls.csv",index=False)
    summary["by_target"]={}
    for cid,group in result.groupby("candidate_id"):
        summary["by_target"][str(cid)]={
            "n_field_controls":len(group),
            "n_position_matched":int(group.catalog_position_match.sum()),
            "n_photometric_galaxy_candidates":int(group.galaxy_photometric_candidate.sum()),
            "n_specz_galaxy_candidates":int((group.galaxy_photometric_candidate & group.spectroscopic_redshift_support).sum()),
            "n_confirmed_stationary":0
        }
    (out/"summary.json").write_text(json.dumps(summary,indent=2,default=str))
    print(json.dumps(summary,indent=2,default=str))

if __name__=="__main__":main()
