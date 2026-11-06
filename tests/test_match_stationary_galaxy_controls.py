from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits

from scripts.match_stationary_galaxy_controls import vet, find_name, name_norm

def test_name_matching():
    assert name_norm("Z_SPEC")=="zspec"
    assert find_name(["RA","DEC","Z_SPEC"],"z_spec")=="Z_SPEC"

def test_catalog_crossmatch_does_not_claim_stationarity(tmp_path):
    phot=fits.BinTableHDU.from_columns([
        fits.Column(name="RA",format="D",array=np.array([150.1,150.2])),
        fits.Column(name="DEC",format="D",array=np.array([2.0,2.1])),
        fits.Column(name="ID",format="K",array=np.array([101,102])),
        fits.Column(name="flag_blend",format="K",array=np.array([0,0])),
    ],name="PHOTOMETRY")
    leph=fits.BinTableHDU.from_columns([
        fits.Column(name="ID",format="K",array=np.array([101,102])),
        fits.Column(name="TYPE",format="K",array=np.array([0,1])),
        fits.Column(name="zfinal",format="D",array=np.array([1.12,0.0])),
        fits.Column(name="z_spec",format="D",array=np.array([1.1,np.nan])),
    ],name="LEPHARE")
    catalog=tmp_path/"test.fits"
    fits.HDUList([fits.PrimaryHDU(),phot,leph]).writeto(catalog)
    controls=pd.DataFrame([
        {"candidate_id":387132,"control_ra_deg":150.1,"control_dec_deg":2.0,
         "east_mas":2.0,"north_mas":-1.0},
        {"candidate_id":1651050,"control_ra_deg":150.2,"control_dec_deg":2.1,
         "east_mas":4.0,"north_mas":3.0},
    ])
    results,summary=vet(catalog,controls,max_sep_arcsec=0.25)
    assert summary["matched"]==2
    assert summary["possible_photometric_galaxies"]==1
    assert summary["specz_galaxy_candidates"]==1
    assert summary["confirmed_stationary_galaxies"]==0
    assert results.iloc[0].vetting_level=="REVIEW_SPECZ_GALAXY"
    assert results.iloc[0].stationarity_status=="NOT_INDEPENDENTLY_CONFIRMED"
    assert results.iloc[1].vetting_level=="UNKNOWN_OR_REJECTED"
