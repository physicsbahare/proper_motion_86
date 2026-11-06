from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from match_stationary_galaxy_controls import vet, find_name, name_norm

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


def test_official_cosmos2025_row_order_without_lephare_id(tmp_path):
    """The real COSMOS2025 LEPHARE HDU usually has no ID column."""
    phot=fits.BinTableHDU.from_columns([
        fits.Column(name="ra",format="D",array=np.array([150.1,150.2])),
        fits.Column(name="dec",format="D",array=np.array([2.0,2.1])),
        fits.Column(name="id",format="K",array=np.array([101,102])),
        fits.Column(name="flag_blend",format="L",array=np.array([False,False]))
    ],name="PHOTOMETRY HOTCOLD AND SE++")
    lep=fits.BinTableHDU.from_columns([
        fits.Column(name="type",format="K",array=np.array([0,1])),
        fits.Column(name="zfinal",format="D",array=np.array([1.12,0.0]))
    ],name="LEPHARE")
    path=tmp_path/"cosmos2025.fits"
    fits.HDUList([fits.PrimaryHDU(),phot,lep]).writeto(path)
    sources=pd.DataFrame([
        {"candidate_id":387132,"control_ra_deg":150.1,"control_dec_deg":2.0},
        {"candidate_id":1651050,"control_ra_deg":150.2,"control_dec_deg":2.1}
    ])
    result,summary=vet(path,sources)
    assert summary["lephare_join_method"]=="COSMOS2025_OFFICIAL_ROW_ORDER_EQUAL_LENGTH"
    assert summary["n_lephare_joined"]==2
    assert summary["possible_photometric_galaxies"]==1
    assert summary["specz_galaxy_candidates"]==0
    assert summary["confirmed_stationary_galaxies"]==0
    assert result.iloc[0].lephare_zphot==1.12

def test_reordered_lephare_with_unique_source_ids(tmp_path):
    phot=fits.BinTableHDU.from_columns([
        fits.Column(name="ra",format="D",array=np.array([150.1,150.2])),
        fits.Column(name="dec",format="D",array=np.array([2.0,2.1])),
        fits.Column(name="id",format="K",array=np.array([101,102]))
    ],name="PHOTOMETRY HOTCOLD AND SE++")
    lep=fits.BinTableHDU.from_columns([
        fits.Column(name="id",format="K",array=np.array([102,101])),
        fits.Column(name="type",format="K",array=np.array([1,0])),
        fits.Column(name="zfinal",format="D",array=np.array([0.0,1.42]))
    ],name="LEPHARE")
    path=tmp_path/"reordered.fits"
    fits.HDUList([fits.PrimaryHDU(),phot,lep]).writeto(path)
    source=pd.DataFrame([{"candidate_id":387132,"control_ra_deg":150.1,"control_dec_deg":2.0}])
    matched,summary=vet(path,source)
    assert summary["lephare_join_method"]=="UNIQUE_ID_JOIN"
    assert summary["possible_photometric_galaxies"]==1
    assert matched.iloc[0].lephare_zphot==1.42
