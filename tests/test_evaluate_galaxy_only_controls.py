import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
import pandas as pd
from evaluate_galaxy_only_controls import galaxy_loo, positive_bool


def dummy_group(n=9):
    rng=np.random.default_rng(456)
    src=rng.normal(0,3,size=(n,2))
    shift=np.array([0.012,-0.006])
    dst=src+shift+rng.normal(0,.001,size=(n,2))
    frame=pd.DataFrame({
        "candidate_id":[387132]*n,
        "local_control_index":range(n),
        "catalog_position_match":[True]*n,
        "galaxy_photometric_candidate":[True]*(n-1)+[False],
        "retained_for_covariance":[True]*n,
        "snr_early":[25.0]*n,"snr_late":[20.0]*n,
        "early_east_arcsec":dst[:,0],"early_north_arcsec":dst[:,1],
        "late_east_arcsec":src[:,0],"late_north_arcsec":src[:,1],
        "east_mas":[1.0]*n, "north_mas":[-2.0]*n,
        "lephare_type":[0]*(n-1)+[1],
        "lephare_zphot":[1.1]*n,
        "id":range(1,n+1)
    })
    return frame


def test_galaxy_loo_uses_only_galaxies_and_no_claim_of_stationarity():
    results,report=galaxy_loo(dummy_group())
    assert len(results)==8
    assert report["status"]=="PRELIMINARY_GALAXY_ONLY_PAIR_NULL"
    assert report["n_photometric_galaxies"]==8
    assert report["spectroscopically_confirmed_stationary"]==0
    assert report["full_epoch_systematics_calibrated"] is False
    assert len(report["covariance_mas2"])==2


def test_fails_closed_on_too_few_galaxies():
    group=dummy_group()
    group.loc[2:,"galaxy_photometric_candidate"]=False
    results,report=galaxy_loo(group)
    assert results.empty
    assert report["status"]=="INSUFFICIENT_GALAXY_CONTROLS_FOR_LOO"


def test_bool_values_from_csv():
    vals=pd.Series(["True","False","true","0","1"])
    assert positive_bool(vals).tolist()==[True,False,True,False,True]
