import numpy as np
import pandas as pd
from pm86.galaxy_full_epoch import original_registration_galaxy_null,summarize_full_epoch_null

def setup_data():
    rng=np.random.default_rng(4319)
    ids=np.arange(101,111)
    xy=rng.uniform(-7,7,(len(ids),2))
    measures=[];raw={};vetted={}
    obs=[("e1","early",61000.,[0,0]),("e2","early",61000.2,[.017,-.008]),
         ("e3","early",61000.4,[-.009,.011]),("e4","early",61000.6,[.014,.004]),
         ("l1","late",61500.,[-.015,.018]),("l2","late",61500.1,[.008,-.011])]
    for name,epoch,mjd,sh in obs:
        shift=np.array(sh);coords=xy+shift+rng.normal(0,.0001,xy.shape)
        c=pd.DataFrame({
            "filename":name,"galaxy_id":ids,"snr":np.full(len(ids),45.0),
            "flux":np.full(len(ids),100.0),
            "east_arcsec_raw_wcs":coords[:,0],"north_arcsec_raw_wcs":coords[:,1]})
        raw[name]=c.copy()
        vetted[name]=c.copy()
        true_motion=np.array([-.030,.010]) if epoch=="late" else np.array([0.,0.])
        t=np.array([1.3,2.2])+shift+true_motion
        measures.append({
          "filename":name,"epoch_label":epoch,"mjd":mjd,
          "filter":"F444W" if epoch=="early" else "F410M",
          "east_2dg_arcsec_raw_wcs":t[0],"north_2dg_arcsec_raw_wcs":t[1],
          "east_com_arcsec_raw_wcs":t[0],"north_com_arcsec_raw_wcs":t[1],
          "sigma_x_pix":.06,"sigma_y_pix":.06,"pixel_scale_mas":63.,
          "astrometric_snr":40.,"clean_snr_err":40.})
    return pd.DataFrame(measures),raw,vetted


def test_mixed_frame_complete_heldout_galaxy_null_and_reference_independence():
    measurements,raw,galaxy=setup_data()
    r,q,s,summary=original_registration_galaxy_null(measurements,raw,galaxy,heldout_id=101)
    assert len(r)==len(measurements)==6
    assert len(s)==6
    assert summary["n_early_registered"]==4 and summary["n_late_registered"]==2
    assert summary["galaxy_appears_in_all_original_target_usable_exposures"]
    assert summary["heldout_galaxy_null"]["status"]=="FULL_EPOCH_GALAXY_NULL_MEASURED"
    assert abs(summary["heldout_galaxy_null"]["delta_east_mas"])<2
    assert abs(summary["heldout_galaxy_null"]["delta_north_mas"])<2
    assert summary["registration_method"]=="ORIGINAL_MIXED_FIELD_CONTROLS_WITH_GALAXY_LEFT_OUT"
    assert all(101 not in c.galaxy_id.values for c in [])  # the exclusion is exercised in registration


def test_original_target_fit_near_expected_motion_and_heldout_sources_not_independent_dithers():
    measurements,raw,galaxy=setup_data()
    allr,_,_,base=original_registration_galaxy_null(measurements,raw,galaxy)
    assert len(allr)==6
    assert abs(base["target_refit_pm"]["mu_alpha_cosdec_masyr"]+21.92)<1.0
    outputs=[original_registration_galaxy_null(measurements,raw,galaxy,int(g))[3]
             for g in range(101,111)]
    summ=summarize_full_epoch_null(outputs)
    assert summ["n_full_exposure_galaxies"]==10
    assert summ["null_calibrated"] is False
    assert summ["minimum_tail_resolution"]==1/11


def test_original_reference_galaxy_lost_in_one_late_exposure_is_not_full_coverage():
    measurements,raw,galaxy=setup_data()
    galaxy["l2"]=galaxy["l2"][galaxy["l2"].galaxy_id.ne(101)]
    _,_,_,summary=original_registration_galaxy_null(measurements,raw,galaxy,101)
    assert summary["n_registered_exposures"]==6
    assert summary["galaxy_appears_in_all_original_target_usable_exposures"] is False
    assert summary["heldout_galaxy_null"]["status"]=="FULL_EPOCH_GALAXY_NULL_MEASURED"
