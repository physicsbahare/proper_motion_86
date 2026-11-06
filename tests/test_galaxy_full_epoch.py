import numpy as np
import pandas as pd

from pm86.galaxy_full_epoch import (
    calculate_affine, registered_solution, summarize_full_epoch_null
)

def make_case(n=8, late_frames=1):
    rng=np.random.default_rng(20261106)
    ids=list(range(101,101+n))
    xy=rng.uniform(-8,8,(n,2))
    epochs=[("e0","early",61000.,[0,0]),("e1","early",61000.1,[.015,-.010])]
    for j in range(late_frames):
        epochs.append((f"l{j}","late",61500.+j*.1,[-.020,.012]))
    records=[];control_map={}
    for name,epoch,mjd,shift in epochs:
        eps=np.array(shift)
        pts=xy+eps+rng.normal(0,.00015,size=(n,2))
        frame=pd.DataFrame({
            "galaxy_id":ids,"east_arcsec_raw_wcs":pts[:,0],
            "north_arcsec_raw_wcs":pts[:,1],"snr":np.full(n,30.),
        })
        control_map[name]=frame
        pm=(np.array([-.030,.010]) if epoch=="late" else np.zeros(2))
        center=np.array([1.2,2.1])+eps+pm
        records.append({
            "filename":name,"epoch_label":epoch,"mjd":mjd,"filter":"F444W" if epoch=="early" else "F356W",
            "east_2dg_arcsec_raw_wcs":center[0],
            "north_2dg_arcsec_raw_wcs":center[1],
            "east_com_arcsec_raw_wcs":center[0],
            "north_com_arcsec_raw_wcs":center[1],
            "sigma_x_pix":.06,"sigma_y_pix":.06,"pixel_scale_mas":63.0,
            "astrometric_snr":35.0,"clean_snr_err":35.0
        })
    return pd.DataFrame(records),control_map


def test_affine_uses_only_other_sources():
    frames,controls=make_case()
    coef,summary=calculate_affine(controls["e0"],controls["e1"],exclude_id=101)
    assert summary["status"]=="OK"
    assert summary["n_common"]==7
    assert 101 not in summary["galaxy_ids_fit"]
    assert np.all(np.isfinite(coef))


def test_full_epoch_heldout_null_preserves_actual_dither_count():
    frames,controls=make_case(late_frames=2)
    r,q,heldout,summary=registered_solution(frames,controls,"e0",heldout=101)
    assert len(r)==len(frames)
    assert summary["n_early_registered"]==2
    assert summary["n_late_registered"]==2
    assert summary["galaxy_appears_in_all_original_target_usable_exposures"]
    assert summary["heldout_galaxy_null"]["status"]=="FULL_EPOCH_GALAXY_NULL_MEASURED"
    assert abs(summary["heldout_galaxy_null"]["delta_east_mas"])<5
    assert abs(summary["heldout_galaxy_null"]["delta_north_mas"])<5
    assert abs(summary["target_refit_pm"]["mu_alpha_cosdec_masyr"])>10


def test_single_late_frame_is_not_hidden_and_zero-motion_false_tails_not_quoted():
    frames,controls=make_case(late_frames=1)
    outputs=[]
    for galaxy_id in controls["e0"].galaxy_id.values:
        _,_,_,summary=registered_solution(frames,controls,"e0",heldout=int(galaxy_id))
        outputs.append(summary)
        assert summary["n_late_registered"]==1
    result=summarize_full_epoch_null(outputs)
    assert result["n_full_exposure_galaxies"]==8
    assert result["n_at_least_two_epoch_galaxies"]==8
    assert result["status"]=="PRELIMINARY_FULL_EPOCH_GALAXY_DIAGNOSTIC"
    assert result["null_calibrated"] is False
    assert result["n_confirmed_stationary"]==0


def test_underconstrained_fit_does_not_become_valid_zero_motion():
    frames,controls=make_case(n=6)
    _,_,_,summary=registered_solution(frames,controls,"e0",heldout=101)
    assert summary["n_late_registered"]==0
    assert summary["heldout_galaxy_null"]["status"]=="MISSING_EPOCH"
    result=summarize_full_epoch_null([summary])
    assert result["status"]=="INSUFFICIENT_FULL_EXPOSURE_GALAXY_NULLS"
