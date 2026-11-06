"""Synthetic regression tests for conservative HST duplicate/anomaly triage."""
import pandas as pd
from scripts.triage_hst_1651050 import triage


def _row(root, name, filt, snr, bad=0, sigma=1.2, chi2=1.0, mjd=61055.0):
    return dict(rootname=root, filename=name, filter=filt,
                coverage_status="COVERS_WITH_STAMP", expstart_mjd=mjd,
                conservative_snr=snr, bad_aperture_pixels=bad,
                gauss_sigma_pix=sigma, gauss_chi2_reduced=chi2)


def test_archive_aliases_do_not_become_independent_exposures():
    rows = [_row("abc", "abc_flt.fits", "F160W", 0.8),
            _row("abc", "abc_alt_flt.fits", "F160W", 0.9),
            _row("def", "def_flt.fits", "F160W", 0.5)]
    roots, bands = triage(pd.DataFrame(rows))
    assert len(roots) == 2
    assert bands.loc[0, "unique_exposure_roots"] == 2
    assert not (roots.screening_status == "HIGH_SNR_REQUIRES_PSF_AND_REGISTRATION").any()


def test_high_forced_snr_with_bad_pixels_is_not_detection():
    rows = [_row("abc", "abc_flt.fits", "F300W", 13.47, bad=6,
                 sigma=0.7, chi2=3.87),
            _row("def", "def_flt.fits", "F300W", 0.3)]
    roots, bands = triage(pd.DataFrame(rows))
    assert roots.iloc[0].screening_status == "UNVALIDATED_HIGH_SNR_DQ_OR_FIT_ANOMALY"
    assert bands.loc[0, "roots_with_clean_repeatable_snr_ge_5"] == 0
    assert (bands.repeatability_status == "NO_VALIDATED_ASTROMETRIC_DETECTION").all()


def test_high_snr_without_flags_still_requires_registration():
    roots, bands = triage(pd.DataFrame([_row("abc", "abc_flt.fits", "F160W", 15)]))
    assert roots.iloc[0].screening_status == "HIGH_SNR_REQUIRES_PSF_AND_REGISTRATION"
    assert bands.iloc[0].repeatability_status == "NO_VALIDATED_ASTROMETRIC_DETECTION"


def test_outside_wcs_not_counted():
    df = pd.DataFrame([_row("abc", "abc_flt.fits", "F160W", 3),
                       _row("other", "other_flt.fits", "F160W", 99)])
    df.loc[1, "coverage_status"] = "OUTSIDE_DETECTOR"
    roots, _ = triage(df)
    assert len(roots) == 1
