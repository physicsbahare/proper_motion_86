import pandas as pd
from pm86.refit_reconciliation import reconcile_frozen_refit


def test_frozen_refit_comparison_does_not_invent_independent_sigma():
    frozen = {
        "baseline_days": 497.6405,
        "mu_alpha_cosdec_masyr": -55.2303,
        "sigma_mu_alpha_cosdec_masyr": 11.5336,
        "mu_delta_masyr": 16.4082,
        "sigma_mu_delta_masyr": 9.9982,
    }
    refit = {
        "pm_status": "PM_MEASURED", "classification": "AMBIGUOUS",
        "baseline_days": 497.6405,
        "mu_alpha_cosdec_masyr": -61.74,
        "sigma_mu_alpha_cosdec_masyr": 13.0,
        "mu_delta_masyr": 17.15,
        "sigma_mu_delta_masyr": 11.0,
    }
    m = pd.DataFrame([{
        "filename": "a_cal.fits", "epoch_label": "early", "filter": "F444W",
        "mjd": 60000., "astrometric_snr": 4.0,
        "east_2dg_arcsec_raw_wcs": 0.01,
        "north_2dg_arcsec_raw_wcs": 0.02,
        "astrometry_source": "catalog_position_2dg",
    }, {
        "filename": "b_cal.fits", "epoch_label": "late", "filter": "F356W",
        "mjd": 60497.64, "astrometric_snr": 2.0,
        "east_2dg_arcsec_raw_wcs": 0.03,
        "north_2dg_arcsec_raw_wcs": 0.04,
        "astrometry_source": "forced_gaussian",
    }])
    r = reconcile_frozen_refit(frozen, refit, m, m.iloc[:1])
    assert abs(r["delta_mu_east_masyr"] - (-6.5097)) < 1e-4
    assert r["quadrature_difference_significance_invalid"]
    assert "difference_sigma" not in r
    assert r["usable_by_epoch"] == {"early": 1, "late": 0}
    assert r["registered_by_epoch"] == {"early": 1, "late": 0}
