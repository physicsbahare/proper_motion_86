import numpy as np
from scripts.audit_hst_flt_1651050 import forced_aperture, gaussian_screen


def test_forced_aperture_synthetic_source():
    yy, xx = np.indices((37, 37))
    sci = 2 + 10 * np.exp(-0.5 * ((xx - 18)**2 + (yy - 18)**2) / 1.2**2)
    err = np.ones_like(sci)
    dq = np.zeros_like(sci, dtype=np.uint16)
    out = forced_aperture(sci, err, dq, 18, 18)
    assert out["status"] == "FORCED_PHOTOMETRY_ONLY"
    assert out["conservative_snr"] > 5
    fit = gaussian_screen(sci, err, dq, 18, 18)
    assert fit["fit_status"] == "EXPLORATORY_GAUSSIAN_NOT_ASTROMETRY"
    assert fit["gauss_offset_from_wcs_pix"] < 0.1


def test_no_source_is_not_positive_detection():
    sci = np.ones((37, 37))
    err = np.ones_like(sci)
    dq = np.zeros_like(sci, dtype=np.uint16)
    out = forced_aperture(sci, err, dq, 18, 18)
    assert abs(out["conservative_snr"]) < 1e-6


def test_bad_pixels_fail_closed():
    sci = np.ones((37, 37))
    err = np.ones_like(sci)
    dq = np.ones_like(sci, dtype=np.uint16)
    out = forced_aperture(sci, err, dq, 18, 18)
    assert out["status"] == "INSUFFICIENT_VALID_PIXELS"
