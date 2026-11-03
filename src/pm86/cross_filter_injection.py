"""Cross-filter injection primitives. No detection-significance claims.

This module changes in-memory CAL cutouts only. A PSF template must be supplied
separately for each filter; it is *not* inferred from the target and does not
purport to calibrate chromatic PSF effects. All positions are sky coordinates.
"""
from __future__ import annotations

from dataclasses import replace
import numpy as np


def normalized_psf(psf):
    p = np.asarray(psf, dtype=float)
    if p.ndim != 2 or min(p.shape) < 3 or not np.all(np.isfinite(p)):
        raise ValueError("PSF must be a finite two-dimensional stamp")
    if np.any(p < 0) or not np.isfinite(p.sum()) or p.sum() <= 0:
        raise ValueError("PSF must be nonnegative and have positive integral")
    return p / p.sum()


def inject_at_sky(exposure, ra_deg, dec_deg, psf, flux, *, rng=None,
                  poisson_sigma=None, strict=True):
    """Return an independent in-memory exposure with a sky-placed fake source.

    Parameters
    ----------
    exposure : pm86.archive.ExposureCutout
        CAL SCI/ERR/DQ cutout and original gWCS. Flux is in SCI image units.
    psf : array
        Filter-specific, normalized or unnormalized nonnegative PSF sampled on
        the *native detector grid*. An externally supplied ePSF shifted to the
        requested subpixel phase is required for quantitative science.
    poisson_sigma : float or None
        Optional *per-pixel* additional injected-noise sigma in SCI units.
        Without it, the injected signal is deterministic and the unchanged ERR
        does not include its shot noise. This is not an absolute noise model.
    strict : bool
        Reject truncated stamps and masked pixels, rather than silently
        changing injected flux or contaminating invalid pixels.
    """
    if not np.isfinite(flux) or flux <= 0:
        raise ValueError("positive finite injected flux required")
    p = normalized_psf(psf)
    # gWCS uses full detector pixel coordinates, not cutout coordinates.
    xfull, yfull = exposure.gwcs.invert(float(ra_deg), float(dec_deg))
    x = float(np.asarray(xfull).squeeze()) - exposure.x0
    y = float(np.asarray(yfull).squeeze()) - exposure.y0
    if not (np.isfinite(x) and np.isfinite(y)):
        raise ValueError("sky coordinate outside valid WCS")
    ny, nx = exposure.sci.shape
    # Integer stamp placement is a *prototype only*: its subpixel
    # quantization makes it unsuitable for a mas-scale null test.
    cx, cy = int(np.floor(x + 0.5)), int(np.floor(y + 0.5))
    sy, sx = p.shape
    x0, y0 = cx - sx // 2, cy - sy // 2
    x1, y1 = x0 + sx, y0 + sy
    if x0 < 0 or y0 < 0 or x1 > nx or y1 > ny:
        raise ValueError("PSF stamp extends outside CAL cutout")
    dq = np.asarray(exposure.dq)
    if strict and np.any(dq[y0:y1, x0:x1] != 0):
        raise ValueError("PSF stamp overlaps DQ-flagged pixels")
    if strict and (not np.all(np.isfinite(exposure.sci[y0:y1, x0:x1]))
                   or not np.all(np.isfinite(exposure.err[y0:y1, x0:x1]))):
        raise ValueError("PSF stamp overlaps nonfinite pixels")
    sci = np.array(exposure.sci, dtype=float, copy=True)
    err = np.array(exposure.err, dtype=float, copy=True)
    signal = float(flux) * p
    sci[y0:y1, x0:x1] += signal
    if poisson_sigma is not None:
        if not np.isfinite(poisson_sigma) or poisson_sigma < 0:
            raise ValueError("invalid poisson_sigma")
        if rng is None:
            raise ValueError("seeded RNG required for noisy injection")
        # Noise realization and ERR update are separate from the real CAL
        # noise; this is a user-specified approximation, not an e-/s model.
        sci[y0:y1, x0:x1] += rng.normal(0, poisson_sigma, size=p.shape)
        err[y0:y1, x0:x1] = np.hypot(err[y0:y1, x0:x1], poisson_sigma)
    return replace(exposure, sci=sci, err=err, ra_deg=float(ra_deg),
                   dec_deg=float(dec_deg), x_full=float(xfull),
                   y_full=float(yfull)), {
        "x_full": float(xfull), "y_full": float(yfull),
        "x_cutout": x, "y_cutout": y,
        "subpixel_phase_x": x - cx, "subpixel_phase_y": y - cy,
        "flux_requested": float(flux), "flux_injected": float(signal.sum()),
        "subpixel_psf_applied": False,
        "valid_for_mas_astrometry": False,
    }


def empirical_tail(null_statistic, observed):
    """One-sided empirical exceedance with +1 correction; never a Gaussian sigma."""
    vals = np.asarray(null_statistic, dtype=float)
    vals = vals[np.isfinite(vals)]
    if not len(vals) or not np.isfinite(observed):
        raise ValueError("finite null samples and observed statistic required")
    return {"n_null": int(len(vals)),
            "n_exceed": int(np.count_nonzero(vals >= observed)),
            "p_empirical": float((np.count_nonzero(vals >= observed) + 1) / (len(vals) + 1)),
            "p95": float(np.quantile(vals, 0.95)),
            "p99": float(np.quantile(vals, 0.99)),
            "minimum_resolvable_p": float(1 / (len(vals) + 1))}
