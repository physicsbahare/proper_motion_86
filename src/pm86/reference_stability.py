"""Reference-source geometry and covariance stability diagnostics.

These checks do not classify galaxies or infer a calibrated systematic floor.
They quantify whether the small-N covariance is dominated by a single control.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def covariance_stability(residuals: pd.DataFrame, min_controls: int = 7) -> dict:
    required = {"east_mas", "north_mas", "retained_for_covariance"}
    if not required.issubset(residuals.columns):
        return {"status": "MISSING_REQUIRED_COLUMNS"}
    retained = residuals.loc[residuals.retained_for_covariance.astype(bool)].copy()
    xy = retained[["east_mas", "north_mas"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    xy = xy[np.isfinite(xy).all(axis=1)]
    n = len(xy)
    if n < min_controls + 1:
        return {"status": "INSUFFICIENT_FOR_JACKKNIFE", "n": n}
    full = np.cov(xy, rowvar=False)
    if not np.isfinite(full).all() or np.linalg.eigvalsh(full).min() <= 0:
        return {"status": "DEGENERATE_COVARIANCE", "n": n}
    jack = np.array([np.cov(np.delete(xy, i, axis=0), rowvar=False) for i in range(n)])
    if not np.isfinite(jack).all():
        return {"status": "INVALID_JACKKNIFE", "n": n}
    base = np.sqrt(np.diag(full))
    scatter = np.sqrt(np.maximum(np.diagonal(jack, axis1=1, axis2=2), 0))
    # This is a stability diagnostic, not the standard error of a systematic.
    ranges = np.stack([scatter.min(axis=0), scatter.max(axis=0)], axis=1)
    max_ratio = float(np.max(np.abs(scatter / base - 1)))
    result = {
        "status": "DIAGNOSTIC_ONLY",
        "n_retained": n,
        "full_scatter_east_north_mas": base.tolist(),
        "jackknife_scatter_range_east_mas": ranges[0].tolist(),
        "jackknife_scatter_range_north_mas": ranges[1].tolist(),
        "maximum_fractional_scatter_change": max_ratio,
        "single_source_sensitive_over_20pct": bool(max_ratio > 0.20),
        "not_a_systematic_uncertainty_estimate": True,
    }
    pos = {"early_east_arcsec", "early_north_arcsec"}
    if pos.issubset(retained.columns):
        points = retained[["early_east_arcsec", "early_north_arcsec"]].to_numpy(float)
        if np.isfinite(points).all() and len(points) >= 3:
            centered = points - points.mean(axis=0)
            s = np.linalg.svd(centered, compute_uv=False)
            result["reference_geometry_minor_major_singular_ratio"] = (
                float(s[-1] / s[0]) if s[0] > 0 else 0.0
            )
            result["reference_geometry_near_collinear"] = bool(s[0] == 0 or s[-1] / s[0] < 0.1)
    return result
