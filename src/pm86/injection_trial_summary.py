"""Trial-level summary for Shabgard-style cross-filter injection/recovery.

Pure numerical diagnostics; these DO NOT establish that any CAL image has
been injected, or that the PSF/registration model is scientifically validated.
Each row must correspond to an independent injected sky position, not a dither.
"""
import numpy as np
from pm86.cross_filter_injection import empirical_tail


def summarize_trials(stationary_pm, moving_pm, true_moving_pm, observed_pm):
    """Summarize fixed-size (N,2) PM trial arrays in mas/yr.

    Rows with either component nonfinite are failed recoveries. The denominator
    for recovery fraction includes failures; only recovered rows enter percentiles.
    Null amplitudes are computed from the *vector*, not quadrature component S/N.
    """
    null = np.asarray(stationary_pm, dtype=float)
    move = np.asarray(moving_pm, dtype=float)
    truth = np.asarray(true_moving_pm, dtype=float)
    observed = np.asarray(observed_pm, dtype=float)
    if null.ndim != 2 or move.ndim != 2 or null.shape[1] != 2 or move.shape[1] != 2:
        raise ValueError("trial PM arrays must have shape (N,2)")
    if not len(null) or not len(move) or truth.shape != (2,) or observed.shape != (2,):
        raise ValueError("nonempty trials and two-component true/observed PM required")
    if not np.all(np.isfinite(truth)) or not np.all(np.isfinite(observed)):
        raise ValueError("true/observed PM must be finite")
    good_null = np.all(np.isfinite(null), axis=1)
    good_move = np.all(np.isfinite(move), axis=1)
    if not good_null.any() or not good_move.any():
        raise ValueError("at least one valid recovery in each trial population required")
    n = null[good_null]
    m = move[good_move]
    amplitudes = np.linalg.norm(n, axis=1)
    tail = empirical_tail(amplitudes, np.linalg.norm(observed))
    residual = m - truth
    return {
        "n_null_attempted": int(len(null)),
        "n_null_recovered": int(good_null.sum()),
        "null_recovery_fraction": float(good_null.mean()),
        "n_moving_attempted": int(len(move)),
        "n_moving_recovered": int(good_move.sum()),
        "moving_recovery_fraction": float(good_move.mean()),
        "null_bias_east_north_masyr": np.mean(n, axis=0).tolist(),
        "null_covariance_masyr2": (np.cov(n, rowvar=False).tolist() if len(n) > 1 else None),
        "moving_recovery_bias_east_north_masyr": np.mean(residual, axis=0).tolist(),
        "moving_residual_p16_p84_masyr": np.quantile(residual, [0.16, 0.84], axis=0).tolist(),
        "stationary_amplitude_tail": tail,
        "tail_conditional_on_success": True,
        "valid_for_publication": False,
        "warning": "Recovery failure/selection, PSF mismatch, source independence, "
                   "and full-epoch registration must be audited before interpretation.",
    }
