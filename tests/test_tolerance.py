from pathlib import Path

import numpy as np
import pandas as pd

from pm86.tolerance import (
    add_tolerance_flags,
    collect_tolerance_rows,
    summarize_tolerance,
    tolerance_sweep,
)


def sample_detail():
    return pd.DataFrame({
        "candidate_id": [1, 1, 2, 3],
        "clean_snr_err": [10.0, 4.0, 2.5, 8.0],
        "method_sep_2dg_com_mas": [12.0, 72.0, 100.0, 59.0],
        "registration_shift_2dg_mas": [5.0, 20.0, 90.0, 61.0],
    })


def test_60mas_method_and_registration_flags_only_use_snr3_rows():
    flagged = add_tolerance_flags(sample_detail(), 60.0)
    assert flagged["method_sensitive_at_tolerance"].tolist() == [False, True, False, False]
    assert flagged["registration_sensitive_at_tolerance"].tolist() == [False, False, False, True]


def test_summary_counts_method_specific_rows():
    summary = summarize_tolerance(sample_detail(), 60.0)
    assert summary["eligible_snr3_rows"] == 3
    assert summary["method_comparable_rows"] == 3
    assert summary["method_rows_over_tolerance"] == 1
    assert summary["registration_rows_over_tolerance"] == 1
    assert summary["assessment"] == "60MAS_METHOD_OR_FRAME_SENSITIVITY_DETECTED"


def test_tolerance_sweep_is_monotonic():
    sweep = tolerance_sweep(sample_detail(), thresholds=[20, 60, 100])
    m = sweep["method_fraction_within"].to_numpy(float)
    r = sweep["registration_fraction_within"].to_numpy(float)
    assert np.all(np.diff(m) >= 0)
    assert np.all(np.diff(r) >= 0)


def test_reference_centroids_are_included(tmp_path: Path):
    ref = tmp_path / "ref.csv"
    pd.DataFrame({
        "candidate_id": [282040],
        "productFilename": ["x.fits"],
        "filter": ["F444W"],
        "detector": ["NRCBLONG"],
        "mjd": [60000.0],
        "snr_err": [10.0],
        "ra_gauss": [150.0],
        "dec_gauss": [2.0],
        "ra_com": [150.0 + 70.0 / (3.6e6 * np.cos(np.deg2rad(2.0)))],
        "dec_com": [2.0],
    }).to_csv(ref, index=False)
    (tmp_path / "candidates").mkdir()
    detail = collect_tolerance_rows(tmp_path, ref)
    assert len(detail) == 1
    assert detail.iloc[0]["source"] == "REFERENCE_282040"
    assert 69.5 < detail.iloc[0]["method_sep_2dg_com_mas"] < 70.5
