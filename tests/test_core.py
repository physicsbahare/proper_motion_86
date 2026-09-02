import numpy as np
import pandas as pd

from pm86.archive import choose_epoch_pair, normalise_filter
from pm86.astrometry import _fit_affine


def inventory(rows):
    out = pd.DataFrame(rows)
    out["obs_collection"] = "JWST"
    out["instrument_name"] = "NIRCAM/IMAGE"
    out["t_min_numeric"] = out["t_min"]
    out["filter_norm"] = out["filters"].map(normalise_filter)
    return out


def test_filter_normalisation():
    assert normalise_filter("F444W") == "F444W"
    assert normalise_filter("CLEAR;F410M") == "F410M"


def test_same_filter_preferred_over_nearby_cross_filter():
    df = inventory([
        {"filters": "F444W", "t_min": 60000.0},
        {"filters": "F444W", "t_min": 60700.0},
        {"filters": "F410M", "t_min": 60800.0},
    ])
    pair = choose_epoch_pair(df)
    assert pair["status"] == "PAIR_FOUND"
    assert pair["pair_type"] == "SAME_FILTER_JWST"
    assert pair["filter_early"] == "F444W"
    assert pair["filter_late"] == "F444W"


def test_red_cross_filter_beats_blue_same_filter_for_f444_selected_sample():
    df = inventory([
        {"filters": "F115W", "t_min": 60000.0},
        {"filters": "F115W", "t_min": 60700.0},
        {"filters": "F444W", "t_min": 60020.0},
        {"filters": "F356W", "t_min": 60720.0},
    ])
    pair = choose_epoch_pair(df)
    assert pair["status"] == "PAIR_FOUND"
    assert pair["pair_type"] == "CROSS_FILTER_JWST"
    assert {pair["filter_early"], pair["filter_late"]} == {"F444W", "F356W"}


def test_same_night_is_not_independent_epoch():
    df = inventory([
        {"filters": "F444W", "t_min": 60000.0},
        {"filters": "F444W", "t_min": 60000.2},
    ])
    pair = choose_epoch_pair(df)
    assert pair["status"] == "NO_INDEPENDENT_JWST_EPOCH_PAIR"


def test_affine_recovery():
    rng = np.random.default_rng(5)
    src = rng.uniform(-10, 10, size=(30, 2))
    dst = np.column_stack([
        0.002 + 1.001 * src[:, 0] - 0.003 * src[:, 1],
        -0.004 + 0.002 * src[:, 0] + 0.999 * src[:, 1],
    ])
    coef, good, resid, sx, sy = _fit_affine(src, dst)
    assert good.sum() == 30
    assert np.nanmax(np.abs(resid)) < 1e-8
