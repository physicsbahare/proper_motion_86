import numpy as np
import pytest
from pm86.injection_trial_summary import summarize_trials


def test_null_tail_and_failed_recoveries_are_separate():
    stationary = np.array([[0., 0.], [3., 4.], [np.nan, 0.], [6., 8.]])
    moving = np.array([[10., 0.], [12., -2.], [np.nan, np.nan]])
    s = summarize_trials(stationary, moving, [10., 0.], [20., 0.])
    assert s["n_null_attempted"] == 4
    assert s["n_null_recovered"] == 3
    assert s["null_recovery_fraction"] == pytest.approx(.75)
    assert s["moving_recovery_fraction"] == pytest.approx(2/3)
    assert s["moving_recovery_bias_east_north_masyr"] == pytest.approx([1., -1.])
    assert s["stationary_amplitude_tail"]["p_empirical"] == pytest.approx(1/4)
    assert s["valid_for_publication"] is False


def test_reject_malformed_or_unrecovered_trials():
    with pytest.raises(ValueError):
        summarize_trials(np.zeros((2, 2)), np.zeros((2, 3)), [0, 0], [0, 0])
    with pytest.raises(ValueError):
        summarize_trials(np.full((2, 2), np.nan), np.zeros((2, 2)), [0, 0], [0, 0])


def test_500_trials_do_not_validate_five_sigma():
    s = summarize_trials(np.zeros((500, 2)), np.zeros((500, 2)),
                         [0, 0], [100, 100])
    assert s["stationary_amplitude_tail"]["minimum_resolvable_p"] == pytest.approx(1/501)
    assert s["stationary_amplitude_tail"]["p_empirical"] == pytest.approx(1/501)
