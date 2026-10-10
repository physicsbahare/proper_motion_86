import numpy as np
import pytest
from scripts.jwst_parallax_1651050 import basis_icrs,vectors_to_parallax_factors,conditional_pm

def test_parallax_vectors_have_correct_basis_and_sign():
    east,north=basis_icrs(0,0)
    assert east == pytest.approx([0,1,0],abs=1e-14)
    assert north == pytest.approx([0,0,1],abs=1e-14)
    got=vectors_to_parallax_factors([[0,1,0],[0,0,1]],0,0)
    assert np.allclose(got, [[-1,0],[0,-1]],atol=1e-14)

def test_conditional_motion_calculation():
    r=conditional_pm(np.array([50.,10.]),np.array([5.,5.]),np.array([.1,-.2]),2.,10.)
    assert r["conditional_pm_east_masyr"]==pytest.approx(45)
    assert r["conditional_pm_north_masyr"]==pytest.approx(20)
    assert 0 < r["residual_p_2d_zero_motion"] < 1

def test_same_epoch_factors_cannot_create_spurious_parallax():
    f=vectors_to_parallax_factors([[.7,.3,-.1]]*2)
    assert f[1]-f[0] == pytest.approx([0.,0.])
