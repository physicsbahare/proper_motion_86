import numpy as np
import pandas as pd
import pytest
from scripts.epsf_1651050_original_cal import isolated_reference_sources,psf_target_fit

def test_avoid_insufficient_stars_instead_of_gaussian_fallback():
    with pytest.raises(RuntimeError,match="cannot construct"):
        isolated_reference_sources(pd.DataFrame(
            [{"x_local":25., "y_local":25.,"snr":5.}]),
            np.ones((101,101)),np.ones((101,101)),np.zeros((101,101),dtype=np.uint32),
            51,51)

def test_reference_star_selection_excludes_close_neighbors_and_bad_dq():
    sci=np.zeros((111,111));err=np.ones_like(sci);dq=np.zeros_like(sci,dtype=np.uint32)
    coords=[(21,21),(79,21),(21,79),(79,79),(55,55),(60,55)]
    r=[]
    for x,y in coords:
        yy,xx=np.indices(sci.shape)
        sci+=100*np.exp(-.5*((xx-x)/1.1)**2-.5*((yy-y)/1.1)**2)
        r.append({"x_local":float(x),"y_local":float(y),"snr":50.,"flux":300.})
    # contaminated close sources get rejected, leaving four isolated refs
    out=isolated_reference_sources(pd.DataFrame(r),sci,err,dq,55,10,min_n=4)
    assert len(out)==4

def test_no_fit_possible_with_all_bad_pixels():
    class DummyPSF:
        pass
    n=31
    with pytest.raises(RuntimeError,match="Too few"):
        psf_target_fit(DummyPSF(),np.ones((n,n)),np.ones((n,n)),
                       np.ones((n,n),dtype=np.uint32),15,15)
