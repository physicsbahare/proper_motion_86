import numpy as np
import pandas as pd
from scripts.check_1651050_gaussian_window import fit_once, summarize_f356w, WINDOW_HALFSIZES

def test_recover_synthetic_gaussian():
    yy,xx=np.indices((69,69))
    sci=.3*np.exp(-.5*((xx-34.3)/.9)**2-.5*((yy-33.1)/.95)**2)+.01
    d={"SCI":sci,"ERR":np.ones_like(sci)*.012,
       "DQ":np.zeros_like(sci,dtype=np.uint64),
       "catalog_x_local":34.0,"catalog_y_local":34.0,
       "pixel_window_x0":0,"pixel_window_y0":0,
       "gaussian_x_local":np.nan,"gaussian_y_local":np.nan}
    fit=fit_once(d,half=5)
    assert fit["optimizer_success"]
    assert not fit["fit_at_bound"]
    assert abs(fit["center_x_local"]-34.3)<.01
    assert abs(fit["center_y_local"]-33.1)<.01

def test_multiple_windows_obligatory():
    assert WINDOW_HALFSIZES==(3,4,5,6,7)
    df=pd.DataFrame([dict(filename="jw05893016008_05101_00001_nrcblong.npz",
                          center_x_local=20+i*.002,center_y_local=19+i*.003,
                          fit_at_bound=False) for i in range(5)])
    dx,dy,r=summarize_f356w(df)
    assert dx<.01 and dy<.02 and r<.03
