import numpy as np
from scripts.audit_miri_third_epoch_1651050 import measure_forced_aperture

def test_detect_centered_synthetic_source():
    y,x=np.indices((61,61))
    sci=10*np.exp(-.5*((x-30)/1.4)**2-.5*((y-30)/1.4)**2)+2
    err=np.ones_like(sci)
    dq=np.zeros_like(sci,dtype=np.uint32)
    res=measure_forced_aperture(sci,err,dq,30,30,"F1000W")
    assert res["status"]=="FORCED_APERTURE_ONLY"
    assert res["snr_conservative"]>5

def test_flagged_core_reduces_aperture_support():
    sci=np.ones((61,61));err=np.ones_like(sci);dq=np.zeros_like(sci,dtype=np.uint32)
    dq[30,30]=1
    out=measure_forced_aperture(sci,err,dq,30,30,"F770W")
    assert out["bad_dq_pixels_aperture"]>=1
