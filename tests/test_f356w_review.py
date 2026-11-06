"""Offline tests for the publication F356W review geometry."""
import numpy as np
from scripts.inspect_1651050_f356w import window2d, centroid_distance_pixels, local_peaks

def test_window_preserves_original_pixel_mapping():
    a=np.arange(101*101).reshape(101,101)
    s,x0,y0=window2d(a,50.1,49.8,12)
    assert s.shape==(25,25)
    assert s[12,12]==a[50,50]
    assert (x0,y0)==(38,38)

def test_centroid_disagreement_reported_not_silently_discarded():
    d={"x_2dg":200.17375,"y_2dg":198.941867,
       "x_com":202.217118,"y_com":198.179901}
    assert abs(centroid_distance_pixels(d)-2.1809)<.003
    assert np.isnan(centroid_distance_pixels({"x_2dg":np.nan,"y_2dg":1.0}))

def test_peak_finder_excludes_main_source_and_keeps_neighbor():
    a=np.zeros((31,31))
    a[15,15]=10
    a[15,21]=7
    a[5,5]=9
    peaks=local_peaks(a,15,15,radius=10,threshold=4)
    assert len(peaks)==1
    assert peaks[0][:2]==(21.,15.)

def test_crop_at_edge_is_bounded():
    a=np.ones((10,10))
    z,x0,y0=window2d(a,0,0,5)
    assert z.shape==(6,6)
    assert x0==0 and y0==0
