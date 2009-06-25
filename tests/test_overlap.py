import pyslic.features
import numpy as np
def test_overlap_simple():
    p = np.zeros((10,10),np.uint8)
    d = np.zeros((10,10),np.uint8)
    p[2:8,2:8] = np.arange(36).reshape((6,6))
    pp = p * (p > 15)
    d[2:8,2:8] = np.arange(36).reshape((6,6))//2
    pd = d * (d > 3)
    vs = pyslic.features.overlap.overlapfeatures(p,d,pp,pd)
    assert vs[0] == vs[1]
    assert np.all(vs >= 0)
    assert np.all(vs[:8] <= 1)

