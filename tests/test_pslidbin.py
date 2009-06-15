import pyslic.features.pslidbinformat
import numpy as np
import os
def test_plsidbin():
    f = np.arange(12).reshape((1,12))
    slf = ['SLF%i' for i in xrange(12)]
    rslf = ['RSLF%i' for i in xrange(12)]
    names = ['feat_%i' for i in xrange(12)]
    imageurls = ['image']
    maskurls = ['mask']
    settype = 2
    channels = [1 for i in xrange(12)]
    pyslic.features.pslidbinformat.writepslidbin('test.bin',f,rslf,slf,names,imageurls,maskurls,settype,channels)
    f2, r2, s2, n2, i2, m2 = pyslic.features.pslidbinformat.readpslidbin('test.bin')
    assert r2 == rslf
    assert np.all(f == f2)
    assert slf == s2
    assert names == n2
    os.unlink('test.bin')
