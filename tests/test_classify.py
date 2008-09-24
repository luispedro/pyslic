import pyslic
import numpy
L=numpy.r_[[0]*40,[1]*60]
D=numpy.random.rand(100,10)
D[:40] += numpy.random.rand(40,10)**2

def test_nfoldcrossvalidation():
    cmat,Ln = pyslic.classify.nfoldcrossvalidation(D,L)
    assert cmat.shape == (2,2)

def test_stringlabels():
    labelnames=['one','two']
    new_L = map(labelnames.__getitem__,L)
    cmat,Lo = pyslic.classify.nfoldcrossvalidation(D,new_L)
    assert Lo[0] in labelnames
    assert Lo[1] in labelnames
    assert Lo[0] != Lo[1] in labelnames

