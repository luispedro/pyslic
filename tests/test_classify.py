import pyslic
import numpy
L=numpy.r_[[0]*40,[1]*60]
D=numpy.random.rand(100,10)
D[:40] += numpy.random.rand(40,10)**2

def slow(f):
    f.slow=True
    return f

@slow
def test_nfoldcrossvalidation():
    cmat,Ln = pyslic.classify.nfoldcrossvalidation(D,L)
    assert cmat.shape == (2,2)

@slow
def test_stringlabels():
    labelnames=['one','two']
    new_L = map(labelnames.__getitem__,L)
    cmat,Lo = pyslic.classify.nfoldcrossvalidation(D,new_L)
    assert Lo[0] in labelnames
    assert Lo[1] in labelnames
    assert Lo[0] != Lo[1] in labelnames


def test_nnclassifier():
    labels=[0,1]
    data=[[0.,0.],[1.,1.]]
    C=pyslic.classify.NNClassifier()
    C.train(data,labels)
    assert C.apply(data[0]) == 0
    assert C.apply(data[1]) == 1
    assert C.apply([.01,.01]) == 0
    assert C.apply([.99,.99]) == 1
    assert C.apply([100,100]) == 1
    assert C.apply([-100,-100]) == 0
    assert C.apply([.9,.9]) == 1
    middle = C.apply([.5,.5])
    assert (middle == 0) or (middle == 1)

