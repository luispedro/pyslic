import numpy
import pyslic.features.texture
computecooccurence = pyslic.features.texture.computecooccurence

def test_computecoocurrence():
    A=numpy.zeros((5,5))
    for i in xrange(5): A[i] = i 
    all_dirs=range(4)
    
    assert numpy.all(computecooccurence(A,0)==numpy.eye(4)/4.)
    assert numpy.all(computecooccurence(A,0,remove_zeros=0)==numpy.eye(5)/5.)

    assert numpy.all(computecooccurence(A.T,2)==numpy.eye(4)/4.)
    assert numpy.all(computecooccurence(A.T,2,remove_zeros=0)==numpy.eye(5)/5.)

    
    for dir in all_dirs:
        assert (computecooccurence(A,dir).sum()-1.)**2 < 1.e-8
        assert (computecooccurence(A,0,remove_zeros=0).sum()-1.)**2 < 1.e-8


def test_single_point():
    A=numpy.zeros((5,5))
    A[2,2]=12
    assert not numpy.any(numpy.isnan(pyslic.features.texture.computecooccurence(A,0)))

