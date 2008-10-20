import pyslic
import numpy

def test_scale():
    A=numpy.arange(10)
    scaled=pyslic.preprocess.preprocess._scale(A)
    assert scaled.shape == A.shape
    assert scaled.min() <= 1
    assert scaled.max() >= 254
def test_scale_neg_numbers():
    A=numpy.arange(-10,10)
    scaled=pyslic.preprocess.preprocess._scale(A)
    assert scaled.shape == A.shape
    assert scaled.min() <= 1
    assert scaled.max() >= 254

