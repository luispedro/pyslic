import numpy
from pyslic.imageprocessing.basics import fullhistogram

def test_fullhistogram():
    A100 = numpy.arange(100).reshape((10,10))
    assert fullhistogram(A100).shape == (100,)
    assert numpy.all(fullhistogram(A100) == numpy.ones(100))

    A1s = numpy.ones((12,12))
    assert fullhistogram(A1s).shape == (2,)
    assert numpy.all(fullhistogram(A1s) == numpy.array([0,144]))

    A1s[0] = 0
    A1s[1] = 2
    assert fullhistogram(A1s).shape == (3,)
    assert numpy.all(fullhistogram(A1s) == numpy.array([12,120,12]))
    
