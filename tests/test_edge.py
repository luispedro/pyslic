import pyslic
from pyslic.imageprocessing.edge import sobeledge
import numpy

def test_sobel():
    A=numpy.zeros((100,100))
    assert sobeledge(A).shape == A.shape
    A=numpy.zeros((15,100))
    assert sobeledge(A).shape == A.shape
    assert sobeledge(A).sum() == 0

