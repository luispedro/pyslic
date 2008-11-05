from os.path import exists
from os import unlink
import numpy
from pyslic import readimg, writeimg
from nose.tools import with_setup

def write_readback(C,fname):
    C=C.copy()
    writeimg(C,fname)
    C2=readimg(fname)
    return C, C2


_PNG = '/tmp/pyslic_test_file.png'
_TIFF = '/tmp/pyslic_test_file.tiff'
def delete_tmps():
    if exists(_PNG): unlink(_PNG)
    if exists(_TIFF): unlink(_TIFF)

@with_setup(None,delete_tmps)
def test_uint8_color():
    A=numpy.array([
        [[0,1,2],[1,2,3],[2,3,4]],
        [[5,5,5],[6,6,6],[7,7,7]]
        ],
        numpy.uint8)
    B=numpy.array([
        [0,1,2,3,4],
        [5,5,5,5,5]
        ],
        numpy.uint8)
    C=numpy.array([
        [[0,1,2],[1,2,3],[2,3,4]],
        [[5,5,5],[6,6,6],[7,7,7]],
        [[257,257,257],[256,256,256],[707,707,707]]
        ],
        numpy.uint16)
    D=numpy.array([
        [0,1,2,3,4],
        [5,5,5,5,5],
        [258,258,258,258,258]
        ],
        numpy.uint16)
    def test_one(C,fname,correct):
        writeimg(C,fname)    
        C,C2=write_readback(C,fname)
        if correct: C2 //= 257
        assert C.shape == C2.shape
        assert numpy.all(C == C2)
    yield test_one, A, _PNG, True
    yield test_one, B, _PNG, True
    yield test_one, C, _TIFF, False
    yield test_one, D, _TIFF, False

