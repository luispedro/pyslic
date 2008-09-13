from ..imageprocessing.bbox import bbox
import numpy
def test_bbox():
    img=numpy.zeros((10,10))
    
    a1,b1,a2,b2=bbox(img)
    assert a1==b1
    assert a2 ==b2

    img[4,2]=1
    a1,b1,a2,b2=bbox(img)
    assert a1==4
    assert b1==5
    assert a2==2
    assert b2==3

    img[6,8]=1
    a1,b1,a2,b2=bbox(img)
    assert a1==4
    assert b1==7
    assert a2==2
    assert b2==9
