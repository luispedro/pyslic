import numpy
import numpy as np
import pyslic
import pyslic.imageprocessing.thresholding
from os.path import dirname
img = pyslic.readimg(dirname(__file__) + '/data/protimg.jp2')


def test_threshold():
    assert pyslic.imageprocessing.thresholding.threshold(img,None)  == -1
    assert pyslic.imageprocessing.thresholding.threshold(img,'otsu') == pyslic.imageprocessing.thresholding.otsu(img)
    assert pyslic.imageprocessing.thresholding.threshold(img,'rc') == pyslic.imageprocessing.thresholding.rc(img)
    assert pyslic.imageprocessing.thresholding.threshold(img,pyslic.imageprocessing.thresholding.rc) == pyslic.imageprocessing.thresholding.rc(img)
    assert pyslic.imageprocessing.thresholding.threshold(img,pyslic.imageprocessing.thresholding.rc) == pyslic.imageprocessing.thresholding.rc(img)
    assert pyslic.imageprocessing.thresholding.threshold(img,12) == 12
    assert pyslic.imageprocessing.thresholding.threshold(img,12.) == 12.
    assert pyslic.imageprocessing.thresholding.threshold(img,np.uint16(12)) == 12
    
