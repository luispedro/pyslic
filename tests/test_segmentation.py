import pyslic
import pyslic.segmentation.thresholding
import pyslic.segmentation.watershed
from readmagick import readimg
def __fakeimg():
    dna = readimg('tests/data/dnaimg.jp2')
    img = pyslic.Image()
    img.channeldata['dna'] = dna
    img.loaded = True
    return img
_fakeimg = __fakeimg()
_shape = _fakeimg.channeldata['dna'].shape

def test_thresholding():
    dna = readimg('tests/data/dnaimg.jp2')
    segmented = pyslic.segmentation.thresholding.threshold_segment(dna)
    assert dna.shape == segmented.shape


def test_watershed():
    segmented = pyslic.segmentation.watershed.watershed_segment(_fakeimg)
    assert _shape == segmented.shape

def test_watershed_threshold():
    segmented = pyslic.segmentation.watershed.watershed_segment(_fakeimg,thresholding='mean')
    assert _shape == segmented.shape

