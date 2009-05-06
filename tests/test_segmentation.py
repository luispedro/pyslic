import pyslic.segmentation.thresholding
from readmagick import readimg
def test_thresholding():
    dna = readimg('tests/data/dnaimg.jp2')
    segmented = pyslic.segmentation.thresholding.threshold_segment(dna)
    assert dna.shape == segmented.shape

