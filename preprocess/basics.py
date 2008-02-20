from scipy.ndimage import histogram, convolve
from numpy import asarray, ones, uint32

__all__ = ['fullhistogram','majority_filter']
def fullhistogram(img):
    """
    Return a histogram containing all the values in the image
    """
    maxt=img.max()
    return histogram(img,0,maxt,maxt+1)

def majority_filter(bwimg, N = 3):
    """
    Implement a majority filter of size 0.

    bwimg is interpreted as a binary image (non-zero pixels correspond to one/True).
    """
    assert N < 2 ** 16
    bwimg=asarray(bwimg > 0,uint32)
    F=ones((N,N))
    bwimg=convolve(bwimg,F)
    return bwimg > (N**2/2)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
