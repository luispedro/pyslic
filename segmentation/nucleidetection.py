from __future__ import division
from scipy.ndimage import median_filter, label, center_of_mass
from ..preprocess.thresholding import otsu
def detectnuclei(dnaimg):
    '''
    Returns a set of nuclear centres.
    '''
    dnaimg=median_filter(dnaimg,4)
    T=otsu(dnaimg)
    labeled,N=label(dnaimg > T)
    cofs=[]
    for obj in xrange(1,N+1):
        cofs.append(center_of_mass(labeled == obj))
    return cofs

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
