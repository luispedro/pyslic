from __future__ import division
from numpy import *
from scipy.ndimage import median_filter, label, center_of_mass, gaussian_filter
from ..imageprocessing.thresholding import otsu, murphy_rc

__all__ = ['labelnuclei','nucleicof']

def labelnuclei(dnaimg,**options):
    '''
    labeled,N = labeled_nuclei(dnaimg, **options)

    N equals the number of nuclei
    labeled is of the same shape as dnaimg and contains the labels of the nuclei
    '''
    sigma=options.get('sigma',10)
    dnaimg=gaussian_filter(dnaimg,sigma)
    thresholding=options.get('thresholding','otsu')
    if thresholding == 'otsu':
        T=otsu(dnaimg)
    elif thresholding == 'murphy_rc':
        T=murphy_rc(dnaimg)
    else:
        raise AttributeError, "Unknown thresholding method (%s)" % thresholding
    return label(dnaimg > T)

def nucleicof(dnaimg,options=None):
    '''
    Returns a set of nuclear centres.
    '''
    labeled,N=labelnuclei(dnaimg)
    cofs=center_of_mass(dnaimg,labeled,range(1,N+1))
    return cofs

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
