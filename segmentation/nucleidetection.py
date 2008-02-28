from __future__ import division
from numpy import *
from scipy.ndimage import median_filter, label, center_of_mass
from ..imageprocessing.thresholding import otsu

__all__ = ['labelnuclei','nucleicof']

def labelnuclei(dnaimg,options=None):
    '''
    labeled,N = labeled_nuclei(dnaimg,options=None)

    N equals the number of nuclei
    labeled is of the same shape as dnaimg and contains the labels of the nuclei
    '''
    dnaimg=median_filter(dnaimg,4)
    T=otsu(dnaimg)
    return label(dnaimg > T)

def nucleicof(dnaimg,options=None):
    '''
    Returns a set of nuclear centres.
    '''
    labeled,N=labelnuclei(dnaimg)
    cofs=center_of_mass(dnaimg,labeled,range(1,N+1))
    return cofs

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
