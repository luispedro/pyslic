from __future__ import division
from numpy import *

__all__ = ['accuracy','waccuracy']

def accuracy(cmatrix):
    '''
    acc = accuracy(cmatrix)

    Accuracy of cmatrix
    '''
    return cmatrix.trace()/cmatrix.sum()

def waccuracy(cmatrix):
    '''
    wacc = waccuracy(cmatrix)

    Weighted accuracy of cmatrix
    '''
    return (cmatrix.diagonal() / cmatrix.sum(1)).mean()
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
