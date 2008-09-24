# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# 
# Written by Luís Pedro Coelho <lpc@cmu.edu>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# For additional information visit http://murphylab.web.cmu.edu or
# send email to murphy@cmu.edu

from __future__ import division
import numpy
from numpy import *
from bisect import bisect_right
try:
    import ncreduce
    std=ncreduce.std
except:
    std=numpy.std

__all__ = ['zscore','zscore_normalise','interval_normalise','chkfinite','icdf_normalise']

def zscore(features):
    """
    features = zscore(features)

    Returns a copy of features which has been normalised to zscores 
    """
    mu=features.mean(0)
    sigma=std(features,0)
    sigma[sigma == 0] = 1
    return (features - mu) / sigma



class zscore_normalise(object):
    '''
    Normalise to z-scores

    A preprocessor that normalises features to z scores.
    '''
    __slots__=['mu','sigma']
    def __init__(self,features=None):
        if features:
            self.train(features)

    def train(self,features,labels):
        self.mu=features.mean(0)
        self.sigma=std(features,0)
        self.sigma[self.sigma == 0.]=1 # This makes the division a null op.

    def apply(self,features):
        return (features - self.mu)/self.sigma

class interval_normalise(object):
    '''
    Linearly scale to the interval [-1,1] (per libsvm recommendation)

    '''
    __slots__=['mu','range']
    def __init__(self,features=None):
        if features:
            self.train(features)

    def train(self,features,labels):
        self.mu=features.mean(0)
        self.range=features.max(0)-features.min(0)
        self.range[self.range == 0.]=1 # This makes the division a null op.

    def apply(self,features):
        return (features - self.mu)/self.range

    def __getstate__(self):
        return self.mu,self.range
    def __setstate__(self,state):
        self.mu,self.range = state

class chkfinite(object):
    '''
    Fill NaN & Inf values

    Replaces NaN & Inf values with zeros.
    '''
    __slots__ = []
    def __init__(self):
        pass

    def train(self,features,labels):
        pass

    def apply(self,features):
        nans=isnan(features)+isinf(features)
        if nans.any():
            features=features.copy()
            features[nans]=0
        return features
        
def _do_one_icdf_normalise(values):
    values=values.copy()
    N=len(values)
    argsorted=values.argsort()
    corrected=numpy.empty(N)
    try:
        from scipy import weave
        from scipy.weave import converters
        code='''
        for (int i = 0; i != N; ++i) {
            corrected(argsorted(i)) = i/double(N);
        }
        '''
        weave.inline(code,
            ['corrected','argsorted','N'],
            type_converters=converters.blitz)
    except:
        corrected[argsorted]=arange(N)/N
    return corrected

def icdf_normalise(fmatrix):
    '''
    icdf_normalise(fmatrix)

    Normalise each feature so that each value is replaced by the inverse of the feature's cummulative distribution function.
    '''
    if len(fmatrix.shape) == 1:
        return _do_one_icdf_normalise(fmatrix)
    elif len(fmatrix.shape) == 2:
        _,q=fmatrix.shape
        for i in xrange(q):
            fmatrix[:,i]=_do_one_icdf_normalise(fmatrix[:,i])
        return fmatrix
    else:
        raise ValueError, "normalise_to_01 only works with 1 or 2 dimensional arrays!"

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
