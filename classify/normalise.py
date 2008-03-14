from __future__ import division
from numpy import *
__all__ = ['zscore','zscore_normalise','interval_normalise','chkfinite']
def zscore(features):
    """
    features = zscore(features)

    Returns a copy of features which has been normalised to zscores 
    """
    mu=features.mean(0)
    sigma=features.std(0)
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
        self.sigma=features.std(0)
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
        

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
