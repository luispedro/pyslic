from __future__ import division
from numpy import *
def normaliselabels(labels):
    labelnames={}
    normalised=[]
    names=[]
    N=0
    for L in labels:
        nr=labelnames.get(L,N)
        if nr == N:
            labelnames[L]=N
            names.append(L)
            N += 1
        normalised.append(nr) 
    return array(normalised),names

class classifier(object):
    __slots__ = ['_trained','_labelnames','_nclasses']
    def __init__(self):
        self._trained = False

    def train(self,features,labels):
        nlabels,self._labelnames=normaliselabels(labels)
        self._nclasses = nlabels[-1]
        features=asanyarray(features)
        self._dotrain(features,nlabels)
        self._trained = True
    
    def apply(self,features):
        assert self._trained
        features=asanyarray(features)
        if features.ndim == 1:
            return self._labelnames[int(self._doapply(features))]
        else:
            return array([self._labelnames[int(self._doapply(features[i]))] for i in xrange(features.shape[0])])

