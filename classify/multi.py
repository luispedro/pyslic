from __future__ import division
from numpy import *
from classifier import classifier
import random

class one_against_rest(classifier):
    __slots__ = ['_classifiers','_base']

    def __init__(self,base):
        classifier.__init__(self)
        self._classifiers = None
        self._base = base

    def _dotrain(self,features,labels):
        if self._nclasses == 2:
            self._base.train(features,labels)
            return
        #else
        self._classifiers=[]
        for i in xrange(self._nclasses):
            s=self._base.__class__()
            s.train(features,labels==i)
            self._classifiers.append(s)

    def _doapply(self,feats):
        if self._classifiers is None:
            return self._base.apply(feats)
        vals=array([c.apply(feats) for c in self._classifiers])
        idxs, = where(vals == 1)
        if len(idxs) == 1:
            return idxs[0]
        elif len(idxs) == 0:
            return random.randint(0,self._nclasses)
        else:
            # choose at random
            return random.choice(idxs)

class one_against_one(classifier):
    __slots__ = ['_classifiers','_base','_labelnames']

    def __init__(self,base):
        classifier.__init__(self)
        self._classifiers = None
        self._base = base

    def _dotrain(self,features,labels):
        nc=self._nclasses
        if nc == 2:
            self._base.train(features,olabels)
            return
        #else
        self._classifiers=[[None for i in xrange(nc)] for j in xrange(nc)]
        for i in xrange(nc):
            for j in xrange(i+1,nc):
                s=self._base.__class__()
                idxs=(labels == i) | (labels == j)
                s.train(features[idxs],labels[idxs]==i)
                self._classifiers[i][j]=s

    def _doapply(self,feats):
        if self._classifiers is None:
            return self._base.apply(feats)
        nc=self._nclasses
        votes=zeros(nc)
        for i in xrange(nc):
            for j in xrange(i+1,nc):
                c=self._classifiers[i][j].apply(feats)
                if c:
                    votes[i] += 1
                else:
                    votes[j] += 1
        return votes.argmax(0)
    
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
