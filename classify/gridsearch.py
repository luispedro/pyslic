from __future__ import division
from nfoldcrossvalidation import *

def _allassignments(D):
    def allassignmentslist(L):
        if len(L) == 0:
            yield []
            return
        key0,vals0=L[0]
        for v in vals0:
            for rest in allassignmentslist(L[1:]):
                yield [(key0,v)]+rest
    for A in allassignmentslist(list(D.items())):
        yield A

def _set_assignment(obj,assignments):
    if not hasattr(obj,'set_option'):
        raise "Don't know how to set options"
    for k,v in assignments:
        obj.set_option(k,v)

class grid_search(object):
    def __init__(self,base,**kwargs):
        self.params = kwargs
        self.base = base
        self.best = None

    def is_multi_class(self):
        return self.base.is_multi_class()

    def train(self,features,labels):
        best_trace=-1
        for assignement in _allassignments(self.params):
            _set_assignment(self.base,assignement)
            S=nfoldcrossvalidation(features,labels,classifier=self.base)
            if S.trace() > best_trace:
                self.best=assignement
                best_trace=S.trace()
        _set_assignment(self.base,self.best)

    def apply(self,features):
        return self.base.apply(features)
