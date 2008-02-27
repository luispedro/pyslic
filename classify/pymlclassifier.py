from numpy import array
from PyML import *

__ALL__=['PyMLSVM']
class PyMLSVM(object):
    def __init__(self):
        pass

    def train(self,feats,labels):
        allabels=dict([(L,None) for L in labels])
        if len(allabels) > 2:
            self.s = multi.OneAgainstRest(svm.SVM())
        else:
            self.s = svm.SVM()
        nlabels=[]
        for L in labels:
            nlabels.append(str(L))
        data=datafunc.VectorDataSet(feats,L=nlabels)
        self.s.train(data)
        self.trained=True

    def apply(self,feats):
        assert self.trained
        if feats.ndim == 1:
            data=datafunc.VectorDataSet([feats])
            return array(self.s.classify(data,0))[0]
        else:
            return array([self.apply(feats[i]) for i in xrange(feats.shape[0])])

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
