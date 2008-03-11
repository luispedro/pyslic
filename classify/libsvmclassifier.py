from __future__ import division
from numpy import array, asanyarray
from classifier import classifier
try:
    from svm import *
except:
    from libsvm.svm import *

__ALL__=['libsvmClassifier']
class libsvmClassifier(classifier):
    def __init__(self):
        classifier.__init__(self)
        self.param = None

    def _dotrain(self,features,labels):
        problem=svm_problem(labels,features)
        if self.param:
            param=self.param
        else:
            param=svm_parameter(kernel_type = RBF, C=1)
        self.model=svm_model(problem,param)

    def _doapply(self,feats):
        return self.model.predict(feats)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
