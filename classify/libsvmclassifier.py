from __future__ import division
from numpy import array, asanyarray
from classifier import classifier, classification_result
try:
    from svm import *
except:
    from libsvm.svm import *

__ALL__=['libsvmClassifier']
class libsvmClassifier(classifier):
    def __init__(self):
        classifier.__init__(self)
        self.param = svm_parameter(kernel_type = RBF)
    
    def set_option(self,optname,value):
        setattr(self.param,optname,value)

    def _dotrain(self,features,labels):
        problem=svm_problem(labels,features)
        self.model=svm_model(problem,self.param)

    def _doapply(self,feats):
        return self.model.predict(feats)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
