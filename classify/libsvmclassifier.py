from numpy import array, asanyarray
from svm import *

__ALL__=['libsvmClassifier']
class libsvmClassifier(object):
    def __init__(self):
        pass

    def train(self,features,labels):
        problem=svm_problem(labels,features)
        param=svm_parameter(C = 10,nr_weight = 2,weight_label = [1,0],weight = [10,1])
        param.kernel_type = RBF
        self.model=svm_model(problem,param)
        self.trained=True

    def apply(self,feats):
        assert self.trained
        feats=asanyarray(feats)
        if feats.ndim == 1:
            return self.model.predict(feats)
        else:
            return array([self.apply(feats[i]) for i in xrange(feats.shape[0])])

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
