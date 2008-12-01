import pyslic
import numpy
class count_classifier(object):
    def __init__(self):
        self.trained_on = []
        self.tested_on = []

    def train(self,feats,_):
        self.trained_on.extend(feats)

    def apply(self,feats):
        self.tested_on.extend(feats)
        return [0]*len(feats)


def test_nfoldcrossvalidation_use_all():
    features=range(100)
    labels=numpy.zeros(100)
    labels[:30]=1
    C=count_classifier()
    pyslic.classify.nfoldcrossvalidation(features,labels,classifier=C)
    assert len(numpy.unique(C.tested_on))==100
    assert len(C.trained_on)==100*9
    assert len(numpy.unique(C.trained_on))==100

