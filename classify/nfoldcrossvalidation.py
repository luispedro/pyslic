from __future__ import division
import numpy
from numpy import *
from classify import defaultclassifier
from classifier import normaliselabels

__all__=['nfoldcrossvalidation']
def nfoldcrossvalidation(features,labels,nfolds=10,classifier=None):
    '''
    Perform n-fold cross validation

    cmatrix,labelnames = nfoldcrossvalidation(features, labels, nfolds=10, classifier=None)

    cmatrix will be a N x N matrix, where N is the number of classes
    cmatrix[i,j] will be the number of times that an element of class i was classified as class j

    labelnames[i] will correspond to the label name of class i

    @param features: a feature matrix or list of feature vectors
    @param labels: an array of labels, where label[i] is the label corresponding to features[i]

    classifier should implement the train() and apply() methods
    '''
    if classifier is None:
        classifier = defaultclassifier()
    labels,labelnames=normaliselabels(labels)

    features = numpy.asanyarray(features)

    classcounts={}
    for L in labels:
        classcounts[L] = classcounts.get(L,0) + 1
    
    nclasses=len(classcounts)
    cmatrix=zeros((nclasses,nclasses))
    for fold in xrange(nfolds):
        testingset=zeros(len(labels),bool)

        for L,C in classcounts.items():
            idxs,=where(labels==L)
            N=len(idxs)
            perfold=N/nfolds
            start=floor(perfold*fold)
            end=floor(perfold*(fold+1))
            idxs=idxs[start:end]
            testingset[idxs]=True
        trainingset= ~testingset
        
        classifier.train(features[trainingset],labels[trainingset])
        prediction=classifier.apply(features[testingset])
        for p, r in zip(prediction,labels[testingset]):
            cmatrix[r,p] += 1

    return cmatrix, labelnames

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
