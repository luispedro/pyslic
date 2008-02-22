from __future__ import division
from numpy import *
from classify import *

__all__=['nfoldcrossvalidation']
def nfoldcrossvalidation(features,labels,nfolds=10,train=learnclassifier,test=applyclassifier):
    '''
    Perform n-fold cross validation

    cmatrix = nfoldcrossvalidation(features, labels, nfolds=10, train=learnclassifier, test=applyclassifier)
    '''
    labels=asarray(labels)
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
        
        model=train(features[trainingset],labels[trainingset])
        prediction=test(features[testingset],model)
        for p, r in zip(prediction,labels[testingset]):
            cmatrix[r,p] += 1

    return cmatrix

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
