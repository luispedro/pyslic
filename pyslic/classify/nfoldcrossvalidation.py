# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# 
# Written by Luís Pedro Coelho <lpc@cmu.edu>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# For additional information visit http://murphylab.web.cmu.edu or
# send email to murphy@cmu.edu

from __future__ import division
import numpy
import numpy as np
from numpy import *
from classify import defaultclassifier
from classifier import normaliselabels

__all__=['nfoldcrossvalidation']
def nfoldcrossvalidation(features, labels, nfolds=None, classifier=None, return_predictions=False, max_train=None):
    '''
    Perform n-fold cross validation

    cmatrix,labelnames = nfoldcrossvalidation(features, labels, nfolds=10, classifier=None, return_predictions=False)
    cmatrix,labelnames,predictions = nfoldcrossvalidation(features, labels, nfolds=10, classifier=None, return_predictions=False)

    cmatrix will be a N x N matrix, where N is the number of classes
    cmatrix[i,j] will be the number of times that an element of class i was classified as class j

    labelnames[i] will correspond to the label name of class i

    @param features: a feature matrix or list of feature vectors
    @param labels: an array of labels, where label[i] is the label corresponding to features[i]
    @param nfolds: Nr of folds.
    @param return_predictions: whether to return predictions

    classifier should implement the train() and apply() methods
    '''
    assert len(features) == len(labels), 'nfoldcrossvalidation: len(features) should match len(labels)'
    if classifier is None:
        classifier = defaultclassifier()
    labels,labelnames = normaliselabels(labels)
    predictions = labels*0-1

    features = numpy.asanyarray(features)

    classcounts={}
    for L in labels:
        classcounts[L] = classcounts.get(L,0) + 1

    min_class_count = min(classcounts.values())
    if nfolds is None:
        nfolds = min(10,min_class_count)
    elif min_class_count < nfolds:
        from warnings import warn
        warn('pyslic.classify.nfoldcrossvalidation: Reducing the nr. of folds to %s (minimum class size).' % min_class_count)
        nfolds = min_class_count
    
    nclasses = len(classcounts)
    cmatrix = zeros((nclasses,nclasses))
    for fold in xrange(nfolds):
        testingset=zeros(len(labels),bool)

        for L,C in classcounts.items():
            idxs, = where(labels==L)
            N = len(idxs)
            perfold = N/nfolds
            start = floor(perfold*fold)
            end = floor(perfold*(fold+1))
            idxs = idxs[start:end]
            testingset[idxs] = True
        trainingset= ~testingset
        if max_train is not None:
            for c in xrange(nclasses):
                if (labels[trainingset] == c).sum() > max_train:
                    idxs, = np.where( (labels == L) & trainingset )
                    trainingset[idxs[max_train:]] = False
        
        classifier.train(features[trainingset],labels[trainingset])
        prediction=classifier.apply(features[testingset])
        predictions[testingset]=prediction
        for p, r in zip(prediction,labels[testingset]):
            cmatrix[r,p] += 1

    if return_predictions:
        return cmatrix, labelnames, predictions
    return cmatrix, labelnames

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
