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

__ALL__ = ['learnclassifier','fastclassifier','applyclassifier','defaultclassifier']

def defaultclassifier():
    '''
    C = defaultclassifier()

    Returns the default classifier
    '''
    from normalise import zscore_normalise,chkfinite,interval_normalise
    from gridsearch import gridsearch
    from classifywrap import pretransformclassifier
    from libsvmclassifier import libsvmClassifier
    from featureselection import sda_filter, remove_repeated_features, remove_linear_dependent_features
    from numpy import arange
    return pretransformclassifier(
                    [chkfinite(),interval_normalise(),remove_linear_dependent_features(),sda_filter()], # pre-process
                    gridsearch(libsvmClassifier(),C=2.**arange(-4,10), gamma=2.**arange(-7,2)) # classify
                    )

def fastclassifier():
    '''
    C = fastclassifier()

    Returns the classifier, which is faster than the defaultclassifier()

    Currently, this means that this classifier does not perform a grid search for the best libsvm parameters
    '''
    from normalise import zscore_normalise,chkfinite,interval_normalise
    from gridsearch import gridsearch
    from classifywrap import pretransformclassifier
    from libsvmclassifier import libsvmClassifier
    from pymlclassifier import PyMLSVM
    from featureselection import sda_filter, remove_repeated_features, remove_linear_dependent_features
    from numpy import arange
    return pretransformclassifier([chkfinite(),interval_normalise(),remove_linear_dependent_features(),sda_filter()],libsvmClassifier())

def learnclassifier(featmatrix,y):
    """
    classifier=learnclassifier(featmatrix,y)

    Learn a classifier for the problem y[i]=f(featmatrix[i])
    """
    classifier=defaultclassifier()
    classifier.train(featmatrix,y)
    return classifier

def applyclassifier(feats,model):
    """
    y = applyclassifier(feats,models)

    output a label for feats according to the model

    feats can be either a feature array or a sequence of feature arrays
    """
    return model.apply(feats)



# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
