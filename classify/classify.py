from __future__ import division

__ALL__ = ['learnclassifier','applyclassifier','defaultclassifier']

def defaultclassifier():
    '''
    C = defaultclassifier()

    Returns the default classifier
    '''
    from normalise import zscore_normalise,chkfinite,interval_normalise
    from gridsearch import parameter_grid
    from classifywrap import pretransformclassifier
    from libsvmclassifier import libsvmClassifier
    from pymlclassifier import PyMLSVM
    from featureselection import sda_filter
    from numpy import arange
    return pretransformclassifier([chkfinite(),interval_normalise(),sda_filter()],fit_parameters(libsvmClassifier(),C=2.**arange(-4,10), gamma=2.**arange(-4,10)))

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
