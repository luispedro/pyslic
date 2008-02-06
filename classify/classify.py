from pymlclassifier import PyMLSVM

__ALL__ = ['learnclassifier','applyclassifier']
C=PyMLSVM

def learnclassifier(featmatrix,y):
    classifier=C()
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
