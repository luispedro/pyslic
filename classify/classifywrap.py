__all__ = ['pretransformclassifier']
class pretransformclassifier(object):
    def __init__(self,transformer,classifier):
        self.transformer=transformer
        self.classifier=classifier

    def train(self,features,labels):
        self.transformer.train(features,labels)
        features=self.transformer.apply(features)
        self.classifier.train(features,labels)
    
    def apply(self,features):
        features=self.transformer.apply(features)
        return self.classifier.apply(features)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
