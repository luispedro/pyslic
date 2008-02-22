__all__ = ['pretransformclassifier','concattransformers']

class concattransformers(object):
    def __init__(self,transformers):
        self.transformers = transformers

    def train(self,features,labels):
        for t in self.transformers:
            t.train(features,labels)
            features=t.apply(features)

    def apply(self,features):
        for t in self.transformers:
            features=t.apply(features)
        return features

class pretransformclassifier(object):
    def __init__(self,transformer,classifier):
        if type(transformer) == list:
            transformer = concattransformers(transformer)
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
