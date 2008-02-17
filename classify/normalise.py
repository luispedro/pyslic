def normalise(features):
    """
    features = normalise(features)

    Returns a copy of features which has been normalised to zscores 
    """
    mu=features.mean(0)
    features-= mu
    sigma=features.std(0)
    sigma[sigma == 0] = 1
    features /= sigma[sigma != 0]
    return features

class normaliser(object):
    __slots__=['mu','sigma']
    def __init__(self,features=None):
        if features:
            self.train(features)

    def train(self,features,labels):
        self.mu=features.mean(0)
        self.sigma=features.std(0)
        self.sigma[self.sigma == 0.]=1

    def apply(self,features):
        return (features - self.mu)/self.sigma

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
