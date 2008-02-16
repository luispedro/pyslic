def normalise(features):
    """
    features = normalise(features)

    Returns a copy of features which has been normalised to zscores and
    with constant features removed    
    """
    sigma=features.std(0)
    features=features[:,sigma != 0]
    mu=features.mean(0)
    features-= mu
    features /= sigma[sigma != 0]
    return features

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
