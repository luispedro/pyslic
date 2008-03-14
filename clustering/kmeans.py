from __future__ import division
from numpy import *
import scipy

__all__ = ['kmeans']

def _euclidean(fmatrix,x):
    return sqrt( ( (fmatrix - x)**2 ).sum(1) )

def _mahalabonis(fmatrix,x,icov):
    diff=(fmatrix-x)
    icov=(icov)
    # The expression below seems to be faster than looping over the elements and summing 
    return sqrt( dot(diff,dot(icov,diff.T)).diagonal() )

def kmeans(fmatrix,K,distance='euclidean',max_iter=1000,**kwargs):
    '''
    assignmens, centroids = kmean(fmatrix, K, distance, icov=None, covmat=None)

    K-Means Clustering

    @param distance can be one of:
        'euclidean' : euclidean distance (default)
        'mahalanobis' : mahalanobis distance.
                This can make use of the following keyword arguments:
                    'icov' (the inverse of the covariance matrix), 
                    'covmat' (the covariance matrix)
                If neither is passed, then the function computes the covariance from the feature matrix
    @param max_iter: Maximum number of iteration
    '''
    if distance == 'euclidean':
        distfunction=_euclidean
    elif distance == 'mahalanobis':
        icov = kwargs.get('icov',None)
        if icov is None:
            covmat=kwargs.get('covmat',None)
            if covmat is None:
                covmat=cov(fmatrix.T)
            icov=linalg.inv(covmat)
        distfunction=lambda f,x: _mahalabonis(f,x,icov)
    else:
        raise 'Distance argument unknown (%s)' % distance

    N,q = fmatrix.shape
    assignments = random.randint(0,K,(N,))
    prev = assignments
    for i in xrange(max_iter):
        centroids = array([fmatrix[assignments == C].mean(0) for C in xrange(K)])
        dists = array([distfunction(fmatrix,C) for C in centroids])
        assignments = dists.argmin(0)
        if (assignments == prev).all():
            break
        prev = assignments
    return assignments, centroids
        
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
