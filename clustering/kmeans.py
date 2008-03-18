from __future__ import division
from numpy import array, zeros, sqrt, inf
import random
import scipy

__all__ = ['kmeans','repeated_kmeans']

def _euclidean(fmatrix,x):
    return sqrt( ( (fmatrix - x)**2 ).sum(1) )

def _mahalabonis(fmatrix,x,icov):
    diff=(fmatrix-x)
    icov=(icov)
    # The expression below seems to be faster than looping over the elements and summing 
    return sqrt( dot(diff,dot(icov,diff.T)).diagonal() )

def centroid_errors(fmatrix,assignments,centroids):
    errors=[]
    for k in xrange(len(centroids)):
        errors.extend(fmatrix[assignments == k] - centroids[k])
    return array(errors)

def residual_sum_squares(fmatrix,assignments,centroids,distance='euclidean',**kwargs):
    if distance != 'euclidean':
        raise NotImplemented, "residual_sum_squares only implemented for 'euclidean' distance"
    return (centroid_errors(fmatrix,assignments,centroids)**2).sum()

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
    centroids = random.sample(fmatrix,K)
    prev = zeros(N)
    for i in xrange(max_iter):
        dists = array([distfunction(fmatrix,C) for C in centroids])
        assignments = dists.argmin(0)
        if (assignments == prev).all():
            break
        centroids = array([fmatrix[assignments == C].mean(0) for C in xrange(K)])
        prev = assignments
    return assignments, centroids
        
def repeated_kmeans(fmatrix,k,iterations,distance='euclidean',max_iter=1000,**kwargs):
    '''
    assignments,centroids = repeated_kmeans(fmatrix, k, repeats, distance='euclidean',max_iter=1000,**kwargs)

    Runs kmeans repeats times and returns the best result as evaluated according to distance

    @see kmeans
    '''
    if distance != 'euclidean':
        raise NotImplemented, "repeated_kmeans is only implemented for 'euclidean' distance"
    best=+inf
    for i in xrange(iterations):
        A,C=kmeans(fmatrix,k,distance,max_iter=max_iter,**kwargs)
        rss=residual_sum_squares(fmatrix,A,C,distance,**kwargs)
        if rss < best:
            Ab,Cb=A,C
            best=rss
    return Ab,Cb
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
