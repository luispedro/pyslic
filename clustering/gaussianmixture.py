from __future__ import division
import numpy
from numpy import log, pi, array
from numpy.linalg import det, inv
from kmeans import residual_sum_squares
import scipy

__all__ = ['BIC','AIC']


def logP_onediagonalcovariance(fmatrix,assignments,centroids):
    N,q=fmatrix.shape
    Rss=residual_sum_squares(fmatrix,assignments,centroids)
    sigma2=Rss/(N-q)
    return -.5*N *(q*log(sigma2) + log(2*pi)) + Rss/2./sigma2

def logP_fullcovariance(fmatrix,assignments,centroids,covs=None,**kwargs):
    N,q=fmatrix.shape
    res=-N*q/2.*log(2*pi)

    for k in xrange(len(centroids)):
        diff=(fmatrix[assignments == k] - centroids[k])
        if covs is None:
            covm=cov(diff.T)
        else:
            covm=covs[k]
        if covm.shape == ():
            covm=mat([[covm]])
        icov=mat(inv(covm))
        diff=mat(diff)
        Nk = diff.shape[0]
        res += -Nk/2.*log(det(covm)) + \
             -.5 * (diff * icov * diff.T).diagonal().sum() 
    return res
    
def nrparameters_fullcovariance(fmatrix,k):
    N,q=fmatrix.shape
    return k*(N+q*q)

def nrparameters_diagonalcovariance(fmatrix,k):
    N,q=fmatrix.shape
    return k*(q+N)

def nrparameters_onediagonalcovariance(fmatrix,k):
    N,q=fmatrix.shape
    return k*(1+N)

def BIC_onediagonalcovariance(fmatrix,assignements,centroids):
    N,q=fmatrix.shape
    k=len(centroids)
    L=logP_onediagonalcovariance(fmatrix,assignements,centroids)
    nrP=nrparameters_onediagonalcovariance(fmatrix,k)
    return -2*L+nrP*log(N)

def AIC_onediagonalcovariance(fmatrix,assignements,centroids):
    N,q=fmatrix.shape
    k=len(centroids)
    L=logP_onediagonalcovariance(fmatrix,assignements,centroids)
    nrP=nrparameters_onediagonalcovariance(fmatrix,k)
    return -2*L+2*nrP

def BIC(fmatrix,assignements,centroids,model='onevariance'):
    '''
    B = BIC(fmatrix,assignements,centroids,model)

    Compute Bayesian Information Criterion

    model can be one of:
        * 'onevariance': All features share the same variance parameter sigma2

    @see AIC
    '''
    if model == 'onevariance':
        return BIC_onediagonalcovariance(fmatrix,assignements,centroids)
    else:
        raise NotImplementedError, "BIC only supports model 'onevariance'"

def AIC(fmatrix,assignements,centroids,model='onevariance'):
    '''
    A = AIC(fmatrix,assignements,centroids,model)

    Compute Akaike Information Criterion

    model can be one of:
        * 'onevariance': All features share the same variance parameter sigma2

    @see BIC
    '''
    if model == 'onevariance':
        return AIC_onediagonalcovariance(fmatrix,assignements,centroids)
    else:
        raise NotImplementedError, "BIC only supports model 'onevariance'"

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
