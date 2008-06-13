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
import numpy
from numpy import log, pi, array
from numpy.linalg import det, inv
from kmeans import residual_sum_squares, centroid_errors
import scipy

__all__ = ['BIC','AIC','log_likelihood','nr_parameters']

def log_likelihood(fmatrix,assignments,centroids,model='one_variance',covs=None):
    N,q=fmatrix.shape
    k=len(centroids)
    if model == 'one_variance':
        Rss=residual_sum_squares(fmatrix,assignments,centroids)
        #sigma2=Rss/N
        return -N/2.*log(2*pi*Rss/N)-N/2
    elif model == 'diagonal_covariance':
        errors=centroid_errors(fmatrix,assignments,centroids)
        errors *= errors
        errors = errors.sum(1)
        Rss=numpy.zeros(k)
        counts=numpy.zeros(k)
        for i in xrange(fmatrix.shape[0]):
            c=assignments[i]
            Rss[c] += errors[i]
            counts[c] += 1
        sigma2s=Rss/(counts+(counts==0))
        return -N/2.*log(2*pi)-N/2.-1/2.*numpy.sum(counts*numpy.log(sigma2s))
    elif model == 'full_covariance':
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

    raise ValueError, "log_likelihood: cannot handle model '%s'" % model

    
def nr_parameters(fmatrix,k,model='one_variance'):
    N,q=fmatrix.shape
    if model == 'one_variance':
        return k*q+1
    elif model == 'diagonal_covariance':
        return k*(q+1)
    elif model == 'full_covariance':
        return k*+q*q

    raise ValueError, "nr_parameters: cannot handle model '%s'" % model

_BIC = 0
_AIC = 1
def _compute(type,fmatrix,assignements,centroids,model='one_variance',covs=None):
    N,q=fmatrix.shape
    k=len(centroids)
    log_like = log_likelihood(fmatrix,assignements,centroids,model,covs)
    n_param = nr_parameters(fmatrix,k,model)
    if type == _BIC:
        return -2*log_like + n_param * log(N)
    elif type == _AIC:
        return -2*log_like + 2 * n_param

def BIC(fmatrix,assignements,centroids,model='one_variance',covs=None):
    '''
    B = BIC(fmatrix,assignements,centroids,model)

    Compute Bayesian Information Criterion

    model can be one of:
        * 'one_variance': All features share the same variance parameter sigma^2
        * 'full_covariance': Estimate a full covariance matrix or use covs[i] for centroid[i]

    @see AIC
    '''
    return _compute(_BIC,fmatrix,assignements,centroids,model,covs)

def AIC(fmatrix,assignements,centroids,model='one_variance',covs=None):
    '''
    A = AIC(fmatrix,assignements,centroids,model)

    Compute Akaike Information Criterion

    model can be one of:
        * 'one_variance': All features share the same variance parameter sigma^2
        * 'full_covariance': Estimate a full covariance matrix or use covs[i] for centroid[i]

    @see BIC
    '''
    return _compute(_AIC,fmatrix,assignements,centroids,model,covs)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
