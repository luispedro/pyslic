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
from numpy import *
from numpy.linalg import det
from classifier import normaliselabels, classifier
import scipy.stats
import warnings

_TOLERANCE = 0
_SIGNIFICANCE_IN = .15
_SIGNIFICANCE_OUT = .15

def _sweep(A,k,flag):
    N,_=A.shape
    Akk=A[k,k]
    B=zeros_like(A)
    try:
        from scipy import weave
        from scipy.weave import converters
        k=int(k)
        code = '''
#line 32 "featureselection.py"
        for (int i = 0; i != N; ++i) {
            for (int j = 0; j != N; ++j) {
                if (i == k) {
                    if (j == k) {
                        B(i,j) =  - 1./A(k,k);
                    } else {
                        B(i,j)=flag*A(i,j)/A(k,k);
                    }
                } else if (j == k) {
                    B(i,j)=flag*A(i,j)/A(k,k);
                } else { 
                    B(i,j)=A(i,j) - A(i,k)*A(k,j)/A(k,k);
                }
            }
        }
        '''
        weave.inline(
                code,
                ['A','B','k','Akk','flag','N'],
                type_converters=converters.blitz)
    except:
        for i in xrange(N):
            for j in xrange(N):
                if i == k:
                    if j == k:
                        B[i,j] =  - 1./Akk
                    else:
                        B[i,j]=flag*A[i,j]/Akk
                elif j == k:
                    B[i,j]=flag*A[i,j]/Akk
                else:
                    B[i,j]=A[i,j] - A[i,k]*A[k,j]/Akk
    return B

def sda(features,labels):
    '''
    features_idx = sda(features,labels)

    Perform Stepwise Discriminant Analysis for feature selection

    Pre filter the feature matrix to remove linearly dependent features
    before calling this function. Behaviour is undefined otherwise.

    This implements the algorithm described in 
    Jennrich, R.I. (1977), "Stepwise Regression" & "Stepwise Discriminant Analysis,"
    both in Statistical Methods for Digital Computers, eds. 
    K. Enslein, A. Ralston, and H. Wilf, New York; John Wiley & Sons, Inc.
    '''


    N, m = features.shape
    labels,labelsu = normaliselabels(labels)
    q=len(labelsu)

    mus=array([features[labels==i,:].mean(0) for i in xrange(q)])
    mu=features.mean(0)
    
    W=zeros((m,m))
    T=zeros((m,m))
    try:
        from scipy import weave
        from scipy.weave import converters
        code='''
#line 68 "featureselection.py"
        for (int i = 0; i != m; ++i) {
            for (int j = 0; j != m; ++j) {
                for (int n = 0; n != N; ++n) {
                    int g=labels(n);
                    W(i,j) += (features(n,i)-mus(g,i))*(features(n,j)-mus(g,j));
                    T(i,j) += (features(n,i)-mu(i))*(features(n,j)-mu(j));
                }
            }
        }
        '''
        weave.inline(
                code,
                ['N','m','W','T','features','mu','mus','labels'],
                type_converters=converters.blitz)
    except:
        warnings.warn('scipy.weave failed. Resorting to (slow) Python code')
        for i in xrange(m):
            for j in xrange(m):
                for n in xrange(N):
                    g=labels[n]
                    W[i,j] += (features[n,i]-mus[g,i])*(features[n,j]-mus[g,j])
                    T[i,j] += (features[n,i]-mu[i])*(features[n,j]-mu[j])
    ignoreidx = ( W.diagonal() == 0 )
    if ignoreidx.any():
        idxs, = where(~ignoreidx)
        F=sda(features[:,~ignoreidx],labels)
        return idxs[F]
    output=[]
    D=W.diagonal()
    df1 = q-1
    while True:
        V = W.diagonal()/T.diagonal() 
        W_d = W.diagonal()
        V_neg = (W_d < 0)
        p=V_neg.sum()
        df2 = N-p-q+1
        if V_neg.any(): 
            V_m = V[V_neg].min()
            k,=where(V == V_m)
            k=k[0]
            Fremove = (N-p-q+1)/(q-1)*(V_m-1)
            PrF = 1 - scipy.stats.f.cdf(Fremove,df1,df2)
            if PrF > _SIGNIFICANCE_OUT:
                #print 'removing ',k, 'V(k)', 1./V_m, 'Fremove', Fremove
                W=_sweep(W,k,1)
                T=_sweep(T,k,1)
                continue
        ks = ( (W_d / D) > _TOLERANCE)
        if ks.any():
            V_m=V[ks].min()
            k,=where(V==V_m)
            k=k[0]
            Fenter = (N-p-q)/(q-1) * (1-V_m)/V_m
            PrF = 1 - scipy.stats.f.cdf(Fenter,df1,df2)
            if PrF < _SIGNIFICANCE_IN:
                #print 'adding ',k, 'V(k)', 1./V_m, 'Fenter', Fenter
                W=_sweep(W,k,-1)
                T=_sweep(T,k,-1)
                if PrF < .0001:
                    output.append((Fenter,k))
                continue
        break

    output.sort(reverse=True)
    return array([x[1] for x in output])

def repeatedfeatures(featmatrix):
    '''
    idxs = repeatedfeatures(feature_matrix)

    Returns a set of indices corresponding to repeated features in featmatrix.

    @see remove_repeated_features
    '''
    featsig={}
    results=[]
    feat_sum=(featmatrix**2).sum(0)
    ordering=feat_sum.argsort()
    i = 0
    #while i < len(ordering) and feat_sum[ordering[i]] == 0:
    #    results.append(ordering[i])
    #    i += 1

    while i < len(ordering) - 1:
        if feat_sum[ordering[i]] == feat_sum[ordering[i+1]]:
            if (featmatrix[:,ordering[i]] == featmatrix[:,ordering[i+1]]).all():
                results.append(ordering[i+1])
        i += 1
    return array(results)


def _rank(A,tol=1e-8):
    s = linalg.svd(A,compute_uv=0)
    return (s > tol).sum()

def linearindependentfeatures(featmatrix):
    '''
    Returns the indices of a set of linearly independent features (columns).

    indices = linearindependentfeatures(features)
    '''
    independent=[]
    R=_rank(featmatrix)
    i=0
    offset=0
    while i < featmatrix.shape[1]:
        R_=_rank(delete(featmatrix,i,1))
        if R_ == R:
            featmatrix=delete(featmatrix,i,1)
            offset+=1
        else:
            independent.append(i+offset)
            i += 1
    return array(independent)
        
class sda_filter(object):
    __slots__ = ['idxs']

    def is_multi_class(self):
        return True

    def __init__(self):
        pass

    def train(self,features,labels):
        self.idxs=sda(features,labels)
        if len(self.idxs) == 0:
            self.idxs = [0]

    def apply(self,features):
        return features[:,self.idxs]

    def __getstate__(self):
        return (self.idxs,)
    def __setstate__(self,state):
        self.idxs, = state


class remove_repeated_features(object):
    __slots__ = ['repeats']

    def is_multi_class(self):
        return True
    def __init__(self):
        pass
    def train(self,features,labels):
        self.repeats=repeatedfeatures(features)

    def apply(self,features):
        return delete(features,self.repeats,1)
    def __getstate__(self):
        return (self.repeats,)
    def __setstate__(self,state):
        self.repeats, = state

        
class remove_linear_dependent_features(object):
    __slots__ = ['independent']

    def is_multi_class(self):
        return True
    def __init__(self):
        pass
    def train(self,features,labels):
        self.independent=linearindependentfeatures(features)

    def apply(self,features):
        return features[:,self.independent]
    def __getstate__(self):
        return (self.independent,)
    def __setstate__(self,state):
        self.independent, = state


