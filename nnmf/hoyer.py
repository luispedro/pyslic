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
import numpy as N
import scipy
from ..utils import get_random

__all__ = ['hoyer_sparse_nnmf']

def _norm2(x):
    return (x**2).sum()

def _solve_alpha(s,m,L2):
    sm=s-m
    s2=_norm2(s)
    sm2=_norm2(sm)
    m2=_norm2(m)
    dot=(m*sm).sum()
    alpha = (-dot + N.sqrt(dot**2 - sm2*(m2-L2)))/sm2
    return alpha

def _project(x,L1,L2):
    '''
    Implement projection onto sparse space
    '''
    x = N.asanyarray(x)
    n=len(x)

    s = x + (L1 - x.sum())/n
    Z=N.zeros(n,bool)
    while True:
        m = (~Z) * L1/(n-Z.sum())
        alpha=_solve_alpha(s,m,L2)
        s = m + alpha * (s - m)
        negs = (s < 0)
        if not negs.any():
            return s
        Z |= negs
        s[Z] = 0
        c = (s.sum() - L1)/(n-Z.sum())
        s = s - c*(~Z)

def _L1for(s,x,L2):
    '''
    Solve for L1 in

    s = [ sqrt(n) - L1/sqrt(L2)] / [sqrt(n) - 1]
    '''
    L2=N.sqrt(L2)
    sn=N.sqrt(len(x))
    return L2*s*((sn-1)-1)

def sparse_nnmf(V,r,sparsenessW = None,sparsenessH = None, max_iter=10000,R=None):
    '''
    W,H = hoyer_sparse_nnmf(V,r,sparsenessW = None, sparsenessH = None,max_iter=10000,R=None)

    Implement sparse nonnegative matrix factorisation.

    Reference:
    "Non-negative Matrix Factorisation with Sparseness Constraints"
    by Patrik Hoyer
    in Journal of Machine Learning Research 5 (2004) 1457--1469
    '''
        
    n,m=V.shape
    R=get_random(R)
    mu_W = .15
    mu_H = .15
    eps = 1e-8
    W=R.standard_normal((n,r))**2
    H=R.standard_normal((r,m))**2

    def fixW():
        for i in xrange(r):
            col=W[:,i]
            L2=_norm2(col)
            W[:,i]=_project(col,_L1for(sparsenessW,col,L2),L2)

    def fixH():
        for i in xrange(r):
            row=H[i,:]
            L2=_norm2(row)
            H[i,:]=_project(row,_L1for(sparsenessH,col,L2),L2)

    if sparsenessW is not None: fixW()
    if sparsenessH is not None: fixH()
    for i in xrange(max_iter):
        if sparsenessW is not None:
            W -= mu_W * N.dot(N.dot(W,H)-V,H.T)
            fixW()
        else:
            updateW = N.dot(V,H.T)/(N.dot(W,N.dot(H,H.T))+eps)
            W *= updateW
        if sparsenessH is not None:
            H -= mu_H * N.dot(W.T,N.dot(W,H)-V)
            fixH()
        else:
            updateH = N.dot(W.T,V)/(N.dot(N.dot(W.T,W),H)+eps)
            H *= updateH

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
