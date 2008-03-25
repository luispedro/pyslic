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
from numpy import dot
from ..utils import get_random
import scipy

__all__ = ['nnmf']

def nnmf(V,r,cost='norm2',max_iter=1e4,tol=1e-8,R=None):
    '''
    A,S = nnmf(X,r,cost='norm2',tol=1e-8,R=None)

    Implement Lee & Seung's algorithm

    @param cost can be one of:
        'norm2' : minimise || X - AS ||_2
        'i-div' : minimise D(X||AS), where D is I-divergence (generalisation of K-L divergence)

    @param max_iter Maximum number of iterations
    @param tol tolerance Threshold for early exit (when the update factor is with tol of 1., the function exits)
    @param R random seed @see get_random

    Reference:
    "Algorithms for Non-negative Matrix Factorization"
    by Daniel D Lee, Sebastian H Seung
    (available at http://citeseer.ist.psu.edu/lee01algorithms.html)
    '''
    # Nomenclature in the function follows lee & seung, while outside nomenclature follows 
    eps = 1e-8
    n,m = V.shape
    R=get_random(R)
    W=R.standard_normal((n,r))**2
    H=R.standard_normal((r,m))**2
    for i in xrange(max_iter):
        if cost == 'norm2':
            updateH = dot(W.T,V)/(dot(dot(W.T,W),H)+eps)
            H *= updateH
            updateW = dot(V,H.T)/(dot(W,dot(H,H.T))+eps)
            W *= updateW
        elif cost == 'i-div':
            raise NotImplementedError,'I-Div not implemented in lee_seung.nnmf'
        if True or (i % 10) == 0:
            max_update = max(updateW.max(),updateH.max())
            if abs(1.-max_update) < tol:
                break
    return W,H,D

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
