# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# 
# Written by Luis Pedro Coelho <lpc@cmu.edu>
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
from ..imageprocessing import thresholding
from ..utils import get_random
from collections import defaultdict
from scipy import ndimage
def _unhaar(X):
    w,h=X.shape
    X2=numpy.empty((w*2,h*2))
    X2[::2,::2]=X
    X2[1::2,::2]=X
    X2[::2,1::2]=X
    X2[1::2,1::2]=X
    return X2

def _haar(X):
    return X[::2,::2]/4.+X[1::2,::2]/4.+X[::2,1::2]/4.+X[1::2,1::2]/4.

def _converge(P,R,b,max_iters_converge=1000):
    Pm=numpy.empty(P.shape,numpy.float32) # Pre-allocate: saves time
    argmax=numpy.zeros_like(P)
    maxval=numpy.zeros(P.shape,numpy.float32)
    for i in xrange(max_iters_converge):
        maxval*=0
        maxval-=1e8
        for m in xrange(int(P.max())+1):
            Pm[:,:]=(P==m)
            Pm=ndimage.gaussian_filter(Pm,b)+R[m] # Pm=ndimage.convolve(Pm,numpy.ones((2*b+1,2*b+1)))+R[m]
            argmax[maxval<Pm]=m
            maxval=maxval*(argmax!=m)+(argmax==m)*Pm
        if (P == argmax).all(): break
        P=argmax.copy()
    mis=defaultdict(xrange(1,int(P.max())+1).__iter__().next)
    mis[0]=0 # set 0 to 0, because it is special
    for i in xrange(P.size):
        P.flat[i]=mis[P.flat[i]]
    return P

def _sigmoid(x):
    return -2./(1.+numpy.exp(-x))+1

def active_masks(f,Kmax,Kmin,A,Mmax,alpha,beta,gamma,b=2,R=None):
    '''
    Regions = active_masks(f,Kmax,Kmin,A,Mmax,alpha,beta,gamma,b=2,R=None)

    Inputs:
    f:      the image
    Kmax:   Maximum k
    Kmin:   Minimum k
    Mmax:   Maximum M
    A:      A[k] should be the schedule of A's
    alpha, beta, gamma: Necessary for smoothing f according to alpha * sigmoid(beta*( (f x g)  - gamma ))
    b:      How much to smooth the regions before voting
    R:      Any object which supports R.randint (default: numpy.random). To control 
            the random initialization

    In this implementation, regions 0 is special and corresponds to the background
    '''
    R = get_random(R)
    f=numpy.asanyarray(f,numpy.float32)
    h,w=f.shape
    need_expand=False
    if h % (2**Kmax):
        nh=(h//2**Kmax+1)*2**Kmax
        need_expand=True
    if w % (2**Kmax):
        nw=(w//2**Kmax+1)*2**Kmax
        need_expand=True
    if need_expand:
        nf=numpy.empty((nh,nw),numpy.float32)
        nf[:h,:w]=f
        f=nf
    else:
        nh,nw=h,w

    nh//=2**Kmax
    nw//=2**Kmax
    P=numpy.array([[R.randint(0,Mmax-1) for i in xrange(nw)] for j in xrange(nh)]) # P is \Psi
    for k in xrange(Kmax,Kmin-1,-1):
        for a in A[k]:
            G1=alpha*_sigmoid(beta*(ndimage.gaussian_filter(f,a)-gamma))
            R0=G1
            for i in xrange(k):
                R0=_haar(R0)
            P=_converge(P,[R0]+[R0*0 for m in xrange(int(P.max()))],b)
        P=_unhaar(P)
    return P[:h,:w]

def active_masks_dna(f,thresh=None,R=None):
    '''
    Regions = active_masks_dna(f,thresh=None,R=None)

    Calls active_masks(f,...) with appropriate parameters for DNA images.
    
    Inputs:
    f:      The image
    thresh: The threshold to use (default: automatically determined)
    R:      Any object with support for R.randint (default: numpy.random)
    '''
    if thresh is None:
        thresh=thresholding.murphy_rc(f)
    return active_masks(f,3,1,{3: [8.,4.,2.], 2 : [8.,4.,2.], 1: [2.]},256,1.2,2*f[f<thresh].std(),thresh,2,R)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
