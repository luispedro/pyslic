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
import numpy as np
from numpy import linalg
from collections import defaultdict
try:
    import ncreduce
    fast_all = ncreduce.all
except:
    fast_all = np.all
from ..imageprocessing import thresholding
from mahotas.polygon import fill_convexhull as convexhull
from mahotas.bbox import bbox, croptobbox
from .. import features
from scipy import ndimage
import pymorph
import mahotas

__all__ = ['roysam_watershed']
def roysam_watershed(dna,thresh=None,blur_factor=3):
    '''
    Run watershed on mixed gradient & intensity image as suggested by Lin et al.

    -Input
    dna:            DNA image
    thresh:         Gray value threshold (default: computed using Murphy's RC)
    blur_factor:    Blur factor (default: 3)
    
    
    REFERENCE
    Gang Lin, Umesh Adiga, Kathy Olson, John F. Guzowski, Carol A. Barnes, and Badrinath Roysam
    "A Hybrid 3-D Watershed Algorithm Incorporating Gradient Cues & Object Models for Automatic
        Segmentation of Nuclei in Confocal Image Stacks"
     Vol. 56A, No. 1, pp. 23-36 Cytometry Part A, November 2003.
    '''
    if thresh is None:
        thresh = 'murphy_rc'
    M = (ndimage.gaussian_filter(dna,4) > thresholding.threshold(dna,thresh))
    G = pymorph.gradm(dna)
    D = ndimage.distance_transform_edt(M)
    D *= np.exp(1-G/float(G.max()))
    T = ndimage.gaussian_filter(D.max() - D,blur_factor)
    if T.max() < 256:
        T = pymorph.to_uint8(T)
    else:
        T = pymorph.to_uint8(T*(256.0/T.max()))
    T *= M
    R = pymorph.regmin(T)
    R *= M
    R,N = ndimage.label(R)
    R[(R==0)&(M==0)] = N+1
    W,WL = mahotas.cwatershed(T,R,return_lines=True)
    W *= M
    return W,WL

def border(W,Bc=None):
    '''
    B, neighbours, border_id = border(Regs,Bc=None)

    Label the borders between regions.

    B[y,x] = border_id[r0,r1]
    
    if (y,x) is on the border between regions (r0,r1).
    neighbours is a dict : (region) -> (list of regions)

    -Input:
    Regs:   A labeled image
    Bc:     A structuring element (default: 3x3 cross)
    '''
    bg = W.max()
    B = np.zeros_like(W)
    neighbours = defaultdict(set)
    border_id = defaultdict(xrange(1,W.max()**2).__iter__().next)
    for obji in xrange(1,W.max()):
        B_obji = (mahotas.dilate(W == obji,Bc)-(W==obji))
        for neighbour in np.unique(W[B_obji]):
            if neighbour == 0 or neighbour == bg: continue
            a1,a2 = obji,neighbour
            if a2 < a1: a1,a2 = a2,a1
            B[B_obji & (W==neighbour)] = border_id[a1,a2]
            neighbours[a1].add(a2)
            neighbours[a2].add(a1)
    return B, dict(neighbours), dict(border_id)

def _compute_features(img):
    bimg = (img > 0)
    s00,s01,s10,s11 = bbox(bimg)
    if s00 > 0: s00 -= 1
    if s10 > 0: s10 -= 1
    bimg = bimg[s00:s01+1,s10:s11+1]
    hull = convexhull(bimg - pymorph.erode(bimg))
    Allfeats = np.r_[features.hullfeatures.hullfeatures(img,hull),features.hullfeatures.hullsizefeatures(img,hull)]
    return Allfeats[np.array([0,1,2,3,5,6,7],int)]

class Merger(object):
    '''
    Implements Roysam's region merging algorithm.

    REFERENCE
    Gang Lin, Umesh Adiga, Kathy Olson, John F. Guzowski, Carol A. Barnes, and Badrinath Roysam
    "A Hybrid 3-D Watershed Algorithm Incorporating Gradient Cues & Object Models for Automatic
        Segmentation of Nuclei in Confocal Image Stacks"
     Vol. 56A, No. 1, pp. 23-36 Cytometry Part A, November 2003.
    '''

    def __init__(self,dna,mu,iSigma,thresh=None):
        '''
        M = Merged(dna)
        '''
        self.mu = mu
        self.iSigma = iSigma
        self.C = pymorph.gradm(dna)
        self.W,self.WL = roysam_watershed(dna,thresh)
        self.B,self.neighbours,self.border_id = border(self.W)
        self.border_regions=dict((y,x) for x,y in self.border_id.iteritems())

    def _merge(self,c0,c1):
        '''Merge region c0 & c1'''
        self.W[self.W==c1] = c0
        self.neighbours[c0].remove(c1)
        self.neighbours[c1].remove(c0)
        self.neighbours[c0].update(self.neighbours[c1])
        del self.neighbours[c1]
        for v in self.neighbours.values():
            if c1 in v:
                v.remove(c1)
                v.add(c0)
        b = self.border_id[min(c0,c1),max(c0,c1)]
        del self.border_id[min(c0,c1),max(c0,c1)]
        self.B[self.B==b] = 0
        for (r0,r1),b in self.border_id.items():
            if r0 == c1 or r1 == c1:
                del self.border_id[r0,r1]
                if r1 == c1:
                    r1 = r0
                d0,d1 = min(r1,c0),max(r1,c0)
                if (d0,d1) in self.border_id:
                    self.B[self.B==b] = self.border_id[d0,d1]
                else:
                    self.border_id[d0,d1] = b
                
    def _Rw(self,c0,c1):
        '''Implement Rw in the paper.'''
        def RSw(c0,c1):
            def logS(img):
                F = _compute_features(img)
                return -.5*np.sqrt( np.dot(np.dot(F-self.mu,self.iSigma),F-self.mu) )
            S0 = logS(self.W == c0)
            S1 = logS(self.W == c1)
            Sc = logS((self.W == c0)|(self.W == c1))
            logR = np.log(2)+Sc-S0-S1
            if abs(logR) < 100:
                return np.exp(logR)
            return np.exp(np.sign(logR)*100)
        def RGw(c0,c1):
            b = self.border_id[min(c0,c1),max(c0,c1)]
            if np.all(self.B != b):
                return .01 # A small number: these probably should not be merged
            return (self.C[self.W==c0].mean()+self.C[self.W==c1].mean())/2./self.C[self.B==b].mean()
        rs = RSw(c0,c1)
        rg = RGw(c0,c1)
        return rs * rg

    def greedy(self,beta=1.2):
        '''
        Merge regions greedily.
        '''
        values = dict(((c0,c1),self._Rw(c0,c1)) for (c0,c1),b in self.border_id.iteritems())
        while True:
            if not values: break
            val,best = max((y,x) for x,y in values.iteritems())
            if val < beta: break
            b0,b1 = best
            affected = self.neighbours[b0]|self.neighbours[b1]
            del values[b0,b1]
            for c0,c1 in values.keys():
                if c0 in affected or c1 in affected:
                    del values[c0,c1]
            self._merge(b0,b1)
            for c0,c1 in self.border_id:
                if (c0,c1) not in values:
                    values[c0,c1] = self._Rw(c0,c1)
        self.W[self.W == self.W.max()] = 0
        return self.W


def train_classifier(labeled_dnas):
    '''
    classifier = train_classifier(labeled_dnas)
    
    Learn a classifier for Roysam's Algorithm.
        (to be used for greedy_roysam_merge)
    '''
    F = []
    for dna in labeled_dnas:
        for obji in xrange(1,dna.max()+1):
            F.append(_compute_features(dna==obji))
    F = np.array(F)
    return np.mean(F,0),linalg.inv(np.cov(F.T))

def greedy_roysam_merge(dna, classifier, thresh=None):
    '''
    labeled = greedy_roysam_merge(dna, classifier, thresh=None)

    Implement Roysam's Greedy Merging Algorithm

    Parameters
    -----------
    
        * dna: DNA image
        * classifier: shape properties
            (should be of the form returned by train_classifier)
        * thresh: Thresholding value or method
    '''
    mu,iSigma = classifier
    M = Merger(dna,mu,iSigma,thresh=thresh)
    M.greedy()
    return M.W

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
