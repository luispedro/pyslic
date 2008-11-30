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

__all__ = ['classification_result','normaliselabels','classifier']


class classification_result(object):
    __slots__ = ['label','p','funcvalue']
    def __init__(self):
        pass

def normaliselabels(labels):
    labelnames={}
    normalised=[]
    names=[]
    N=0
    for L in labels:
        nr=labelnames.get(L,N)
        if nr == N:
            labelnames[L]=N
            names.append(L)
            N += 1
        normalised.append(nr) 
    return array(normalised),names

class classifier(object):
    def __init__(self):
        self._trained = False

    def train(self,features,labels):
        nlabels,self._labelnames=normaliselabels(labels)
        self._nclasses = nlabels.max()+1
        features=asanyarray(features)
        self._dotrain(features,nlabels)
        self._trained = True
    
    def apply(self,features):
        assert self._trained
        features=asanyarray(features)
        if features.ndim == 1:
            return self._labelnames[int(self._doapply(features))]
        else:
            return array([self._labelnames[int(self._doapply(features[i]))] for i in xrange(features.shape[0])])

class results_to_labels(classifier):
    def __init__(self):
        pass
    def _dotrain(self,features,labels):
        pass

    def _doapply(self,res):
       return res.label 


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:

