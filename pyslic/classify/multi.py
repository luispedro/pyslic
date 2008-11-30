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
from classifier import classifier
import random

class one_against_rest(classifier):
    __slots__ = ['_classifiers','_base','_base2']

    def __init__(self,base):
        classifier.__init__(self)
        self._classifiers = None
        self._base = base
        self._base2 = None

    def is_multi_class(self):
        return True

    def _dotrain(self,features,labels):
        if labels.max() + 1 == 2:
            self._base2 = self._base()
            self._base2.train(features,labels)
            return
        #else
        self._classifiers=[]
        for i in xrange(labels.max()+1):
            s = self._base()
            s.train(features,labels==i)
            self._classifiers.append(s)

    def _doapply(self,feats):
        if self._classifiers is None:
            return self._base2.apply(feats)
        vals=array([c.apply(feats) for c in self._classifiers])
        idxs, = where(vals == 1)
        if len(idxs) == 1:
            return idxs[0]
        elif len(idxs) == 0:
            return random.randint(0,self._nclasses)
        else:
            # choose at random
            return random.choice(idxs)

class one_against_one(classifier):
    __slots__ = ['_classifiers','_base','_labelnames']

    def __init__(self,factory):
        classifier.__init__(self)
        self._classifiers = None
        self._base = factory

    def _is_multiclass(self):
        return True

    def _dotrain(self,features,labels):
        nc = labels.max()+1
        self._nclasses = nc
        if nc == 2:
            self._base = self._base()
            self._base.train(features,labels)
            return
        #else
        self._classifiers=[[None for i in xrange(nc)] for j in xrange(nc)]
        for i in xrange(nc):
            for j in xrange(i+1,nc):
                s=self._base()
                idxs=(labels == i) | (labels == j)
                s.train(features[idxs],labels[idxs]==i)
                self._classifiers[i][j]=s

    def _doapply(self,feats):
        if self._classifiers is None:
            return self._base.apply(feats)
        nc=self._nclasses
        votes=zeros(nc)
        for i in xrange(nc):
            for j in xrange(i+1,nc):
                c=self._classifiers[i][j].apply(feats)
                if c:
                    votes[i] += 1
                else:
                    votes[j] += 1
        return votes.argmax(0)
    
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
