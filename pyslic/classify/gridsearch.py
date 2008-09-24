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
from nfoldcrossvalidation import *

def _allassignments(D):
    def allassignmentslist(L):
        if len(L) == 0:
            yield []
            return
        key0,vals0=L[0]
        for v in vals0:
            for rest in allassignmentslist(L[1:]):
                yield [(key0,v)]+rest
    for A in allassignmentslist(list(D.items())):
        yield A

def _set_assignment(obj,assignments):
    if not hasattr(obj,'set_option'):
        raise "Don't know how to set options"
    for k,v in assignments:
        obj.set_option(k,v)

class gridsearch(object):
    def __init__(self,base,**kwargs):
        self.params = kwargs
        self.base = base
        self.best = None

    def is_multi_class(self):
        return self.base.is_multi_class()

    def train(self,features,labels):
        best_trace=-1
        for assignement in _allassignments(self.params):
            _set_assignment(self.base,assignement)
            S,_=nfoldcrossvalidation(features,labels,classifier=self.base)
            if S.trace() > best_trace:
                self.best=assignement
                best_trace=S.trace()
        _set_assignment(self.base,self.best)

    def apply(self,features):
        return self.base.apply(features)
