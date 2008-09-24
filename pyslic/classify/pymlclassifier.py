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

from numpy import array
from PyML import *

__ALL__=['PyMLSVM']
class PyMLSVM(object):
    def __init__(self):
        pass

    def train(self,feats,labels):
        allabels=dict([(L,None) for L in labels])
        if len(allabels) > 2:
            self.s = multi.OneAgainstRest(svm.SVM())
        else:
            self.s = svm.SVM()
        nlabels=[]
        for L in labels:
            nlabels.append(str(L))
        data=datafunc.VectorDataSet(feats,L=nlabels)
        self.s.train(data)
        self.trained=True

    def apply(self,feats):
        assert self.trained
        if feats.ndim == 1:
            data=datafunc.VectorDataSet([feats])
            return array(self.s.classify(data,0))[0]
        else:
            return array([self.apply(feats[i]) for i in xrange(feats.shape[0])])

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
