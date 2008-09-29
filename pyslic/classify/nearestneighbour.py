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
from .classifier import classifier, classification_result

__all__ = ['NNClassifier']

class NNClassifier(classifier):
    '''
    Nearest Neighbour Classifier
    ----------------------------

    Implements a nearest neighbour classifier. Currently, uses the simplest
    implementation with O(N) classification time.
    '''
    def __init__(self):
        classifier.__init__(self)
    
    def _dotrain(self,features,labels):
        self.training=features.copy()
        self.labels=labels.copy()

    def _doapply(self,feats):
        dists = ((self.training-feats)**2).sum(1)
        idx = dists.argmin()
        return self.labels[idx]

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
