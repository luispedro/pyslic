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
from numpy import array, asanyarray
from classifier import classifier, classification_result
from tempfile import NamedTemporaryFile
_svm = None
try:
    import svm
    _svm = svm
except:
    pass
try:
    from libsvm import svm
    _svm = svm
except:
    pass

__all__=['libsvmClassifier']
class libsvmClassifier(classifier):
    def __init__(self,probability = False, auto_weighting = True):
        classifier.__init__(self)
        if _svm is None:
            raise RuntimeError('SVM Library not found. Cannot use this classifier.')
        self.param = _svm.svm_parameter(kernel_type = svm.RBF, probability = probability)
        self.output_probability = probability
        self.auto_weighting = auto_weighting
    
    def set_option(self,optname,value):
        setattr(self.param,optname,value)

    def _dotrain(self,features,labels):
        if self.auto_weighting:
            nlabels = labels.max() + 1
            self.param.nr_weight = int(nlabels)
            self.param.weight_label = range(nlabels)
            self.param.weight = [(labels != i).mean() for i in xrange(nlabels)]
        problem = _svm.svm_problem(labels.astype(float), features)
        self.model = _svm.svm_model(problem,self.param)

    def apply(self,feats):
        if len(feats.shape) == 2: return [self.apply(f) for f in feats]
        if self.output_probability:
            return self.model.predict_probability(feats)
        return self.model.predict(feats)
    
    def __getstate__(self):
        # This is really really really hacky, but it works
        N=NamedTemporaryFile()
        self.model.save(N.name)
        S=N.read()
        return S,self.output_probability,self._trained,self._labelnames

    def __setstate__(self,state):
        if _svm is None:
            raise RuntimeError('LibSVM Library not found. Cannot use this classifier.')
        S,self.output_probability,self._trained,self._labelnames = state
        N=NamedTemporaryFile()
        N.write(S)
        N.flush()
        self.model = _svm.svm_model(N.name)


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
