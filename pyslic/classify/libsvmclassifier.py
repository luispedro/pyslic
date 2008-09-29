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
    def __init__(self):
        classifier.__init__(self)
        if _svm is None:
            raise RunTimeError('SVM Library not found. Cannot use this classifier.')
        self.param = _svm.svm_parameter(kernel_type = svm.RBF)
    
    def set_option(self,optname,value):
        setattr(self.param,optname,value)

    def _dotrain(self,features,labels):
        problem=_svm.svm_problem(labels,features)
        self.model=_svm.svm_model(problem,self.param)

    def _doapply(self,feats):
        return self.model.predict(feats)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
