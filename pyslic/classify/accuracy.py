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

__all__ = ['accuracy','waccuracy','rowsto1']

def rowsto1(cmatrix):
    return (cmatrix.T / cmatrix.sum(1)).T

def accuracy(cmatrix):
    '''
    acc = accuracy(cmatrix)

    Accuracy of cmatrix
    '''
    return cmatrix.trace()/cmatrix.sum()

def waccuracy(cmatrix):
    '''
    wacc = waccuracy(cmatrix)

    Weighted accuracy of cmatrix
    '''
    return (cmatrix.diagonal() / cmatrix.sum(1)).mean()
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
