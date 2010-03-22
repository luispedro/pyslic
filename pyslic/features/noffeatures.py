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
__all__ = ['noffeatures']
def noffeatures(procimg,nofimg):
    """
    Compute non object fluorescence features
    """
    obj_fluor = procimg.sum()
    nonobj_fluor = nofimg.sum()
    total_fluor = (obj_fluor + nonobj_fluor)
    if total_fluor == 0:
        return 1.
    return nonobj_fluor/total_fluor

noffeatures.names = ['fract_nonobj_fluor']
noffeatures.slf_names = ['SLF7.79']
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
