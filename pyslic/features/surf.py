# -*- coding: utf-8 -*-
# Copyright (C) 2011  Murphy Lab
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
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

from __future__ import division, with_statement
import numpy as np
import mahotas.surf

def surf_ref(f, ref):
    fi = mahotas.surf.integral(f.copy())
    points = mahotas.surf.interest_points(fi, 6, 24, 1, max_points=1024, is_integral=True)
    descs = mahotas.surf.descriptors(fi, points, is_integral=True, descriptor_only=True)
    if ref is None:
        return descs
    ri = mahotas.surf.integral(ref.copy())
    descsref = mahotas.surf.descriptors(ri, points, is_integral=True, descriptor_only=True)
    return np.hstack( (descs, descsref) )

