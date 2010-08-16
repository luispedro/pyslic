# -*- coding: utf-8 -*-
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
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

from mahotas.texture import haralick as haralickfeatures

__all__ = ['haralickfeatures']

haralickfeatures.names=[
    'angular_second_moment',
    'contrast','correlation',
    'sum_of_squares',
    'inverse_diff_moment',
    'sum_avg',
    'sum_var',
    'sum_entropy',
    'entropy',
    'diff_var',
    'diff_entropy',
    'info_measure_corr_1',
    'info_measure_corr_2']
haralickfeatures.slfnames = [ ('SLF3.%d' % n) for n in xrange(66,66+13)]

