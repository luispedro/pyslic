# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
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

from __future__ import division
from mahotas.tas import tas, pftas

def pftasinfo():
    return [('SLF31.%s' % (i+1),'pftas:center_%s' % i,2,1) for i in xrange(9)] + \
            [('SLF31.%s' % (i+10),'npftas:center_%s' % i,2,1) for i in xrange(9)] + \
            [('SLF33.%s' % (i+1),'pftas:mu_margin_%s' % i,2,1) for i in xrange(9)] + \
            [('SLF33.%s' % (i+1+2*9),'npftas:mu_margin_%s' % i,2,1) for i in xrange(9)] + \
            [('SLF33.%s' % (i+1+9),'pftas:mu_%s' % i,2,1) for i in xrange(9)] + \
            [('SLF33.%s' % (i+1+3*9),'npftas:mu_%s' % i,2,1) for i in xrange(9)]

