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
from ..image import Image
import os
from glob import glob
import re
from os.path import exists

def detect_ic100dir(basedir):
    '''
    is_ic100_dir = detect_ic100dir(basedir)

    Returns True if basedir looks like a IC100 basedir
    '''
    return exists('%s/tables' % basedir) and exists('%s/ProtocolArchive' % basedir) and exists('%s/CytoShopLogFile.txt' % basedir)

def read_ic100dir(basedir):
    '''
    imgs = read_ic100dir(basedir)

    Read IC100 output starting on basedir
    '''
    imgs=[]
    wells=glob('%s/tables/well_*' % basedir)
    for well in wells:
        channel_0=glob('%s/channel_0/*.bmp' % well)
        channel_1=glob('%s/channel_1/*.bmp' % well)
        channel_2=glob('%s/channel_2/*.bmp' % well)
        channel_0.sort()
        channel_1.sort()
        channel_2.sort()

        wellpat=re.compile('well__([A-Z])___?([0-9]{1,2})')
        for i,c0,c1,c2 in zip(xrange(len(channel_0)),channel_0,channel_1,channel_2):
            img=Image()
            img.channels[Image.dna_channel]=c0
            img.channels[Image.protein_channel]=c1
            img.channels['autofl']=c2
            name,nr=wellpat.search(well).groups()
            wellname=name + nr
            img.label=wellname
            img.id=(wellname,i)
            imgs.append(img)
    return imgs


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
