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
import re
import os.path
import os
import sys
from ..image import Image
from .read_cellomics_dib import read_cellomics_dib
from warnings import warn

__all__ = ['read_ksr_dir','read_ksrdir','detect_ksrdir', 'detect_ksr_dir']
_ksrpat=re.compile('KSR_.*t([0-9]+)([A-H][0-9]{1,2})f([0-9]+)d([0-9])\.(tif|TIF|DIB|dib)')

def detect_ksrdir(dir):
    '''
    is_ksrdir = detect_ksrdir(dir)

    Returns true if the directory seems to contain ksr files
    '''
    f=os.listdir(dir)[0]
    return _ksrpat.match(f)

def _fixlabel(L):
    '''
    L = _fixlabel(L)

    Turns B02 into B2, while preserving B11
    '''
    if L[1] == '0':
        return L[0]+L[2]
    return L

def read_ksrdir(dir):
    '''
    images = read_ksrdir(dirname)

    Read all the files in dirname and return them as a dictionary:
        (WellName, FieldNr) -> Image
    '''
    assert type(dir) is not unicode, 'pyslic.image.io.read_ksrdir does not work with unicode input' # The problem is that it creates images with unicode paths which cannot be loaded!
    Files=os.listdir(dir)
    channelcode = { 1 : Image.dna_channel , 2 : Image.protein_channel, 3 : Image.autofluorescence_channel } 
    images={}
    for f in Files:
        m=_ksrpat.match(f)
        if not m:
            print "Don't know how to process", f
            continue
        T,Well,Field,Channel,_=m.groups()
        T=int(T)
        Field=int(Field)
        if T != 1:
            warn("Don't know how to handle more than one time point, ignoring timepoints with t > 1")
            continue
        img=images.get((Well,Field),None)
        isdib = f.endswith('.DIB') or f.endswith('.dib')
        if img is None:
            img = Image()
            img.label = _fixlabel(Well)
            img.id = (img.label,Field)
            images[(Well,Field)]=img

        # If there's a DIB and a TIFF, take the DIB!
        # The logic below is simplified, but it works for our dataset:
        if isdib:
            img.load_function = read_cellomics_dib
        if (channelcode[int(Channel)] not in img.channels) or isdib:
            img.channels[channelcode[int(Channel)]]=os.path.abspath(os.path.join(dir,f))
    return list(images.values())


read_ksr_dir = read_ksrdir
detect_ksr_dir = detect_ksrdir

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
