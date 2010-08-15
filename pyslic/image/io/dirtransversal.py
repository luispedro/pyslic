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
from collections import defaultdict
import os.path
from os.path import isdir
from ..image import Image

__all__ = ['loadimages','detect_dirtransversal','dirtransversal']

def _validfilename(f):
    return f.endswith('.png') or f.endswith('.tif') or f.endswith('.bmp')

def detect_dirtransversal(startdir):
    '''
    is_dirtransversal = detect_dirtransversal(startdir)

    Returns true if startdir seems like the start of a dirtransversal
    directory.
    '''
    count_pos = 0
    for d in os.listdir(startdir):
        d = os.path.join(startdir,d)
        if isdir(d):
            subs = os.listdir(d)
            if not 'prot' in subs:
                return False
            for sub in subs:
                if sub not in ('prot','dna','crop'): return False
            count_pos += 1
    return bool(count_pos)

def dirtransversal(startdir):
    """
    images = loadimages(startdir)

    Load images which are stored in files below directory startdir
    This function expects a directory structure like:
    
    startdir/
        class1/
            dna/
                image1
                image2
                ...
            prot/
                image1
                image2
                ...
            crop/
                image1
                image2
                ...
        class2/
            ....
    """
    assert type(startdir) is not unicode, 'pyslic.image.io.dirtransversal does not work with unicode input' # The problem is that it creates images with unicode paths which cannot be loaded!
    images = []
    labelcount = defaultdict(int)
    def _loadimage(d,dna,prot,crop):
#        print 'loadimage(%s,%s,%s,%s)' % (d,dna,prot,crop)
        res = Image()
        res.label = d
        res.id = (d, labelcount[d])
        labelcount[d] += 1
        res.channels['dna'] = os.path.abspath(os.path.join(startdir,d,'dna',dna))
        res.channels['protein'] = os.path.abspath(os.path.join(startdir,d,'prot',prot))
        if crop:
            res.channels['crop']= os.path.abspath(os.path.join(startdir,d,'crop',crop))
        return res
    for d in sorted(os.listdir(startdir)):
        if not isdir(os.path.join(startdir,d)): continue
        D = os.listdir(os.path.join(startdir,d,'dna'))
        D.sort()
        P = os.listdir(startdir+'/'+d+'/prot')
        P.sort()
        C = None
        if isdir(os.path.join(startdir,d,'crop')):
            C = os.listdir(os.path.join(startdir,d,'crop'))
            C.sort()

        Di,Pi,Ci = 0,0,0
        while Di < len(D):
            while not _validfilename(D[Di]): Di += 1
            while not _validfilename(P[Pi]): Pi += 1
            if C:
                while not _validfilename(C[Ci]): Ci += 1
                images.append(_loadimage(d,D[Di],P[Pi],C[Ci]))
            else:
                images.append(_loadimage(d,D[Di],P[Pi],None))
            Di += 1
            Pi += 1
            Ci += 1
    return images

loadimages = dirtransversal

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
