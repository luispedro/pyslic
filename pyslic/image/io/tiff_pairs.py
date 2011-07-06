# -*- coding: utf-8 -*-
# Copyright (C) 2010-2011  Murphy Lab, Carnegie Mellon University
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
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
from os import listdir, path
import re
from collections import defaultdict
import pyslic

_name_pat = re.compile('([A-Z0-9]{5})_Image_([0-9]+)_(prot|dna).tiff')
def _channelname(c):
    if c == 'prot':
        return 'protein'
    return c

def read_tiff_pairs_dir(directory):
    '''
    imgs = read_tiff_pairs_dir(directory)

    Read the files in `directory` if it contains TIFF pairs
        (i.e., image_prot.tiff/image_dna.tiff)
    and returns a list of the images found.
    '''
    imgs = defaultdict(pyslic.image.Image)
    for fname in listdir(directory):
        match = _name_pat.match(fname)
        if match is not None:
            code,nr,type = match.groups()
            img = imgs[code,nr]
            img.channels[_channelname(type)] = path.abspath(path.join(directory,fname))
            img.label = code[3:]
            img.id = (img.label, nr)
    return imgs.values()


def detect_tiff_pairs(directory):
    '''
    '''
    for img in listdir(directory):
        if img.endswith('_prot.tiff'):
            base = img[:-len('_prot.tiff')]
            return path.exists(path.join(directory, base+'_dna.tiff'))
