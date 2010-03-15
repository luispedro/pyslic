# -*- coding: utf-8 -*-
# Copyright (C) 2010  Murphy Lab, Carnegie Mellon University
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
import pyslic

def read_tiff_pairs_dir(directory):
    '''
    imgs = read_tiff_pairs_dir(directory)

    Read the files in `directory` if it contains TIFF pairs
        (i.e., image_prot.tiff/image_dna.tiff)
    and returns a list of the images found.
    '''
    pairs = {}
    for img in listdir(directory):
        if img.endswith('_prot.tiff'):
            pairs[img[:-len('_prot.tiff')], 'protein'] = img
        elif img.endswith('_dna.tiff'):
            pairs[img[:-len('_dna.tiff')], 'dna'] = img

    imgs = []
    for label,channel in pairs:
        if channel != 'protein': continue
        if (label, 'dna') not in pairs:
            raise IOError("Unmatched protein pair: %s" % label)
        img = pyslic.Image()
        img.channels['dna'] = path.join(directory, pairs[label, 'dna'])
        img.channels['protein'] = path.join(directory, pairs[label, 'protein'])
        img.label = label
        imgs.append(img)
    return imgs


def detect_tiff_pairs(directory):
    '''
    '''
    for img in listdir(directory):
        if img.endswith('_prot.tiff'):
            base = img[:-len('_prot.tiff')]
            return path.exists(path.join(directory, base+'_dna.tiff'))
