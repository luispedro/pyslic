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
import numpy
import pyslic
import pickle
from pyslic.image.io import readimg
import os
from os.path import dirname
basedir=dirname(__file__)
protimg=(basedir+'/data/protimg.jp2')

def test_readimg():
    img = readimg(protimg)
    assert img.shape == (1024, 1344)

def _buildimg():
    img = pyslic.Image()
    img.channels[img.protein_channel] = protimg
    return img

def test_image_load():
    img = _buildimg()
    assert img.protein_channel not in img.channeldata
    img.load()
    assert img.loaded

    assert img.dna_channel not in img.channeldata
    img_direct = readimg(protimg)
    assert numpy.all(img_direct == img.channeldata[img.protein_channel])

    img.unload()
    assert not img.loaded
    assert img.protein_channel not in img.channeldata
    assert img.dna_channel not in img.channeldata

def test_image_pickle():
    img = _buildimg()
    pickle_filename='/tmp/pyslic.test.img.pp'
    pickle.dump(img,file(pickle_filename,'w'))
    img2=pickle.load(file(pickle_filename))
    assert img.channels[img.protein_channel] == img2.channels[img.protein_channel]
    assert not img2.loaded
    img.load()
    pickle.dump(img,file(pickle_filename,'w'))
    img2=pickle.load(file(pickle_filename))
    assert not img2.loaded
    assert img2.protein_channel not in img2.channeldata
    assert img2.dna_channel not in img2.channeldata
    os.unlink(pickle_filename)
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
