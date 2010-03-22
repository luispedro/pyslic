# -*- coding: utf-8 -*-
# Copyright (C) 2009  Murphy Lab
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
from dirtransversal import *
from ic100 import *
from ksr import *
from readtjz import *
from tiff_pairs import *

def auto_detect_load(basedir):
    '''
    imgs = auto_detect_load(basedir)

    Auto-detect the type of directory and load
    images from it.

    Returns None if it cannot read the directory

    Currently supported:
    --------------------
        * dirtransversal
        * ic100
        * ksr
        * TJZ containing directory
        * TIFF pairs directory
    '''
    if detect_dirtransversal(basedir):
        return dirtransversal(basedir)
    if detect_ic100dir_flat(basedir):
        return read_ic100dir_flat(basedir)
    if detect_ic100dir(basedir):
        return read_ic100dir(basedir)
    if detect_ksr_dir(basedir):
        return read_ksr_dir(basedir)
    if detect_tjzdir(basedir):
        return readtjz_recursive(basedir)
    if detect_tiff_pairs(basedir):
        return read_tiff_pairs_dir(basedir)
    return None
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
