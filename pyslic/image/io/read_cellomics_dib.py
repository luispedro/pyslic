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
import numpy

__all__ = ['read_cellomics_dib']

def read_cellomics_dib(input):
    '''
    Reads an image in Cellomics DIB format.

    Returns a numpy uint16 array.

    Cellomics DIB format is slightly different than Window DIB format,
    as the data is 16-bit greyscale (as opposed to colour). Therefore,
    tools such as ImageMagick parse it incorrectly.
    '''

    if type(input) == str:
        input=file(input)

    def get2():
        S=input.read(2)
        return ord(S[0])+(ord(S[1])<<8)
    def get4():
        S=input.read(4)
        return ord(S[0])+(ord(S[1])<<8)+(ord(S[2])<<16)+(ord(S[3])<<24)

    # See http://en.wikipedia.org/wiki/BMP_file_format#Bitmap_information_.28DIB_header.29
    hsize=get4()
    width=get4()
    height=get4()
    planes=get2()
    bitsppixel=get2()
    compression=get4()
    imgsize=get4()
    hres=get4()
    vres=get4()
    colours=get4()
    importantcolours=get4()

    # I don't actually know what the 12 bytes are, but it seems it is some
    # additional header which is safe to ignore
    input.read(12)
    if hsize != 40:
        raise IOError, 'Header size is not 40'
    if bitsppixel != 16:
        raise IOError, 'Bits per pixel is not 16'
    img=numpy.zeros((width,height),numpy.dtype('<H'))
    img.data[:]=input.read(2*width*height) # This is very fast

    # The scan order needs to be fixed
    return img[::-1,:]

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
