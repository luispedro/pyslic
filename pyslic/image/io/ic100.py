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
from ..image import Image
import os
from glob import glob
import re
from collections import defaultdict
from os.path import exists, abspath
import random

def detect_ic100dir(basedir):
    '''
    is_ic100_dir = detect_ic100dir(basedir)

    Returns True if basedir looks like a IC100 basedir
    '''
    return exists('%s/tables' % basedir) and exists('%s/ProtocolArchive' % basedir)

def read_ic100_BMP(input):
    '''
    img = read_ic100_BMP(file_or_filename)

    Reads a BMP from the IC100. If it's 8-bits, then this is just the traditional
    Windows BMP format. However, 16-bit grey-scale BMPs are not correctly handled by
    traditional software and are handled by our code.

    Reference: http://en.wikipedia.org/wiki/BMP_file_format
    '''
    def do_close(): pass
    if type(input) == str:
        input=file(input)
        do_close = input.close
    try:
        def get2():
            S=input.read(2)
            return ord(S[0])+(ord(S[1])<<8)
        def get4():
            S=input.read(4)
            return ord(S[0])+(ord(S[1])<<8)+(ord(S[2])<<16)+(ord(S[3])<<24)

        # See http://en.wikipedia.org/wiki/BMP_file_format
        magick = input.read(2)
        if magick != 'BM': raise IOError, 'read_ic100_BMP: Unknown file format'
        size=get4()
        get2()
        get2()
        offset=get4()
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
        if hsize != 40:
            raise IOError, 'Header size is not 40'
        if bitsppixel == 8:
            img = numpy.empty((height,width),numpy.uint8)
        elif bitsppixel == 16:
            img = numpy.empty((height,width),numpy.dtype('<H'))
        else:
            raise IOError, 'Bits per pixel is not 8 or 16'
        input.seek(offset)
        img.data[:] = input.read() # This is very fast
        # The scan order needs to be fixed
        return img[::-1,:]
    except TypeError:
        raise IOError("Cannot load image '%s' (TypeError)" % input)
    finally:
        do_close()

def read_ic100dir(basedir):
    '''
    imgs = read_ic100dir(basedir)

    Read IC100 output starting on basedir
    '''
    assert type(basedir) is not unicode, 'pyslic.image.io.read_ic100dir does not work with unicode input' # The problem is that it creates images with unicode paths which cannot be loaded!
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
            img.load_function=read_ic100_BMP
            img.channels['dna']=c0
            img.channels['protein']=c1
            img.channels['autofl']=c2
            name,nr=wellpat.search(well).groups()
            wellname=name + nr
            img.label=wellname
            img.id=(wellname,i)
            imgs.append(img)
    return imgs

_flat_pat = re.compile(r'(^|/)[0-9]{6}[A-Z_]+[0-9]__([A-H])___?([0-9]{1,2})_T_001_ch_0([012])_image_0+([1-9][0-9]?)_Z_001\.bmp$')
_channels = ('dna','protein','autofluorescence')

def detect_ic100dir_flat(basedir):
    '''
    is_flat_ic100dir = detect_ic100dir_flat(basedir)

    Returns whether it seems like read_ic100dir_flat can
    parse the structure of basedir.
    '''
    bmps = glob(basedir+'/*.bmp')
    if not bmps: return False
    random.shuffle(bmps)
    for b in bmps[:2+len(bmps)//10]:
        if _flat_pat.search(b) is not None:
            return True
    return False

def read_ic100dir_flat(basedir):
    '''
    imgs = read_ic100dir_flat(basedir)

    Read images from IC100 flat directory.
    '''
    wellfiles = defaultdict(dict)
    for name in glob(basedir+'/*.bmp'):
        match = _flat_pat.search(name)
        if match is not None:
            _,wellrow,wellcol,channelnr,imageid = match.groups()
            channel = _channels[int(channelnr)]
            wellfiles[wellrow+wellcol,imageid][channel] = name
    images = []
    for (well,idx),channels in wellfiles.iteritems():
        img = Image()
        for ch,f in channels.iteritems():
            img.channels[ch] = abspath(f)
        img.id = (well,idx)
        img.label = well
        img.load_function = read_ic100_BMP
        images.append(img)
    images.sort(key=lambda Img: Img.label)
    return images
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
