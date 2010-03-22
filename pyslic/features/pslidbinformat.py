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

from __future__ import division, with_statement
from struct import pack, unpack

__all__ = ['readpslidbin','writepslidbin']
def readpslidbin(filename):
    input = file(filename)
    def read_int32():
        bytes = input.read(4)
        return unpack('>i',bytes)[0]
    def read_float32():
        bytes = input.read(4)
        return unpack('>f',bytes)[0]
    def read_str():
        nbytes = read_int32()
        return input.read(nbytes)
    def read_strlist():
        return read_str().split('@')[:-1]
    def read_intlist():
        n = read_int32()
        return [read_int32() for i in xrange(n)]
    version = read_int32()
    if version != 5:
        raise NotImplementedError('pyslic.readpslidbin: Can only read version 5 files')
    rows = read_int32()
    cols = read_int32()
    features = [read_float32() for i in xrange(cols)]
    featids = read_intlist()
    real_slf_names = read_strlist()
    slf_names = read_strlist()
    names = read_strlist()
    n_sampleids = read_int32()
    samples = [read_int32() for i in xrange(n_sampleids)]
    settype = read_int32()
    imageurls = read_strlist()
    maskurls = read_strlist()
    #ndims = read_int32()  # This is in the documentation but not in the matlab code
    n_channels = read_int32()
    channel_nrs = [read_int32() for i in xrange(n_channels)]
    return features, real_slf_names, slf_names, names, imageurls, maskurls
        
def writepslidbin(output, features, real_slf_names, slf_names, names, imageurls, maskurls, settype, channel_nrs):
    if type(output) in (str,unicode):
        output = file(output, 'w')
    def write_int32(x):
        output.write(pack('>i',x))
    def write_float32(x):
        output.write(pack('>f',x))
    def write_strlist(s):
        s = ''.join(a+'@' for a in s)
        write_int32(len(s))
        output.write(s)
    def write_intlist(l):
        write_int32(len(l))
        for i in l:
            write_int32(i)
    if len(features.shape) == 1:
        features = features.reshape((1,features.size))
    rows,cols = features.shape
    write_int32(5) # Version
    write_int32(rows)
    write_int32(cols)
    for fs in features:
        for f in fs:
            write_float32(f)
    write_intlist([])
    write_strlist(real_slf_names)
    write_strlist(slf_names)
    write_strlist(names)
    write_intlist([])
    write_int32(settype)
    write_strlist(imageurls)
    write_strlist(maskurls)
    write_intlist(channel_nrs)
    output.close()

