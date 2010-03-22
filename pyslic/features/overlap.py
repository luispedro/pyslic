# -*- coding: utf-8 -*-
# Copyright (C) 2009  Murphy Lab, Carnegie Mellon University
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
import numpy as np
from scipy import ndimage

def _corrcoef(a,b):
    return np.corrcoef(a.ravel(), b.ravel())[1,0]


def overlapfeatures(protein, reference, procprotein, procreference):
    """
    values = overlapfeatures(protein, reference, procprotein, procreference)
    """
    assert procreference is not None, 'pyslic.features.overlapfeatures: reference must not be None'
    binprot = (procprotein > 0)
    binref = (procreference > 0)
    dist = ndimage.distance_transform_edt(~binref)
    return np.array([
        binprot.sum()/float(binref.sum()),
        binprot[binref].mean(),
        protein[binref].sum()/protein.sum(),
        procprotein[binref].sum()/procprotein.sum(),
        binref[binprot].mean(),
        _corrcoef(binprot,binref),
        _corrcoef(protein,binref),
        _corrcoef(protein,reference),
        np.median((dist * procprotein)[binprot]),
        np.mean((dist * procprotein)[binprot]),
        ])

def overlapinfo():
    return [
        ('SLF34.1',  'overlap:prot-to-ref-overlap', 2, 2),
        ('SLF34.2',  'overlap:fraction-above-thresh-prot-in-above-thresh-ref', 2, 2),
        ('SLF34.3',  'overlap:fraction-of-protein-in-above-thresh-ref', 2, 2),
        ('SLF34.4',  'overlap:fraction-of-proc-protein-in-above-thresh-ref', 2, 2),
        ('SLF34.5',  'overlap:fraction-above-thresh-ref-in-above-thresh-prot', 2, 2),
        ('SLF34.6',  'overlap:correlation:binprot-binref', 2, 2),
        ('SLF34.7',  'overlap:correlation:prot-binref', 2, 2),
        ('SLF34.8',  'overlap:correlation:prot-ref', 2, 2),
        ('SLF34.9',  'overlap:median-prot-dist-ref', 2, 2),
        ('SLF34.10', 'overlap:mean-prot-dist-ref', 2, 2),
        ]

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
