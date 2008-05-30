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
import random

__all__ = [
    'get_random',
    'get_pyrandom',
    'format_table',
    'format_confusion_matrix',
    ]

def get_random(R):
    '''
    R = get_random(R)

    Returns a numpy.RandomState from R
    @param R can be one of:
        * None          : Returns the default numpy global state
        * integer       : Uses it as a seed for constructing a new random generator
        * RandomState   : returns R
    '''
    if R is None:
        return numpy.random.mtrand._rand
    if type(R) == int:
        return numpy.random.RandomState(R)
    return R

def get_pyrandom(R):
    '''
    R = get_pyrandom(R)

    Returns a random.Random object based on R

    @param R can be one of:
        * None          : Returns the default python Random object
        * integer       : Uses it as seed for constructing a new random generator
        * RandomState   : Uses it as to generate a seed for a new random generator
    '''
    if R is None:
        return random
    if type(R) is int:
        return random.Random(R)
    if type(R) is numpy.random.RandomState:
        return random.Random(R.randint(2**30))
    raise TypeError,"get_pyrandom() does not know how to handle type %s." % type(R)


def format_table(table,collabels,rowlabels,format='latex'):
    '''
    format_table(table,collabels,rowlabels,format='latex')

    Format a table.

    @param format: Currently support "latex"
    '''
    if format != 'latex':
        raise AttributeError, 'format_table: only \'latex\' format supported'

    table=numpy.asanyarray(table)
    r,c=table.shape

    eol="\\\\\n"
    heading=' & ' + ' & '.join(collabels) + eol
    content=''
    for i in xrange(r):
        header = (str(rowlabels[i]) if rowlabels else '')
        row = header + ' & ' + ' & '.join(str(elem) for elem in table[i]) + eol
        content += row

    return '''
    \\begin{tabular}
    \\toprule
        %(heading)s
    \\midrule
        %(content)s
    \\bottomrule
    ''' % locals()

def format_confusion_matrix(cmatrix,labels,format='latex'):
    '''
    @see format_table
    '''
    return format_table(cmatrix,labels,labels,format)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
