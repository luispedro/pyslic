from __future__ import division
from numpy import *
from scipy.ndimage import label
from ..imageprocessing.bwperim import bwperim
from ..imageprocessing.bbox import croptobbox
import warnings

__all__ = ['prefilter']
def prefilter(dnamasks,scale,minarea,maxarea,minroundness):
    '''
    labeled, N = prefilter(dnamasks,scale,minarea,maxarea,minroundness)

    Performs a filtering version of label(), where objects are kept
        only if they 
    '''

    labeled,N=label(dnamasks)

    positives=zeros(N)
    minarea /= scale
    maxarea /= scale
    for obj in xrange(1,N+1):
        objimg = croptobbox(labeled == obj)
        area=objimg.sum()
        print area, maxarea, minarea
        if area > maxarea or area < minarea:
            continue
        perim=bwperim(objimg).sum()
        roundness=perim**2/(4*pi*area)
        print roundness
        if roundness < minroundness:
            continue
        positives[obj-1]=1
    try:
        from scipy import weave
        from scipy.weave import converters
        D1,D2=labeled.shape
        code = '''
        for (int y = 0; y != D1; ++y) {
            for (int x = 0; x != D2; ++x) {
                if (labeled(y,x)) {
                    if (!positives(labeled(y,x)-1)) labeled(y,x) = 0;
                }
            }  
        }'''
        weave.inline(code,
            ['labeled','positives','D1','D2'],
            type_converters=converters.blitz)
    except Exception, e:
        warnings.warn('scipy.weave failed. Resorting to (slow) Python code (Error: %s)' % e)
        for obj in xrange(N):
            if not positives[obj]:
                labeled[labeled == (obj+1)]=0
    return labeled,int(positives.sum())
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
