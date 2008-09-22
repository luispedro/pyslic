from __future__ import division
__all__ = ['noffeatures']
def noffeatures(procimg,nofimg):
    """
    Compute non object fluorescence features
    """
    obj_fluor = procimg.sum()
    nonobj_fluor = nofimg.sum()
    nonobj_feat = nonobj_fluor / (obj_fluor + nonobj_fluor)
    return nonobj_feat

noffeatures.names = 'fract_nonobj_fluor'
noffeatures.slf_names = ['SLF7.79']
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
