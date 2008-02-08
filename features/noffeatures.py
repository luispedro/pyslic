def noffeatures(procimg,nofimg):
    obj_fluor = procimg.sum()
    nonobj_fluor = nofimg.sum()
    nonobj_feat = nonobj_fluor / (obj_fluor + nonobj_fluor)
    #feat_slf = [feat_slf cellstr('SLF7.79')];
    #feat_names = [feat_names cellstr('fract_nonobj_fluor')];
    return nonobj_feat


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
