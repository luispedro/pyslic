class Image(object):
    dna_channel='dna'
    protein_channel='protein'
    autofluorescence_channel='autofluorescence'
    __slots__ = ['label','features','regions','channels']
    def __init__(self):
        self.label=''
        self.features=None
        self.regions=None
        self.channels={}

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
