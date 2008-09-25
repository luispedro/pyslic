import numpy
import pyslic
from os.path import dirname
basedir=dirname(__file__)

def test_imgfeatcalc():
    Img=pyslic.Image()
    Img.channels[Img.protein_channel]=basedir+'/data/protimg.jp2'
    Img.channels[Img.dna_channel]=basedir+'/data/dnaimg.jp2'
    F=pyslic.computefeatures(Img,'SLF7dna')
    assert F.size == 90
