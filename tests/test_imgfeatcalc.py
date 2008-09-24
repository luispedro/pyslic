import numpy
import pyslic
from os.path import dirname
basedir=dirname(__file__)

def test_imgfeatcalc():
    Img=pyslic.Image()
    Img.channels[Img.protein_channel]=basedir+'/data/protimg.bmp'
    Img.channels[Img.dna_channel]=basedir+'/data/dnaimg.bmp'
    F=pyslic.computefeatures(Img,'SLF7dna')
    assert F.size == 90
