from .features import computefeatures
import image
from image import Image
try:
    from readmagick import readimg
except:
    from scipy.misc.pilutil import imread
    readimg = imread
    del imread
import preprocess
import clustering
import utils
import segmentation
from .preprocess import preprocessimage
from .imageprocessing import thresholding
__all__ = [
    'image',
    'Image',
    'preprocess',
    'computefeatures',
    'clustering',
    'utils',
    'segmentation',
    'thresholding',
    ]
