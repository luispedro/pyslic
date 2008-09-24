import setuptools
import numpy.distutils.core as numpyutils

convexhull = numpyutils.Extension('pyslic/imageprocessing/_convexhull', sources = ['pyslic/imageprocessing/convexhull.cpp'])

magick_libraries= ['Magick++','Wand','Magick'] # This is gotten from magick++-config --libs
readimg = numpyutils.Extension('pyslic.image.io.readimg', sources = ['pyslic/image/io/readimg.cpp'], libraries=magick_libraries) 
packages=setuptools.find_packages()
packages.remove('tests')

numpyutils.setup(name='PySLIC',
      version='0.3.1',
      description='Subcellular Location Image Classifier',
      author='Murphy Lab',
      author_email='murphy@mcu.edu',
      url='http://murphylab.cbi.cmu.edu/',
      packages=packages,
      ext_modules = [convexhull,readimg],
      test_suite='nose.collector',
      )


