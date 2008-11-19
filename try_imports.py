def try_imports():
    try:
        import numpy
    except:
        print '''Numpy not present!

    Please install numpy.
    Under linux, the package is often called python-numpy.'''

    try:
        import scipy
    except:
        print '''Scipy not present!

    Please install scipy.
    Under linux, the package is often called python-scipy.'''

    try:
        import Image, ImageDraw
    except:
        print '''Python Imaging Library (PIL) not present!

    Please install PIL.'''



    try:
        import scipy.weave
    except:
        print '''Scipy.weave not present!

    Please install scipy.weave.

    This is not strictly necessary, but it will make feature calculation faster (by 2~3 orders of magnitude).'''
    
    try:
        scipy.weave.inline('')
    except:
        print '''Scipy.weave.inline() failed.

    This is not strictly necessary, but it will make feature calculation faster (by 2~3 orders of magnitude).'''

    try:
        import ncreduce
    except:
        print '''import ncreduce failed.

    Ncreduce is not necessary, but makes classification and clustering faster.'''

    try:
        import svm
    except:
        print '''import svm failed.

    LibSVM classification will not be available (everything else will work).'''

if __name__ == '__main__':
    test_imports()

