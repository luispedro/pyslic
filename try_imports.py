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

if __name__ == '__main__':
    try_imports()

