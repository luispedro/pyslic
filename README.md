# PySLIC: Python Subcellular Location Image Classifier

This is a re-implementation in Python of the Murphy Lab image classification
code.

# INSTALL

    python setup.py install

or

    python setup.py install --prefix=<SOME_PATH_TO_INSTALL_TO>


# DEPENDENCIES

## Hard Dependencies

These are needed for Python 2.5 to run.

The dependencies are fairly standard libraries:

- numpy
- scipy
- [mahotas](http://mahotas.readthedocs.org)

On ubuntu, you can use the following command to install dependencies:

    sudo aptitude install python-numpy python-scipy

# Citation

PySLIC does not have a paper associated with it, but is a wrapper around
[mahotas](http://mahotas.readthedocs.org). Thus, if you use this in a
scientific publication, please cite:

**Luis Pedro Coelho** Mahotas: Open source software for scriptable computer
vision in Journal of Open Research Software, vol 1, 2013.
[DOI](http://dx.doi.org/10.5334/jors.ac)

