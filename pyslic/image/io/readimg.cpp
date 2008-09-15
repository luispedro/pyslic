// Copyright (C) 2008  Murphy Lab
// Carnegie Mellon University
// 
// Written by Luís Pedro Coelho <lpc@cmu.edu>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
//
// For additional information visit http://murphylab.web.cmu.edu or
// send email to murphy@cmu.edu


#include <string>
#include <Magick++.h>


extern "C" {
    #include <Python.h>
    #include <numpy/ndarrayobject.h>
}

using namespace Magick;
namespace {
PyObject* array_from_image(Magick::Image& img) {
    PyArrayObject* ret = 0;
    try {
        Geometry size = img.size();
        unsigned w = size.width();
        unsigned h = size.height();

        int dimensions[3];
        dimensions[0] = h;
        dimensions[1] = w;
        dimensions[2] = 3;

        bool colourimage = (img.type() != GrayscaleType);
        int type;
        StorageType pixeltype;
        if (QuantumDepth == 8 || img.depth() == 8 || img.depth() == 1) {
            type = PyArray_UBYTE;
            pixeltype = CharPixel;
        } else if (QuantumDepth == 16) {
            type = PyArray_USHORT;
            pixeltype = ShortPixel;
        } else {
            PyErr_SetString(PyExc_EOFError,"Magick++ issue: Quantum depth is neither 8 nor 16.\nDon't know how to handle that.");
            return 0;
        }
        ret = (PyArrayObject*)PyArray_FromDims(2+colourimage,dimensions,type);
        if (!ret) {
            PyErr_SetString(PyExc_MemoryError,"Out of Memory");
            return 0;
        }
        const char* write_format = (colourimage ? "RGB" : "I");
        img.write(0,0,w,h,write_format,pixeltype,ret->data);

        return PyArray_Return(ret);
    } catch ( std::exception& error_ ) {
        PyErr_SetString(PyExc_EOFError,error_.what());
        Py_DECREF(ret);
        return 0;
    }
}
PyObject* readimg(PyObject* self, PyObject* args) {
	PyObject* input;
	if (!PyArg_ParseTuple(args,"O",&input) || !PyString_Check(input) || PyString_Size(input) == 0) {
         PyErr_SetString(PyExc_TypeError,"readimg takes a filename as input");
         return 0;
    }
    try {
        std::string fname = PyString_AsString(input);
        Image img;
        try {
            img.read(fname);
        } catch ( WarningCoder &warning ) {
            // Issuing warnings turned the program too verbose
            //
            // int status = PyErr_WarnEx( 0, ( std::string("ImageMagick Warning: ") + warning.what() ).c_str(), 1 );
            // if (status < 0) return 0;
        }
        return array_from_image(img);
    } catch ( std::exception& error_ ) {
        PyErr_SetString(PyExc_EOFError,error_.what());
        return 0;
    }
}
PyObject* readimgfromblob(PyObject* self, PyObject* args) {
	PyObject* input;
	if (!PyArg_ParseTuple(args,"O",&input) || !PyString_Check(input) || PyString_Size(input) == 0) {
         PyErr_SetString(PyExc_TypeError,"readimgfromblob takes image data as input");
         return 0;
    }
    try {
        Image img;
        try {
            img.read(Blob(PyString_AsString(input),PyString_Size(input)));
        } catch ( WarningCoder &warning ) {
            // Issuing warnings turned the program too verbose
            //
            //int status = PyErr_WarnEx( 0, ( std::string("ImageMagick Warning: ") + warning.what() ).c_str(), 1 );
            //if (status < 0) return 0;
        }
        return array_from_image(img);
    } catch ( std::exception& error_ ) {
        PyErr_SetString(PyExc_EOFError,error_.what());
        return 0;
    }
}
        

const char * readimg_doc = 
    "img = readimg(fname)\n"
    "\n"
    "Read an image using Image Magick.\n"
    "\n"
    "The returned array is either a uint8, or a uint16 array, depending on the image type.\n"
    "\n"
    "If the image is reported to be a colour image, then a 3-dimensional array is returned. "
    "The last dimension encodes the red, green, and blue channels.\n"
    "\n"
    "@see readimgfromblob\n"
    ;

const char * readimgfromblob_doc = 
    "img = readimgfromblob(data)\n"
    "\n"
    "Read an image using Image Magick.\n"
    "\n"
    "The returned array is either a uint8, or a uint16 array, depending on the image type.\n"
    "\n"
    "If the image is reported to be a colour image, then a 3-dimensional array is returned. "
    "The last dimension encodes the red, green, and blue channels.\n"
    "\n"
    "@see readimg\n"
    ;

PyMethodDef methods[] = {
  {"readimg",readimg, METH_VARARGS , readimg_doc },
  {"readimgfromblob",readimgfromblob, METH_VARARGS , readimgfromblob_doc },
  {NULL, NULL,0,NULL},
};

const char * module_doc = 
    "readimg: Read an image into a numpy array.\n"
    "\n"
    "The implementation of this module is based on ImageMagick. Therefore, it supports\n"
    "all image formats that ImageMagick supports.\n"
    "\n"
    "\n"
    "This module provides two functions:\n"
    "   - readimg(filename): reads an image from a file\n"
    "   - readimgfromblob(blob): reads an image from a blog (character string)\n"
    ;
}

extern "C"
void initreadimg()
  {
    import_array();
    (void)Py_InitModule3("readimg", methods, module_doc);
  }

