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
    try {
        Geometry size = img.size();
        unsigned w = size.width();
        unsigned h = size.height();

        int dimensions[3];
        dimensions[0] = h;
        dimensions[1] = w;
        dimensions[2] = 3;
        const PixelPacket* pixels = img.getConstPixels(0,0,w,h);
        PyArrayObject* ret = (PyArrayObject*)PyArray_FromDims(3,dimensions,PyArray_USHORT);
        if (!ret) {
            PyErr_SetString(PyExc_MemoryError,"Out of Memory");
            return 0;
        }
        for (int i = 0; i != h; ++i) {    
            for (int j = 0; j != w; ++j) {    
                const PixelPacket& p = pixels[i*w+j];
                *reinterpret_cast<unsigned short*>(ret->data + (i*ret->strides[0]+j*ret->strides[1]+0*ret->strides[2])) = p.red;
                *reinterpret_cast<unsigned short*>(ret->data + (i*ret->strides[0]+j*ret->strides[1]+1*ret->strides[2])) = p.green;
                *reinterpret_cast<unsigned short*>(ret->data + (i*ret->strides[0]+j*ret->strides[1]+2*ret->strides[2])) = p.blue;
            }
        }
        return PyArray_Return(ret);
    } catch ( std::exception& error_ ) {
        PyErr_SetString(PyExc_EOFError,error_.what());
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
        img.read(fname);
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
        Image img(Blob(PyString_AsString(input),PyString_Size(input)));
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
    "Always returns a 3-colour image.\n";

const char * readimgfromblob_doc = 
    "img = readimgfromblob(data)\n"
    "\n"
    "Read an image using Image Magick.\n"
    "Always returns a 3-colour image.\n";

PyMethodDef methods[] = {
  {"readimg",readimg, METH_VARARGS , readimg_doc },
  {"readimgfromblob",readimgfromblob, METH_VARARGS , readimg_doc },
  {NULL, NULL,0,NULL},
};
}

extern "C"
void initreadimg()
  {
    import_array();
    (void)Py_InitModule("readimg", methods);
  }

