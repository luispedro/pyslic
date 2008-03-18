#include <Python.h>
#include <algorithm>
#include <iostream>
namespace {
struct Point {
	long x, y;
};

inline
bool forward_cmp(const Point& a, const Point& b) {
	if (a.x == b.x) return a.y < b.y;
	return a.x < b.x;
}
inline
bool reverse_cmp(const Point& a, const Point& b) {
	if (a.x == b.x) return a.y > b.y;
	return a.x > b.x;
}

inline
double isLeft(Point p0,Point p1, Point p2) { 
	return (p2.x-p0.x)*(p1.y-p0.y) - (p1.x-p0.x)*(p2.y-p0.y);
}

unsigned inPlaceScan(Point* P, unsigned N, bool reverse) {
	if (reverse) {
		std::sort(P,P+N,forward_cmp);
	} else {
		std::sort(P,P+N,reverse_cmp);
	}
	int h = 1;
	for (int i = 1; i != N; ++i) {
		while (h >= 2 && isLeft(P[h-2],P[h-1],P[i]) >= 0) {
			--h;
		}
		std::swap(P[h],P[i]);
		++h;
	}
	return h;
}

unsigned inPlaceGraham(Point* P, unsigned N) {
	int h = inPlaceScan(P,N,false);
	for (int i = 0; i != h - 1; ++i) {
		std::swap(P[i],P[i+1]);
	}
	int h2=inPlaceScan(P+h-2,N-h+2,true);
	return h + h2;
}

PyObject*
compute_convexHull(PyObject* self, PyObject* args) {
	PyObject* input;
	if (!PyArg_ParseTuple(args,"O",&input)) return 0;
	const Py_ssize_t N = PyList_Size(input);
    if (N <= 3) {
        // The python wrapper should have stopped this from happening
        Py_RETURN_NONE; 
        //    Py_INCREF(input);
        //    return input;
    }
	Point* P = new Point[N];
	for (Py_ssize_t i = 0; i != N; ++i) {
		PyObject* tup = PyList_GetItem(input,i);
		long x = PyInt_AsLong(PyTuple_GetItem(tup,0));
		long y = PyInt_AsLong(PyTuple_GetItem(tup,1));
		P[i].x=x;
		P[i].y=y;
	}
	unsigned h = inPlaceGraham(P,N);
	PyObject* output = PyList_New(h);
	if (!output) {
		PyErr_NoMemory();
		return 0;
	}
	for (unsigned i = 0; i != h; ++i) {
		PyObject* tup = PyTuple_New(2);
		PyTuple_SetItem(tup,0,PyInt_FromLong(P[i].x));
		PyTuple_SetItem(tup,1,PyInt_FromLong(P[i].y));
		PyList_SetItem(output,i,tup);
	}
	delete [] P;
	return output;
}


}

PyMethodDef methods[] = {
  {"computeConvexHull",compute_convexHull, METH_VARARGS , "compute convex hull"},
  {NULL, NULL,0,NULL},
};

extern "C"
void init_convexhull()
  {
    (void)Py_InitModule("_convexhull", methods);
  }

