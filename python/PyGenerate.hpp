//		Copyright (C) Songtao Guo, guost@mathu.cn
//		Copyright (C) MathU

//    Generate PyObject* to Matrix or Vector

#ifndef PY_GENERATE_H
#define PY_GENERATE_H

#include <Python.h>

namespace COPT
{

template<class Matrix>
void PyGenerateMatrix(PyObject *A,Matrix& A0)
{
	int m = PyObject_Size(A);
	int n = PyObject_Size(PyList_GetItem(A,0));
    A0.resize(m,n);
	for(int i = 0;i < m;i++)
		for(int j = 0;j < n;j++)
			A0(i,j) = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(A,i),j));
}

template <class Vector>
void PyGenerateVector(PyObject *b,Vector& b0)
{
	int m = PyObject_Size(b);
	b0.resize(m);
	for(int i = 0;i < m;i++)
		b0[i] = PyFloat_AsDouble(PyList_GetItem(b,i));
}

}

#endif