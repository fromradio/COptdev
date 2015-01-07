// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Songtao Guo <guost@mathu.cn>
// Copyright (C) 2015 MathU
//			Reviewed by Ruimin Wang <ruimin.wang13@gmail.com>, <wangrm@mathu.cn>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.



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