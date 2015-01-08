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


#include <Core>
#include <Python.h>
#include <IO>
#include <PyGenerate.hpp>

typedef double 		FT;
typedef COPT::KernelTrait<FT>		kernel;
typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::LeastSquaresProblem<kernel>	problem;
typedef COPT::LeastMeanSquareSolver<problem,COPT::SolverTimeStatistics>	lmssolver;
typedef COPT::LeastSquareSolver<problem,COPT::SolverTimeStatistics>	lssolver;
typedef COPT::RecursiveLeastSquareSolver<problem,COPT::SolverTimeStatistics> 	rlssolver;

/*
		A wrapper of Least Squares methods for Python
*/

/*		Least mean square method
 *		args:
 *		A:			the input matrix
 *		b:			the right hand vector
 *		mu:			the scaling factor
 *		return value: the obtained coefficient
 */
static PyObject* pyLeastMeanSquare(PyObject *self,PyObject *args,PyObject *kw)
{
    PyObject *A,*b;
    static char* kwlist[] = {"A","b","mu",NULL};
    float mu = 0.01;
	if (!PyArg_ParseTupleAndKeywords(args,kw,"OO|f",kwlist,&A,&b,&mu))
		return NULL;
	Matrix A0;
	Vector b0;
	PyGenerateMatrix(A,A0);
	PyGenerateVector(b,b0);
	int n = PyObject_Size(PyList_GetItem(A,0));
	Vector x(n);

	// call COPT interface
	problem pro(A0,b0);
	lmssolver lmssol(pro,mu);
	lmssol.solve();
	x = lmssol.result();

	PyObject *list;
	list = PyList_New(n);
	for(int j = 0;j < n;j++)
		PyList_SetItem(list,j,Py_BuildValue("f",x[j]));
	return Py_BuildValue("O",list);
}


/*		normal least square method
 *		args:
 *		A:			input matrix
 *		b:			right hand vector
 *		return value:	the obtained coefficient
 */
static PyObject* pyLeastSquare(PyObject *self,PyObject *args,PyObject *kw)
{

    PyObject *A,*b;
    static char* kwlist[] = {"A","b",NULL};
	if (!PyArg_ParseTupleAndKeywords(args,kw,"OO",kwlist,&A,&b))
	 	return NULL;
    Matrix A0;
	Vector b0;
	PyGenerateMatrix(A,A0);
	PyGenerateVector(b,b0);
	int n = PyObject_Size(PyList_GetItem(A,0));
	Vector x(n);

	// call COPT interface!
	problem pro(A0,b0);
	lssolver lssol(pro);
	lssol.solve();
	x = lssol.result();

	PyObject *list;
	list = PyList_New(n);
	for(int j = 0;j < n;j++)
		PyList_SetItem(list,j,Py_BuildValue("f",x[j]));
	return Py_BuildValue("O",list);
}

/*		normal least square method
 *		args:
 *		A:			input matrix
 *		b:			right hand vector
 *		lam:		parameter one
 *		delta:		the second parameter
 *		return value:	the obtained coefficient
 */
static PyObject* pyRecursiveLeastSquare(PyObject *self,PyObject *args,PyObject *kw)
{

    PyObject *A,*b;
    static char* kwlist[] = {"A","b","lam","delta",NULL};
    float lam = 1,delta = 250;
	if (!PyArg_ParseTupleAndKeywords(args,kw,"OO|ff",kwlist,&A,&b,&lam,&delta))
	 	return NULL;
	Matrix A0;
	Vector b0;
	PyGenerateMatrix(A,A0);
	PyGenerateVector(b,b0);
	int n = PyObject_Size(PyList_GetItem(A,0));
	Vector x(n);
 
	// call COPT interface
	problem pro(A0,b0);
	rlssolver rlssol(pro,lam,delta);
	rlssol.solve();
	x = rlssol.result();

	PyObject *list;
	list = PyList_New(n);
	for(int j = 0;j < n;j++)
		PyList_SetItem(list,j,Py_BuildValue("f",x[j]));
	return Py_BuildValue("O",list);
}

static PyMethodDef leastSquaresMethods[]=
{
    {"lms",(PyCFunction)pyLeastMeanSquare,METH_VARARGS|METH_KEYWORDS,"Least Mean Square Method"},
    {"ls",(PyCFunction)pyLeastSquare,METH_VARARGS|METH_KEYWORDS,"Least Square Method"},
    {"rls",(PyCFunction)pyRecursiveLeastSquare,METH_VARARGS|METH_KEYWORDS,"Recursive Least Square Method"},
    {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcopt()
{
	Py_InitModule("copt",leastSquaresMethods);
}