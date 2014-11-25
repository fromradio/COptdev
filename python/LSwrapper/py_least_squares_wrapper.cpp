//		Copyright (C) Songtao Guo
//		Copyright (C) MathU
//			Reviewed by Ruimin Wang, ruimin.wang13@gmail.com, wangrm@mathu.cn


#include <Python.h>
#include <Header>

typedef double                          FT;
typedef COPT::Array<FT>                 Array;
typedef COPT::VectorBase<FT>            Vector;
typedef COPT::MatrixBase<FT>            Matrix;
typedef COPT::LeastSquaresSolver<FT>    LeastSquares;

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
static PyObject* pyLeastMeanSquare(PyObject *self,PyObject *args)
{

    PyObject *A,*b;
    float mu;
	if (!PyArg_ParseTuple(args,"OOf",&A,&b,&mu))
	 	return NULL;
	int m = PyObject_Size(A);
	int n = PyObject_Size(PyList_GetItem(A,0));
	Matrix A0(m,n);
	for(int i = 0;i < m;i++)
		for(int j = 0;j < n;j++)
			A0(i,j) = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(A,i),j));
	Vector b0(m);
	for(int i = 0;i < m;i++)
		b0[i] = PyFloat_AsDouble(PyList_GetItem(b,i));
	Vector x(n);

	// call COPT interface
	LeastSquares ls(A0,b0,mu);
	ls.setType(LeastSquares::LMS);
	ls.solve(x);
	x = ls.result();

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
static PyObject* pyLeastSquares(PyObject *self,PyObject *args)
{

    PyObject *A,*b;
	if (!PyArg_ParseTuple(args,"OO",&A,&b))
	 	return NULL;
	int m = PyObject_Size(A);
	int n = PyObject_Size(PyList_GetItem(A,0));
    Matrix A0(m,n);
	for(int i = 0;i < m;i++)
		for(int j = 0;j < n;j++)
			A0(i,j) = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(A,i),j));
	Vector b0(m);
	for(int i = 0;i < m;i++)
		b0[i] = PyFloat_AsDouble(PyList_GetItem(b,i));
	Vector x(n);

	// call COPT interface!
	LeastSquares ls(A0,b0);
	ls.setType(LeastSquares::LS);
	ls.solve(x);
	x = ls.result();

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
static PyObject* pyRecursiveLeastSquare(PyObject *self,PyObject *args)
{

    PyObject *A,*b;
    float lam,delta;
	if (!PyArg_ParseTuple(args,"OOff",&A,&b,&lam,&delta))
	 	return NULL;
	int m = PyObject_Size(A);
	int n = PyObject_Size(PyList_GetItem(A,0));
	Matrix A0(m,n);
	for(int i = 0;i < m;i++)
		for(int j = 0;j < n;j++)
			A0(i,j) = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(A,i),j));
	Vector b0(m);
	for(int i = 0;i < m;i++)
		b0[i] = PyFloat_AsDouble(PyList_GetItem(b,i));
	Vector x(n);

	// call COPT interface
	LeastSquares ls(A0,b0,0.01,lam,delta);
	ls.setType(LeastSquares::RLS);
	ls.solve(x);
	x = ls.result();

	PyObject *list;
	list = PyList_New(n);
	for(int j = 0;j < n;j++)
		PyList_SetItem(list,j,Py_BuildValue("f",x[j]));
	return Py_BuildValue("O",list);
}

static PyMethodDef leastSquaresMethods[]=
{
    {"lms",pyLeastMeanSquare,METH_VARARGS,"Least Mean Square Method"},
    {"ls",pyLeastSquares,METH_VARARGS,"Least Square Method"},
    {"rls",pyRecursiveLeastSquare,METH_VARARGS,"Recursive Least Square Method"},
    {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcopt()
{
	Py_InitModule("copt",leastSquaresMethods);
}