#include <Python.h>
#include <iostream>
#include <Header>

typedef float                      FT;
typedef COPT::Array<FT>             Array;
typedef COPT::VectorBase<FT>        Vector;
typedef COPT::MatrixBase<FT>        Matrix;


static PyObject* LADMethod(PyObject *self,PyObject *args)
{
	PyObject *A,*b,*u,*z;
	float rho;
	int number;
	float e;
	if (!PyArg_ParseTuple(args, "OOfOOif", &A,&b,&rho,&u,&z,&number,&e))
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
	Vector u0(m);
	for(int i = 0;i < m;i++)
		u0[i] = PyFloat_AsDouble(PyList_GetItem(u,i));
	Vector z0(m);
	for(int i = 0;i < m;i++)
		z0[i] = PyFloat_AsDouble(PyList_GetItem(z,i));
	Vector x(n);
	COPT::LeastAbsoluteDeviationMethod(A0,b0,rho,u0,x,z0,number,e);
	PyObject *list;
	list = PyList_New(n);
	for(int j = 0;j < n;j++)
		PyList_SetItem(list,j,Py_BuildValue("f",x[j]));
	std::cout<<number<<std::endl;
	return Py_BuildValue("O",list);
}

static PyMethodDef callMethods[]={
	{ "LADM", LADMethod, METH_VARARGS, "Least Absolute Deviation Method." },
	{ NULL, NULL, 0, NULL }	
};


PyMODINIT_FUNC initcall()
{
	Py_InitModule("call",callMethods);
}