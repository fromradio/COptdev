//		Copyright (C) Songtao Guo
//		Copyright (C) MathU
//			Reviewed by Ruimin Wang, ruimin.wang13@gmail.com, wangrm@mathu.cn


#include <Header>
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
static PyObject* pyLeastMeanSquare(PyObject *self,PyObject *args)
{

    PyObject *A,*b;
    float mu;
	if (!PyArg_ParseTuple(args,"OOf",&A,&b,&mu))
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
static PyObject* pyLeastSquare(PyObject *self,PyObject *args)
{

    PyObject *A,*b;
	if (!PyArg_ParseTuple(args,"OO",&A,&b))
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
static PyObject* pyRecursiveLeastSquare(PyObject *self,PyObject *args)
{

    PyObject *A,*b;
    float lam,delta;
	if (!PyArg_ParseTuple(args,"OOff",&A,&b,&lam,&delta))
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
    {"lms",pyLeastMeanSquare,METH_VARARGS,"Least Mean Square Method"},
    {"ls",pyLeastSquare,METH_VARARGS,"Least Square Method"},
    {"rls",pyRecursiveLeastSquare,METH_VARARGS,"Recursive Least Square Method"},
    {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcopt()
{
	Py_InitModule("copt",leastSquaresMethods);
}