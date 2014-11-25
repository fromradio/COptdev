#include <Python.h>
#include <Header>

typedef double                      FT;
typedef COPT::Array<FT>             Array;
typedef COPT::VectorBase<FT>        Vector;
typedef COPT::MatrixBase<FT>        Matrix;


/*
		A wrapper of Least Squares methods for Python
*/

//		Least mean square method
//		args:
//		A:			the input matrix
//		b:			the right hand vector
//		mu:			the scaling factor
//		return value: the obtained coefficient
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
	COPT::LeastMeanSquareMethod(A0,b0,mu,x);

	PyObject *list;
	list = PyList_New(n);
	for(int j = 0;j < n;j++)
		PyList_SetItem(list,j,Py_BuildValue("f",x[j]));
	return Py_BuildValue("O",list);
}

static PyObject* geLS(PyObject *self,PyObject *args)
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
	//std::cout<<A0<<std::endl;
	Vector b0(m);
	for(int i = 0;i < m;i++)
		b0[i] = PyFloat_AsDouble(PyList_GetItem(b,i));
	Vector x(n);
	COPT::LeastSquareMethod(A0,b0,x);
	//std::cout<<x<<std::endl;

	PyObject *list;
	list = PyList_New(n);
	for(int j = 0;j < n;j++)
		PyList_SetItem(list,j,Py_BuildValue("f",x[j]));
	//std::cout<<list<<std::endl;
	return Py_BuildValue("O",list);
	//return Py_BuildValue("s","successful extension!");
}

static PyObject* geRLS(PyObject *self,PyObject *args)
{

    PyObject *A,*b;
    float lam,delta;
	if (!PyArg_ParseTuple(args,"OOff",&A,&b,&lam,&delta))
	 	return NULL;
	//std::cout<<PyObject_Size(A)<<std::endl;
    //std::cout<<PyObject_Size(PyList_GetItem(A,1))<<std::endl;
	int m = PyObject_Size(A);
	int n = PyObject_Size(PyList_GetItem(A,0));
    //std::cout<<PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(A,1),1))<<std::endl;
	Matrix A0(m,n);
	for(int i = 0;i < m;i++)
		for(int j = 0;j < n;j++)
			A0(i,j) = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(A,i),j));
	//std::cout<<A0<<std::endl;
	Vector b0(m);
	for(int i = 0;i < m;i++)
		b0[i] = PyFloat_AsDouble(PyList_GetItem(b,i));
	Vector x(n);
	COPT::RLS_Method(A0,b0,x,lam,delta);
	//std::cout<<x<<std::endl;

	PyObject *list;
	list = PyList_New(n);
	for(int j = 0;j < n;j++)
		PyList_SetItem(list,j,Py_BuildValue("f",x[j]));
	//std::cout<<list<<std::endl;
	return Py_BuildValue("O",list);
	//return Py_BuildValue("s","successful extension!");
}

static PyMethodDef geLSsMethods[]=
{
    {"geLMS",geLMS,METH_VARARGS,"Least Mean Square Method"},
    {"geLS",geLS,METH_VARARGS,"Least Square Method"},
    {"geRLS",geRLS,METH_VARARGS,"Recursive Least Square Method"},
    {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initgeLSs()
{
	Py_InitModule("geLSs",geLSsMethods);
}