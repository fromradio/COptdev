#include <Python.h>
#include <iostream>
#include <Header>

static PyObject* geVec(PyObject *self, PyObject *args)
{
	int i = 0,j = 0;
	PyObject *list;
	if (!PyArg_ParseTuple(args, "i", &i))
		return NULL;
	COPT::VectorBase<double> vec(i);
    list = PyList_New(i);
    for (j = 0;j < i;j++)
    	PyList_SetItem(list,j,Py_BuildValue("f",vec[j]));
    return Py_BuildValue("O", list);
}

static PyMethodDef geVecMethods[] =
{
	{ "geVec", geVec, METH_VARARGS, "generate a vector." },
	{ NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initgeVec()
{
	Py_InitModule("geVec", geVecMethods);
}
