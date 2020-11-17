/*
 * Implementación de la electrostatica en C
 * $Id$
 */

#include "Python.h"
#include <numpy/arrayobject.h>
#include <string.h>

static PyObject *Volume_Class;

static PyObject *
rotation_calc(PyObject *self, PyObject *args)
{
	PyObject *vol_orig = NULL;
	PyObject *vol_dest = NULL;
	PyObject *rotation = NULL;

	if (!PyArg_ParseTuple(args, "OOO!", &vol_orig, &vol_dest, &PyArray_Type, &rotation))
		return NULL; 

	if (!PyObject_IsInstance(vol_orig, Volume_Class) || !PyObject_IsInstance(vol_dest, Volume_Class))
		return NULL;

        if (PyArray_TYPE(rotation) != NPY_DOUBLE || PyArray_)
		return NULL;

}

/*
 * Declara el calculo de la rotacion
 */

static PyMethodDef rotation_methods[] = {
    {"calc",  rotation_calc, METH_VARARGS,
     "Calcule electrostatic in a grid."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

/* 
 * Iniciación del modulo.
 */

PyMODINIT_FUNC
initelectrostatic_c(void)
{
	PyObject *m, *__dict__, *s;
	PyObject *modGeoVol;

	m = Py_InitModule("rotation_c", rotation_methods);
	import_array();

	/* Importa el volumen y se queda con la clase Volume */
	modGeoVol = PyImport_ImportModule("mtk.geometry.vol");
	__dict__ = PyModule_GetDict(modGeoVol);
	Volume_Class = PyDict_GetItemString(__dict__, "Volume");
	if (!PyClass_Check(Volume_Class)) {
		/* TODO: Must raise an error import exception */
		printf("No es una clase!!\n");
	}
	s = PyObject_Str(Volume_Class);
	Py_INCREF(Volume_Class);
}

