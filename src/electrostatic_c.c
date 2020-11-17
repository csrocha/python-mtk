/*
 * Implementación de la electrostatica en C
 * $Id$
 */

#include <Python.h>
#include <numpy/arrayobject.h>
#include <string.h>
#include "vol_c.h"

struct atom_s {
	real_t x, y, z, vw, c;
};

typedef struct atom_s atom_t;

inline npy_double epsilon(npy_double r) {
	if (r <= 6.0) {
		return 4.0;
	} else
	if (r <= 8.0) {
		return 38.0*r - 224.0;
	} else {
		return 80.0;
	}
}

#define GET_ATTR_ARRAY(OBJ, ATTR) ( PyArray_DATA((PyArrayObject *)PyObject_GetAttrString(OBJ, ATTR) ) )

static PyObject *
electrostatic_calc(PyObject *self, PyObject *args)
{
	PyObject *vol = NULL;
	PyObject *array = NULL;
	real_t *shape, *_min, *delta;
	PyArrayObject *grid = NULL, *molecule = NULL;
	PyArrayIterObject *itr, *itr_atm;
	atom_t *atom;
	int_t  no_atoms, i;
	complex_t *v;
	real_t *distance;
	real_t x[3];
	real_t d[3];
	real_t r, charge;

	// Load python arguments
	if (!PyArg_ParseTuple(args, "OO!", &vol, &PyArray_Type, &array))
		return NULL;

	if (!PyObject_IsInstance(vol, (PyObject *)BasicVolume_Class))
		return NULL;

        if (PyArray_TYPE(array) != NPY_DOUBLE)
		return NULL;

	shape = GET_DATA_ATTR_REAL_IN(vol, "shape");
	_min  = GET_DATA_ATTR_REAL_IN(vol, "min");
	delta = GET_DATA_ATTR_REAL_IN(vol, "delta");
	grid  = GET_ATTR_COMPLEX_INOUT(vol, "_data");

	molecule = NPY_ARRAY_IN_TYPE(array, REAL_T);
	if (molecule == NULL)
		return NULL;
	if (PyArray_DIM(molecule,1) != 5) {
		PyErr_SetString(PyExc_TypeError, "The list of coordinates must be a (nx5) array, each column must contain ( x, y, z, vw, c ) real values.");
		return NULL;
		}

	no_atoms = PyArray_DIM(molecule, 0);
	distance = malloc(no_atoms * sizeof(real_t));

	// Create iterator over the volume

	itr = (PyArrayIterObject *)PyArray_IterNew((PyObject *)grid);
	assert(itr->nd_m1 == 2);

	// Create iterator over atoms of the molecule

	itr_atm = (PyArrayIterObject *)PyArray_IterNew((PyObject *)molecule);
	assert(itr_atm->nd_m1 == 1);

	// printf("\n");

	// Iterate over nodes of the grid

	PyArray_ITER_RESET(itr);
	while( PyArray_ITER_NOTDONE(itr) ) {
		v = (complex_t *) itr->dataptr;
		for (i = 0; i < 3; i++) x[i] = _min[i] + delta[i]*itr->coordinates[i];

		charge = 0;
		PyArray_ITER_RESET(itr_atm);

		// Sum atoms in radius r
		while( PyArray_ITER_NOTDONE(itr_atm) ) {
			atom = (atom_t *) itr_atm->dataptr;
			d[0] = x[0] - atom->x;
			d[1] = x[1] - atom->y;
			d[2] = x[2] - atom->z;
			for (i = 0; i < 3; i++) d[i] *= d[i];
			r = sqrt(d[0] + d[1] + d[2]);
			if (r <= atom->vw) { charge = 0; break; }
			charge += atom->c / (epsilon(r) * r);
			//printf("x: %+2.2f %+2.2f %+2.2f , a: %+2.2f %+2.2f %+2.2f, r: %+2.4f , e: %+2.4f, pc: %+2.4f, c: %+2.4f\n", x[0], x[1], x[2], atom->x, atom->y, atom->z, r, epsilon(r),  atom->c / (epsilon(r) * r), charge);
			_PyAIT(itr_atm)->index+=5;
			_PyArray_ITER_NEXT1(itr_atm);
		}

		// Store the charge
		v->real = charge;
		v->imag = 0;

		_PyAIT(itr)->index++;
		_PyArray_ITER_NEXT3(itr);
	}

	free(distance);
	Py_DECREF(itr);
	Py_DECREF(itr_atm);
	Py_DECREF(molecule);

	Py_RETURN_NONE;
}

/*
 * Declara el calculo de la rotacion
 */
static PyMethodDef electrostatic_methods[] = {
    {"calc",  electrostatic_calc, METH_VARARGS,
     "Calcule electrostatic in a grid."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

/* 
 * Iniciación del modulo.
 */

PyMODINIT_FUNC
initelectrostatic_c(void)
{
	PyObject *m;

	m = Py_InitModule("electrostatic_c", electrostatic_methods);
	import_array();
	import_BasicVolume();
}

