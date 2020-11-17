/*
 * Implementación del volumen en C
 * $Id$
 */
#include <Python.h>
#include "vol_c.h"

/* Structure to control de value assigment in put_ball function */
typedef struct {
	real_t radius;
	complex_t value;
} radius_map_t;

/* Copy Dictionary to radius_map_t */
static radius_map_t *
create_radius_map(PyObject *dict, Py_ssize_t *outsize) {
	PyObject *key, *value;
	Py_ssize_t size;
	Py_ssize_t pos = 0;
	radius_map_t *outmap;
	radius_map_t *item;

	// Check if input is a dictionary
	assert(PyDict_Check(dict));

	// Create map structure
	size = PyDict_Size(dict);
	outmap = malloc(size * sizeof(radius_map_t));
	if (outmap == NULL) return (radius_map_t*)PyErr_NoMemory();

	item = outmap;
	while (PyDict_Next(dict, &pos, &key, &value)) {
		item->radius = PyFloat_AsDouble(key);
		item->value.real = PyComplex_RealAsDouble(value);
		item->value.imag = PyComplex_ImagAsDouble(value);
		item++;
	}

	*outsize = size;
	return outmap;
}

/* Sort radius map */
static void
sort_radius_map(radius_map_t *map, const Py_ssize_t size) {
	int i, j;
	radius_map_t tmp;

	for (i = 0; i < size; i++) {
		for (j = i+1; j < size; j++) {
			if (map[j].radius < map[i].radius) {
				// tmp := map[i]
				tmp.radius = map[i].radius;
				tmp.value.real = map[i].value.real;
				tmp.value.imag = map[i].value.imag;
				// map[i] := map[j]
				map[i].radius = map[j].radius;
				map[i].value.real = map[j].value.real;
				map[i].value.imag = map[j].value.imag;
				// map[j] := tmp
				map[j].radius = tmp.radius;
				map[j].value.real = tmp.value.real;
				map[j].value.imag = tmp.value.imag;
			}
		}
	}
}

/* Get value in sorted map 
 *
 * vdw    9j
 * vdw+e  1
 * r      0   
 *
 *  0, 1 -> 1
 *  0,9j -> 9j
 *  1, 1 -> 1
 *  1,9j -> 9j
 * 9j, 1 -> 9j
 * 9j,9j -> 9j
 *
 * 9j -> 9j
 *  1 -> 
 *
 * */
static complex_t*
get_value_in_map(const real_t radius, radius_map_t *map, const Py_ssize_t size) {
	int i;

	for (i = 0; i < size; i++) {
		if (radius <= map[i].radius) return &(map[i].value);
	}

	return NULL;
}

/*
 * Construye una tabla más eficiente para checkear las
 * prioridades de los valores en el volumen.
 */
static complex_t*
create_order(PyObject *obj, Py_ssize_t *size) {
	complex_t *order, *tmp;
	PyObject *iter, *item;

	*size = PySequence_Size(obj);
	order = (complex_t*)malloc(*size * sizeof(complex_t));
	if (order == NULL) return (complex_t*)PyErr_NoMemory();

	iter = PyObject_GetIter(obj);

	tmp = order;
	while ((item = PyIter_Next(iter))) {
		tmp->real = PyComplex_RealAsDouble(item);
		tmp->imag = PyComplex_ImagAsDouble(item);
		tmp++;
		Py_DECREF(item);
	}
	Py_DECREF(iter);
	return order;
}

/*
 * La funcion is_substitute permite determinar que número complejo es más
 * importante en la asignación. Si A es más importante que B entonces
 * A no puede ser substituida.
 */
static int
is_substitute(complex_t *A, complex_t *B, complex_t *order, Py_ssize_t size) {
	int i;
	int Aorder = -1, Border = -1;

	for (i = 0; i < size; i++) {
		if (A->real == order[i].real && A->imag == order[i].imag) Aorder = i;
		if (B->real == order[i].real && B->imag == order[i].imag) Border = i;
	}
	return (Aorder > Border);
}

static PyObject *
Volume_put_ball(PyObject *self, PyObject *args)
{
	PyObject *ivol = NULL;
	PyObject *icoord = NULL;
	PyObject *imap = NULL;
	PyObject *iorder = NULL;
	PyArrayObject *np_coord, *np_Vmin, *np_Vdelta, *np_Vdata;
	complex_t *v, *sv, *order;
	real_t *Vmin, *Vdelta, *coord;
	npy_intp e_coord[3];
	real_t max_radius, radius, d[3];
	PyArrayIterObject *iter;
	PyArrayNeighborhoodIterObject *neigh_iter;
	int i, j;
	int debug=0;
	npy_intp bounds[6];
	radius_map_t *map;
	Py_ssize_t mapsize=0, ordersize;

	// Setting parameters
	if (!PyArg_ParseTuple(args, "OOO!O|i",
				&ivol,
				&icoord,
				&PyDict_Type, &imap,
				&iorder,
				&debug))
		return NULL;

	if (!PyObject_IsInstance(ivol, (PyObject *)BasicVolume_Class)) {
		PyErr_SetString(PyExc_TypeError, "Arguments one is not mtk.geometry.vol.BasicVolume");
		return NULL;
	}

	if (!(PyArray_NDIM(icoord)==1 &&
	      PyArray_DIM(icoord,0) == 3)) {
		PyErr_SetString(PyExc_TypeError, "Argument two must be a coordinate, tree values.");
		return NULL;
	}

	if (!(PySequence_Check(iorder))) {
		PyErr_SetString(PyExc_TypeError, "Argument four must be a sequence object of complexes numbers.");
		return NULL;
	}
#if 1
	np_coord = NPY_ARRAY_IN_REAL(icoord);
	np_Vmin = GET_ARRAY_IN_REAL(ivol, "min");
	np_Vdelta = GET_ARRAY_IN_REAL(ivol, "delta");
	np_Vdata = GET_ARRAY_INOUT_COMPLEX(ivol, "_data");
	coord = (real_t *)PyArray_DATA(np_coord);
	Vmin = (real_t *)PyArray_DATA(np_Vmin);
	Vdelta = (real_t *)PyArray_DATA(np_Vdelta);
#else
	coord  = NPY_DATA_ARRAY_IN_REAL(icoord);
	Vmin   = GET_DATA_ATTR_REAL_IN(ivol, "min");
	Vdelta = GET_DATA_ATTR_REAL_IN(ivol, "delta");
	np_Vdata  = GET_ATTR_COMPLEX_INOUT(ivol, "_data");
#endif

	map = create_radius_map(imap, &mapsize);
	if (map == NULL) return PyErr_NoMemory();
	sort_radius_map(map, mapsize);
	max_radius = map[mapsize-1].radius;

	order = create_order(iorder, &ordersize);
	if (order == NULL) return PyErr_NoMemory();

	for(i = 0; i < 3; i++) e_coord[i] = (coord[i] - Vmin[i]) / Vdelta[i];
	for(i = 0; i < 3; i++) {
		bounds[i*2] = -max_radius/Vdelta[i]-1;
		bounds[i*2+1] = max_radius/Vdelta[i]+1;
	}

#ifdef NDEBUG
#else
	if (debug) {
	printf("coord  :%.2f %.2f %.2f:\n", (float)coord[0], (float)coord[1], (float)coord[2]);
	printf("e_coord:%i %i %i:\n", (int)e_coord[0], (int)e_coord[1], (int)e_coord[2]);
	printf("Vmin   :%.2f %.2f %.2f:\n", (float)Vmin[0], (float)Vmin[1], (float)Vmin[2]);
	printf("Vdelta :%.2f %.2f %.2f:\n", (float)Vdelta[0], (float)Vdelta[1], (float)Vdelta[2]);
	printf("bounds :%i %i %i %i %i %i:\n", (int)bounds[0], (int)bounds[1], (int)bounds[2], (int)bounds[3], (int)bounds[4], (int)bounds[5]);
	}
#endif
	// Iterating over volume space where ball must be stored
	iter = (PyArrayIterObject *)PyArray_IterNew((PyObject*)np_Vdata);
	PyArray_ITER_GOTO(iter, e_coord);
	neigh_iter = (PyArrayNeighborhoodIterObject*)PyArray_NeighborhoodIterNew(
		      iter, bounds, NPY_NEIGHBORHOOD_ITER_CIRCULAR_PADDING, NULL);
	for (i = 0; i < neigh_iter->size; ++i) {
		v = (complex_t *) neigh_iter->dataptr;
		for (j = 0; j < 3; j++) d[j] = coord[j] - (Vmin[j] + (neigh_iter->coordinates[j] + iter->coordinates[j])*Vdelta[j]);
		radius = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
		sv = get_value_in_map(radius, map, mapsize);
		if (sv != NULL && is_substitute(sv, v, order, ordersize)) {
			v->real = sv->real;
			v->imag = sv->imag;
		}
		PyArrayNeighborhoodIter_Next(neigh_iter);
	}
	// Limpieza de memoria
	free(map);
	free(order);

	Py_DECREF(neigh_iter);
	Py_DECREF(iter);
	Py_DECREF(np_Vdata);
	Py_DECREF(np_Vdelta);
	Py_DECREF(np_Vmin);
	Py_DECREF(np_coord);

	Py_RETURN_NONE;
}

static PyObject *
Volume_transfer(PyObject *self, PyObject *args)
{
	/*
	 * Python function to transfer data from Volume VolA to Volume VolB
	 * with a transformation T.
	 */
	PyObject *volA = NULL;
	PyObject *volB = NULL;
	PyObject *transform = NULL;
	npy_intp i;
	complex_t *v;
	real_t x[4], u[4];
	real_t *Amin,*Ashape,*Asize,*Adelta;
	real_t *Bmin,*Bshape,*Bsize,*Bdelta;
	PyArrayObject *Adata, *Bdata, *T;
	PyArrayIterObject *iter;
	int debug=0;

	if (!PyArg_ParseTuple(args, "OOO!|i",
				&volA,
				&volB,
				&PyArray_Type, &transform,
				&debug))
		return NULL;

	if (!PyObject_IsInstance(volA, (PyObject *)BasicVolume_Class) &&
	    !PyObject_IsInstance(volB, (PyObject *)BasicVolume_Class)) {
		PyErr_SetString(PyExc_TypeError, "Arguments one or two are not mtk.geometry.vol.BasicVolume");
		return NULL;
	}

	if (!(PyArray_SIZE(transform)==4*4 &&
	      PyArray_DIM(transform,0) == 4 &&
	      PyArray_DIM(transform,1) == 4) ) {
		PyErr_SetString(PyExc_TypeError, "Arguments tree is not a afine transformation array (4x4)");
		return NULL;
	}

	T = NPY_ARRAY_IN_TYPE(transform, REAL_T);

	Amin   = GET_DATA_ATTR_REAL_IN(volA, "min");
	Ashape = GET_DATA_ATTR_REAL_IN(volA, "shape");
	Asize  = GET_DATA_ATTR_REAL_IN(volA, "size");
	Adelta = GET_DATA_ATTR_REAL_IN(volA, "delta");
	Adata  = GET_ATTR_COMPLEX_IN(volA, "_data");

	Bmin   = GET_DATA_ATTR_REAL_IN(volB, "min");
	Bshape = GET_DATA_ATTR_REAL_IN(volB, "shape");
	Bsize  = GET_DATA_ATTR_REAL_IN(volB, "size");
	Bdelta = GET_DATA_ATTR_REAL_IN(volB, "delta");
	Bdata  = GET_ATTR_COMPLEX_INOUT(volB, "_data");

	if (Adata == NULL || Bdata == NULL) {
		PyErr_SetString(PyExc_TypeError, "Volume have not complex numbers");
		return NULL;
	}

	assert(Amin != NULL && Ashape != NULL && Asize != NULL && Adelta != NULL);
	assert(Bmin != NULL && Bshape != NULL && Bsize != NULL && Bdelta != NULL);
	assert(Adata != NULL && Bdata != NULL);

#ifdef NDEBUG
#else
	if (debug) {
	printf("Amin  :%.2f %.2f %.2f:\n", (float)Amin[0], (float)Amin[1], (float)Amin[2]);
	printf("Adelta:%.2f %.2f %.2f:\n", (float)Adelta[0], (float)Adelta[1], (float)Adelta[2]);
	printf("Bmin  :%.2f %.2f %.2f:\n", (float)Bmin[0], (float)Bmin[1], (float)Bmin[2]);
	printf("Bdelta:%.2f %.2f %.2f:\n", (float)Bdelta[0], (float)Bdelta[1], (float)Bdelta[2]);
	}
#endif

	iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)Bdata);
	assert(iter->nd_m1 == 2);

	x[3] = 1.;
	PyArray_ITER_RESET(iter);
	while( PyArray_ITER_NOTDONE(iter) ) {
		v = (complex_t *) iter->dataptr;
		for (i = 0; i < 3; i++) x[i] = Bmin[i] + Bdelta[i]*iter->coordinates[i];
		dot(T, x, u);
		get_value(u, Amin, Adelta, Adata, v);
#ifdef NDEBUG
#else
		if(debug) {
		printf("c: %3i %3i %3i: ", (int)iter->coordinates[0], (int)iter->coordinates[1], (int)iter->coordinates[2]);
		printf("x: %-3.3f %-3.3f %-3.3f: ", (float)x[0], (float)x[1], (float)x[2]);
		printf("u: %-3.3f %-3.3f %-3.3f: ", (float)u[0], (float)u[1], (float)u[2]);
		printf("v: %-3.3f + %-3.3fi:\n", (float)v->real, (float)v->imag);
		}
#endif
		_PyAIT(iter)->index++;
		_PyArray_ITER_NEXT3(iter);
	}
	Py_DECREF(iter);
	Py_DECREF(T);
	Py_DECREF(Adata);
	Py_DECREF(Bdata);

	Py_RETURN_NONE;
}

static PyMethodDef vol_methods[] = {
    {"transfer",  Volume_transfer, METH_VARARGS,
     "Transfer a volume data to other with a geometry transformation."},
    {"put_ball",  Volume_put_ball, METH_VARARGS,
     "Put a ball in a volume."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


/*
 * Iniciación del modulo.
 */
PyMODINIT_FUNC
initvol_c(void)
{
	PyObject *m;

	m = Py_InitModule("vol_c", vol_methods);
	import_array();
	import_BasicVolume();
}

