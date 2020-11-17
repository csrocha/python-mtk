#ifndef __VOL_C_H__
#define __VOL_C_H__
#ifdef __cplusplus
extern "C" {
#endif

#include <Python.h>
#include <numpy/arrayobject.h>
#include <complex.h>


#define INT_T NPY_INTP
typedef npy_intp int_t;
#define REAL_T NPY_DOUBLE
typedef npy_double real_t;
#if 1
#define COMPLEX_T NPY_CDOUBLE
typedef npy_cdouble complex_t;
#else
#define COMPLEX_T NPY_CFLOAT
typedef npy_cfloat complex_t;
#endif

#define NPY_ARRAY_IN_TYPE(ARR, TYPE_ID)      (PyArrayObject *)PyArray_FromAny(ARR,PyArray_DescrFromType(TYPE_ID), 0, 0, NPY_IN_ARRAY, NULL)
#define NPY_ARRAY_INOUT_TYPE(ARR, TYPE_ID)   (PyArrayObject *)PyArray_FromAny(ARR,PyArray_DescrFromType(TYPE_ID), 0, 0, NPY_INOUT_ARRAY, NULL)

#define NPY_ARRAY_IN_REAL(ARR)               NPY_ARRAY_IN_TYPE(ARR, REAL_T)
#define NPY_ARRAY_INOUT_REAL(ARR)            NPY_ARRAY_INOUT_TYPE(ARR, REAL_T)

#define NPY_ARRAY_IN_COMPLEX(ARR)               NPY_ARRAY_IN_TYPE(ARR, COMPLEX_T)
#define NPY_ARRAY_INOUT_COMPLEX(ARR)            NPY_ARRAY_INOUT_TYPE(ARR, COMPLEX_T)

#define NPY_DATA_ARRAY_IN_TYPE(ARR, TYPE_ID, TYPE)    (TYPE *)PyArray_DATA(NPY_ARRAY_IN_TYPE(ARR, TYPE_ID))
#define NPY_DATA_ARRAY_INOUT_TYPE(ARR, TYPE_ID, TYPE) (TYPE *)PyArray_DATA(NPY_ARRAY_INOUT_TYPE(ARR, TYPE_ID))

#define NPY_DATA_ARRAY_IN_REAL(ARR)   NPY_DATA_ARRAY_IN_TYPE(ARR, REAL_T, real_t)
#define NPY_DATA_ARRAY_INOUT_REAL(ARR)   NPY_DATA_ARRAY_IN_TYPE(ARR, REAL_T, real_t)

#define NPY_DATA_ARRAY_IN_COMPLEX(ARR)   NPY_DATA_ARRAY_IN_TYPE(ARR, COMPLEX_T, complex_t)
#define NPY_DATA_ARRAY_INOUT_COMPLEX(ARR)   NPY_DATA_ARRAY_IN_TYPE(ARR, COMPLEX_T, complex_t)

#define GET_ATTR_ARRAY_IN(OBJ, ATTR, TYPE_ID)     NPY_ARRAY_IN_TYPE(PyObject_GetAttrString( OBJ, ATTR ), TYPE_ID)
#define GET_ATTR_ARRAY_INOUT(OBJ, ATTR, TYPE_ID)  NPY_ARRAY_INOUT_TYPE(PyObject_GetAttrString( OBJ, ATTR ), TYPE_ID)

#define GET_ARRAY_IN_REAL(OBJ, ATTR) GET_ATTR_ARRAY_IN(OBJ, ATTR, REAL_T)
#define GET_ARRAY_INOUT_REAL(OBJ, ATTR) GET_ATTR_ARRAY_INOUT(OBJ, ATTR, REAL_T)

#define GET_ARRAY_IN_COMPLEX(OBJ, ATTR) GET_ATTR_ARRAY_IN(OBJ, ATTR, COMPLEX_T)
#define GET_ARRAY_INOUT_COMPLEX(OBJ, ATTR) GET_ATTR_ARRAY_INOUT(OBJ, ATTR, COMPLEX_T)

#define GET_DATA_ATTR_ARRAY_IN(OBJ, ATTR, TYPE_ID, TYPE)    (TYPE *)PyArray_DATA(GET_ATTR_ARRAY_IN( OBJ, ATTR, TYPE_ID ))
#define GET_DATA_ATTR_ARRAY_INOUT(OBJ, ATTR, TYPE_ID, TYPE) (TYPE *)PyArray_DATA(GET_ATTR_ARRAY_INOUT( OBJ, ATTR, TYPE_ID ))

#define GET_ATTR_COMPLEX_IN(OBJ, ATTR)         GET_ATTR_ARRAY_IN(OBJ, ATTR, COMPLEX_T)
#define GET_ATTR_COMPLEX_INOUT(OBJ, ATTR)      GET_ATTR_ARRAY_INOUT(OBJ, ATTR, COMPLEX_T)
#define GET_DATA_ATTR_COMPLEX_IN(OBJ, ATTR)    GET_DATA_ATTR_ARRAY_IN(OBJ, ATTR, COMPLEX_T, complex_t)
#define GET_DATA_ATTR_COMPLEX_INOUT(OBJ, ATTR) GET_DATA_ATTR_ARRAY_INOUT(OBJ, ATTR, COMPLEX_T, complex_t)

#define GET_ATTR_REAL_IN(OBJ, ATTR)         GET_ATTR_ARRAY_IN(OBJ, ATTR, REAL_T)
#define GET_ATTR_REAL_OUT(OBJ, ATTR)        GET_ATTR_ARRAY_OUT(OBJ, ATTR, REAL_T)
#define GET_DATA_ATTR_REAL_IN(OBJ, ATTR)    GET_DATA_ATTR_ARRAY_IN(OBJ, ATTR, REAL_T, real_t)
#define GET_DATA_ATTR_REAL_INOUT(OBJ, ATTR) GET_DATA_ATTR_ARRAY_INOUT(OBJ, ATTR, REAL_T, real_t)

#define pmod(X,N) ((X<0)*N+(X%N))

/*
 * Calcula la multiplicación entre una matriz de transformación y un vector.
 *
 * dot( T, x, u )
 *
 * T: Matriz afin 4x4 de transformación
 * x: Punto a transformar
 * u: Punto resultante
 */
static inline
void dot(const PyArrayObject *T, const real_t x[4], real_t u[4]) {
	int i, j;
	assert(PyArray_DIM(T,0) == 4 && PyArray_DIM(T,1) == 4 );
	for (i = 0; i < 4; i++)
		for (j = 0, u[i] = 0; j < 4; j++)
			u[i] += *(real_t*)PyArray_GETPTR2(T, i, j)*x[j];
	for (i = 0; i < 4; i++) u[i] /= u[3];
}

/*
 * Resuelve el valor de un punto en el espacio a partir de los vecinos de la grilla
 * 
 * trilineallint(Xd, v)
 *
 * Xd: Coordenada en el espacio
 * v:  Valores de los vecinos más cercanos que existen en la grilla
 */
static inline
complex trilinealint(const float *Xd, const complex v[8]) {
	float  zd =   Xd[2];
	float lzd = 1-Xd[2];
	float  yd =   Xd[1];
	float lyd = 1-Xd[1];
	float  xd =   Xd[0];
	float lxd = 1-Xd[0];
	complex i1 = v[0]*lzd + v[1]*zd;
	complex i2 = v[2]*lzd + v[3]*zd;
	complex j1 = v[4]*lzd + v[5]*zd;
	complex j2 = v[6]*lzd + v[7]*zd;
	complex w1 =   i1*lyd +   i2*yd;
	complex w2 =   j1*lyd +   j2*yd;
	return w1*lxd+w2*xd;
}

/*
 * Obtiene el valor de cualquier coordenada a partir de una grilla.
 * El espacio es continuo y repetitivo.
 *
 * get_value(X, min, delta, array, v)
 *
 * X: coordenada a obtener el valor
 * min: valor mínimo de la grilla
 * delta: separación entre puntos de la grilla
 * array: arreglo de los datos de la grilla
 * v: valor resultante
 */
static inline
void get_value(const real_t *X, const real_t *min, const real_t *delta, const PyArrayObject *array, complex_t *v) {
	uint i ;
	float  Xo[3];
	float  Xd[3];
	float  Xt[3];
	complex_t vi[8];
	complex vo[8];
	npy_intp *dims = PyArray_DIMS(array);
	npy_intp o[2][3];
	int all = -1;
	for (i = 0; i < 3; i++) Xo[i] = (X[i] - min[i]) / delta[i];
	for (i = 0; i < 3; i++) Xt[i] = floor(Xo[i]);
	for (i = 0; i < 3; i++) Xd[i] = Xo[i] - Xt[i];
	for (i = 0; i < 3; i++) all &= Xd[i] == 0;
	for (i = 0; i < 3; i++) o[0][i] = pmod((int)Xt[i],dims[i]);
	if (all) { 
		*v = *(complex_t*)PyArray_GETPTR3(array, o[0][0], o[0][1], o[0][2]);
		return;
	}
	for (i = 0; i < 3; i++) o[1][i] = pmod((int)(Xt[i]+1),dims[i]);
	vi[0] = *(complex_t*)PyArray_GETPTR3(array, o[0][0], o[0][1], o[0][2]);
	vi[1] = *(complex_t*)PyArray_GETPTR3(array, o[0][0], o[0][1], o[1][2]);
	vi[2] = *(complex_t*)PyArray_GETPTR3(array, o[0][0], o[1][1], o[0][2]);
	vi[3] = *(complex_t*)PyArray_GETPTR3(array, o[0][0], o[1][1], o[1][2]);
	vi[4] = *(complex_t*)PyArray_GETPTR3(array, o[1][0], o[0][1], o[1][2]);
	vi[5] = *(complex_t*)PyArray_GETPTR3(array, o[1][0], o[0][1], o[0][2]);
	vi[6] = *(complex_t*)PyArray_GETPTR3(array, o[1][0], o[1][1], o[1][2]);
	vi[7] = *(complex_t*)PyArray_GETPTR3(array, o[1][0], o[1][1], o[0][2]);
	for(i = 0; i < 8; i++) vo[i] = vi[i].real + vi[i].imag * I;
	complex r = trilinealint(Xd, vo);
	v->real = creal(r);
	v->imag = cimag(r);
}

/*
 * Importa el modulo volumen
 */

static PyClassObject *BasicVolume_Class;

static inline
int import_BasicVolume(void)
{
	PyObject *declarations, *i;

	declarations = PyImport_ImportModule("mtk.geometry.vol");
	if (declarations == NULL)
		return -1;

	i = PyObject_GetAttrString(declarations, "BasicVolume");
	if (i == NULL)
		return -1;

	if (!PyClass_Check(i)) {
		PyErr_SetString(PyExc_TypeError, "mtk.geometry.vol.BasicVolume is not a type");
		return -1;
	}
	BasicVolume_Class = (PyClassObject *)i;

	Py_DECREF(declarations);

	Py_INCREF(BasicVolume_Class);

	return 0;

}

#ifdef __cplusplus
}
#endif
#endif /* __VOL_C_H__ */
