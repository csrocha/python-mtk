# -*- coding: ISO-8859-1 -*-
# $Id: rms.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Calcule rms related functions"""

from numpy import array, shape, dot, transpose, identity
from numpy import sqrt, sum, max
from numpy.linalg import svd, det

def array_rmsd(arr1, arr2):
    return sqrt(sum((arr1 - arr2)**2)/len(arr1))

def rmsd(crds1, crds2):
    """Returns RMSD between 2 sets of [nx{3,4}] numpy array
    url: http://boscoh.com/protein/rmsd-root-mean-square-deviation

    >>> from numpy import array
    >>> a = array([ [ 0, 0, 0, 1 ] ])
    >>> b = array([ [ 1, 0, 0, 1 ] ])
    >>> rmsd(a, b)
    1.0

    >>> from numpy import array
    >>> a = array([ [ 0, 0, 0 ] ])
    >>> b = array([ [ 1, 0, 0 ] ])
    >>> rmsd(a, b)
    1.0

    >>> a = array([ [ 0, 0, 0, 1 ], [ 1, 0, 0, 1 ] ])
    >>> b = array([ [ 1, 0, 0, 1 ], [ 2, 0, 0, 1 ] ])
    >>> rmsd(a, b)
    1.0

    >>> a = array([ [ 0, 0, 0 ], [ 2, 0, 0 ] ])
    >>> b = array([ [ 1, 0, 0 ], [ 2, 0, 0 ] ])
    >>> "%08.4f" % rmsd(a, b)
    '000.7071'

    >>> a = array([ [ i, j, k ] for i in xrange(10) for j in range(10) for k in range(10) ])
    >>> b = array([ [ k, j, i ] for i in xrange(10) for j in range(10) for k in range(10) ])
    >>> "%08.4f" % rmsd(a, b)
    '005.7446'
    """
    assert(crds1.shape[1] > 2)
    assert(crds1.shape == crds2.shape)

    # Checking for affine coordinates
    if crds1.shape[1] > 3:
        crds1 = crds1[:,:3]
        crds2 = crds2[:,:3]

    E0 = sum(crds1 * crds1) + \
         sum(crds2 * crds2)
    S  = sum(crds1 * crds2)
    rms2 = (E0 - 2*S) / float(crds1.shape[0])
    assert(rms2 >= 0.0)

    return sqrt(rms2)

def fit_rotation(crds1, crds2, affine=False):
    """Returns best-fit_rotation rotation matrix as [3x{3,4}] numpy matrix
    url: http://boscoh.com/protein/rmsd-root-mean-square-deviation

    >>> from numpy import allclose
    >>> a = array([ [ 1, 0, 0, 1 ], [ 0, 1, 0, 1 ] ])
    >>> b = array([ [ 1, 0, 0, 1 ], [ 0, 1, 0, 1 ] ])
    >>> fit_rotation(a, b)
    array([[ 1.,  0.,  0.,  0.],
           [ 0.,  1.,  0.,  0.],
           [ 0.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  1.]])

    >>> a = array([ [ 0, 0, 0, 1 ], [ 1, 0, 0, 1 ] ])
    >>> b = array([ [ 0, 0, 0, 1 ], [ 0, 1, 0, 1 ] ])
    >>> R = fit_rotation(a, b)
    >>> rmsd( dot(a, R), b )
    0.0

    >>> a = array([ [ 0, 0, 0, 1 ], [ 0, 1, 0, 1 ], \
                    [ 1, 1, 0, 1 ], [ 1, 0, 0, 1 ] ])
    >>> b = array([ [ 0, 0, 0, 1 ], [ 0,-1, 0, 1 ], \
                    [-1,-1, 0, 1 ], [-1, 0, 0, 1 ] ])
    >>> R = fit_rotation(a, b)
    >>> allclose(rmsd( dot(a, R), b ), .0, atol=1e-7)
    True

    >>> a = array([ [ 0, 0, 0, 1 ], [ 0, 1, 0, 1 ], \
                    [ 1, 1, 0, 1 ], [ 1, 0, 0, 1 ] ])
    >>> b = array([ [ 0, 1, 0, 1 ], [ 1, 1, 0, 1 ], \
                    [ 1, 0, 0, 1 ], [ 0, 0, 0, 1 ] ])
    >>> ca, cb = centre(a), centre(b)
    >>> R = fit_rotation(a - ca, b -cb)
    >>> rmsd( dot(a-ca, R) + ca, b )
    0.0
    """
    assert(crds1.shape[1] > 2)
    assert(crds1.shape == crds2.shape)

    # Checking for affine coordinates
    if crds1.shape[1] > 3:
        crds1 = crds1[:,:3]
        crds2 = crds2[:,:3]
        affine = True

    correlation_matrix = dot(transpose(crds1), crds2)
    v, s, w = svd(correlation_matrix)
    is_reflection = (det(v) * det(w)) < 0.0
    if is_reflection:
        v[-1,:] = -v[-1,:]

    r = dot(v, w)

    if affine:
        I = identity(4)
        I[:3,:3] = r
        return I
    else:
        return r

def fit_transform(crds1, crds2):
    """Return the best fit transformation from crds1 to crds2

    >>> a = array([ [ 0, 0, 0, 1 ], [ 0, 1, 0, 1 ], \
                    [ 1, 1, 0, 1 ], [ 1, 0, 0, 1 ] ])
    >>> b = array([ [ 0, 1, 0, 1 ], [ 1, 1, 0, 1 ], \
                    [ 1, 0, 0, 1 ], [ 0, 0, 0, 1 ] ])
    >>> R = fit_transform(a, b)
    >>> rmsd( dot(a, R), b )
    0.0
    """
    ca, cb = centre(crds1), centre(crds2)
    T = identity(4)
    T[3,:3] = -ca[:3]
    T = dot(T, fit_rotation(crds1-ca, crds2-cb, affine = True))
    C = identity(4)
    C[3,:3] = ca[:3]
    T = dot(T, C)
    return T

def fit(crds1, crds2):
    """Return the best fit rmsd value of crds1 to crds2

    >>> a = array([ [ 0, 0, 0, 1 ], [ 0, 1, 0, 1 ], \
                    [ 1, 1, 0, 1 ], [ 1, 0, 0, 1 ] ])
    >>> b = array([ [ 0, 1, 0, 1 ], [ 1, 1, 0, 1 ], \
                    [ 1, 0, 0, 1 ], [ 0, 0, 0, 1 ] ])
    >>> fit(a, b)
    0.0

    >>> a = array([ [ i, j, k, 1 ] for i in xrange(10) for j in range(10) for k in range(10) ])
    >>> b = array([ [ k, j, i, 1 ] for i in xrange(10) for j in range(10) for k in range(10) ])
    >>> "%08.4f" % fit(a, b)
    '005.7446'
    """
    R = fit_transform(crds1, crds2)
    return rmsd( dot(crds1, R), crds2 )

def centre(crds):
    """Returns the geometric centre of crds

    >>> a = array([ [ 0, 0, 0, 1 ], [ 0, 1, 0, 1 ], \
                    [ 1, 1, 0, 1 ], [ 1, 0, 0, 1 ] ])
    >>> centre(a)
    array([ 0.5,  0.5,  0. ,  1. ])
    """
    n_vect = float(crds.shape[0])
    return array([ sum(crds[:,0])/n_vect, sum(crds[:,1])/n_vect,
                  sum(crds[:,2])/n_vect , 1.0])[:crds.shape[1]]

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

