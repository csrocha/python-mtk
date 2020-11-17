# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Inicializador del modulo que calcula funciones geometricas

>>> from numpy import array, identity

Verificando la inicialización de las implementaciones en C y en Python

>>> import vol
>>> import vol_c
>>> p = list(vol.GridIterator((0,0,0),(1,2,3)))

Checking shape

>>> p = vol.Volume((0,0,0), (10,10,10), 1)
>>> (p.shape == array([11, 11, 11])).all()
True

Checking IndexError message

>>> p[(11,11,11)]
Traceback (most recent call last):
    ...
IndexError: index (11) out of range (0<=index<11) in dimension 0


Checking size without integer range

>>> p = vol.Volume((0,0,0), (9.7,9.3,9.0), 1)
>>> (p.shape == array([11, 11, 10])).all()
True

Checking transfer code.

>>> from numpy import identity, zeros, complex, cfloat
>>> from time import time
>>> transfer = vol_c.transfer
>>> p = vol.BasicVolume((0,0,0), (1,1,1), zeros((10,10,10), dtype=cfloat))
>>> q = vol.BasicVolume((0,0,0), (1,1,1), zeros((10,10,10), dtype=cfloat))
>>> T = identity(4)
>>> tm = time()
>>> transfer(p, q, T)
>>> time()-tm < 1
True

Checking transfer values without rotation

>>> from numpy import array, zeros, identity, all, allclose, reshape, arange
>>> from math import radians
>>> from mtk.geometry.transformation import rotation
>>> from mtk.geometry.vol import BasicVolume
>>> T = identity(4)
>>> A = reshape(arange(0,3**3),(3,3,3))
>>> p = vol.BasicVolume((-1,-1,-1), (1,1,1), A)
>>> q = vol.BasicVolume((-1,-1,-1), (1,1,1), zeros((3,3,3)))
>>> transfer(p, q, T)
>>> all(p._data == A)
True
>>> allclose(q._data, A)
True

Checking transfer values with rotation

>>> from numpy import array, zeros, identity, all, allclose, reshape, arange
>>> from math import radians
>>> from mtk.geometry.transformation import rotation
>>> from mtk.geometry.vol import BasicVolume
>>> T = rotation(radians(90),0,0)
>>> A = reshape(arange(0,3**3),(3,3,3))
>>> p = vol.BasicVolume((-1,-1,-1), (1,1,1), A)
>>> q = vol.BasicVolume((-1,-1,-1), (1,1,1), zeros((3,3,3)))
>>> transfer(p, q, T)
>>> all(p._data == A)
True
>>> R = [ 2, 5, 8, \
          1, 4, 7, \
          0, 3, 6, \
         \
         11, 14, 17, \
         10, 13, 16, \
          9, 12, 15, \
         \
         20, 23, 26, \
         19, 22, 25, \
         18, 21, 24]
>>> allclose(q._data.flat, R)
True

Checking transfer values with rescalation to increment size

>>> from numpy import array, zeros, identity, all, reshape, arange, savetxt
>>> from math import radians
>>> from mtk.geometry.transformation import rotation
>>> from mtk.geometry.vol import BasicVolume
>>> T = rotation(0,0,0)
>>> A = reshape(arange(0,20**3),(20,20,20))
>>> p = vol.BasicVolume((-10,-10,-10), (1,1,1), A)
>>> q = vol.BasicVolume((-10,-10,-10), (.5,.5,.5), zeros((40,40,40)))
>>> tm = time()
>>> transfer(p, q, T)
>>> time()-tm < 3
True
>>> all(p._data == A)
True
>>> allclose(q._data[0,0,:-1], arange(0,20,.5)[:-1])
True
>>> allclose(q._data[0,1,:-1], arange(10,30,.5)[:-1])
True
>>> allclose(q._data[0,2,:-1], arange(20,40,.5)[:-1])
True

Checking transfer values with rescalation to reduce size

>>> from numpy import array, zeros, identity, all, reshape, arange, savetxt
>>> from math import radians
>>> from mtk.geometry.transformation import rotation
>>> from mtk.geometry.vol import BasicVolume
>>> T = rotation(0,0,0)
>>> A = reshape(arange(0,20**3),(20,20,20))
>>> p = vol.BasicVolume((-10,-10,-10), (1,1,1), A)
>>> q = vol.BasicVolume((-10,-10,-10), (2,2,2), zeros((10,10,10)))
>>> tm = time()
>>> transfer(p, q, T)
>>> time()-tm < 3
True
>>> all(p._data == A)
True
>>> allclose(q._data[0,0,:-1], arange(0,20,2)[:-1])
True
>>> allclose(q._data[0,1,:-1], arange(40,60,2)[:-1])
True

Transfer with rotations

>>> import csv
>>> from time import time
>>> from mtk.geometry.transformation import genRotationsByDirections
>>> from numpy import hstack, vstack
>>> R = csv.reader(open("mtk/data/csv/UnitSphere_500.csv"))
>>> directions = array([ map(float, line) for line in R ])
>>> R = genRotationsByDirections(directions, radians(15))
>>> Osize = 200
>>> Rsize = Osize/2
>>> A = reshape(arange(0,Osize**3,dtype=complex),(Osize,Osize,Osize))
>>> p = vol.BasicVolume((-Rsize,-Rsize,-Rsize), (1,1,1), A)
>>> q = vol.BasicVolume((-Rsize,-Rsize,-Rsize), (2,2,2), zeros((Rsize,Rsize,Rsize), dtype=complex))
>>> t = time()
>>> len(R)
6048
>>> for r in R[:1]: q.transfer(p, r)
>>> time() - t < 2
True

  Tiempo testeado para todas las rotaciones con números reales
     4034.2435631752014

  Tiempo esperado para todas las rotaciones con números complejos
     9374.3999999999996

Time results:

6048 rotations in a Thurion 64 run 4034.2435631752014seg ~ 1hs 7min
1 rotation in a AMD Turion(tm) 64 X2 1.2GHz run in 0.60seg
1 rotation in a Intel(R) Xeon(TM) CPU 2.66GHz run in 0.77seg
1 rotation in a Intel(R) Atom(TM) CPU 330 @ 1.60GHz run in 1.49seg

Testing put_ball implemented in c

>>> from numpy import array, sum
>>> from mtk.geometry.vol_c import put_ball
>>> from mtk.geometry.vol import create_from_coords
>>> a = array([ [ 0, 0, 0 ] ])
>>> V = create_from_coords(a, resolution=1, inc=4.0)
>>> put_ball(V, a[0], {1:4, 2:5, 3:6}, [0,4,5,6])
>>> [ sum(V._data == i) for i  in [4,5,6] ]
[7, 26, 90]

"""
from numpy import array

# Coordinate operations
def origin(coords):
    """
    Return the origin coordinate of the minimal cube witch include all coordinates.

    >>> A = array([ [1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1] ])
    >>> origin(A)
    array([-1, -1, -1])
    """
    return array((coords[:,0].min(),coords[:,1].min(),coords[:,2].min()))

def extreme(coords):
    """
    Return the most extreme coordinate of the minimal cube witch include all coordinates.

    >>> A = array([ [1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1] ])
    >>> extreme(A)
    array([1, 1, 1])
    """
    return array((coords[:,0].max(),coords[:,1].max(),coords[:,2].max()))

def size(coords):
    """
    Return the size of the minimal cube witch include all coordinates.

    >>> A = array([ [1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1] ])
    >>> size(A)
    array([2, 2, 2])
    """
    return extreme(coords) - origin(coords)

def centre(coords):
    """
    Return the centre of the minimal cube witch include all coordinates.

    >>> A = array([ [1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1] ])
    >>> centre(A)
    array([0, 0, 0])
    """
    return (origin(coords) + extreme(coords)) / 2

# import vol as volume

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

