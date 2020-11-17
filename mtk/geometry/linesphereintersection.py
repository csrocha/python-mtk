# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
#
from numpy import *
from numpy.linalg import norm
from line import line
from sphere import sphere

def line_sphere_intersect(l, s) :
    """
    Solve two points of intersection of a line and circle

    >>> s = sphere([0, 0, 0], 10.0)
    >>> a = line([0, 0, 0], [1, 0, 0])
    >>> f, g = line_sphere_intersect(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    >>> allclose([s.dist(f), s.dist(g)], [0.0, 0.0])
    True

    >>> a = line([1, 0, 0], [0, 0, 0])
    >>> f, g = line_sphere_intersect(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    >>> allclose([s.dist(f), s.dist(g)], [0.0, 0.0])
    True

    >>> a = line([1, 0, 0],[0, 1, 0])
    >>> f,g = line_sphere_intersect(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    >>> allclose([s.dist(f), s.dist(g)], [0.0, 0.0])
    True

    >>> a = line([1, 1, 0],[0, 1, 0])
    >>> f,g = line_sphere_intersect(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    >>> allclose([s.dist(f), s.dist(g)], [0.0, 0.0])
    True
    """
    o = l.o
    l = l.d
    c = s.o
    r = s.r
    LO = dot(l,o)
    LC = dot(l,c)
    OC = dot(o,c)
    A = LO - LC
    AA = A*A

    LL = dot(l,l)
    OO = dot(o,o)
    CC = dot(c,c)
    RR = r*r

    B = OO + CC - RR - 2*OC

    C = LL

    tsqr = AA - C*B

    if tsqr < 0:
        return False

    tsqr = sqrt(tsqr)
    k1 = (-A + tsqr)/LL
    k2 = (-A - tsqr)/LL

    return (l*k1+o, l*k2+o)

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:

