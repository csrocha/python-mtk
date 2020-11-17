# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
# Primitive Geometry Arc Object
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

from math import radians, degrees

from mtk.geometry.line import segment
from mtk.geometry.triangle import triangle
from numpy import array, dot, arccos, pi, allclose, isnan, cross, sin, cos, min, argmin, max, argmax, zeros, sign
from numpy.linalg import norm

def _mirror_(o, a, b):
        """
        Calcula el punto espejo de b sobre el vector oa.

        # Ejemplo básico
        >>> o = array([0,0,0])
        >>> a = array([10,0,0])
        >>> b = array([0,10,0])
        >>> _mirror_(o,a,b)
        array([  0, -10,   0])

        # Ejemplo simple
        >>> o = array([0,0,0])
        >>> a = array([10,0,0])
        >>> b = array([7.0711,-7.0711,0])
        >>> _mirror_(o,a,b)
        array([ 7.0711,  7.0711,  0.    ])
 
        # Ejemplo símple trasladado
        >>> o = array([20,10,0])
        >>> a = array([30,10,0])
        >>> b = array([27.0711,2.9289,0])
        >>> _mirror_(o,a,b)
        array([ 27.0711,  17.0711,   0.    ])

        # Ejemplo símple trasladado y rotado
        >>> o = array([20,10,0])
        >>> a = array([28.6603,15,0])
        >>> b = array([29.6593,7.4118,0])
        >>> _mirror_(o,a,b)
        array([ 22.58824745,  19.65928729,   0.        ])

        """
        A = a-o
        B = b-o
        C = (dot(B,A)/dot(A,A))*A
        mB = 2*C-B
        return mB+o

class circle:
    def __init__(self, n, o, r):
        """
        Create a primitive graphic object Circle with center o and radius r, coplanar to the plane with normal n and o belong to it.
        """
        self.N = array(n)
        self.o = array(o)
        self.r = r

    def plane(self):
        """
        Return the plane coplanar to the circle
        """
        return plane(self.N, self.o)

    def area(self):
        return pi*self.r**2

    def tangent(self, p):
        """
        >>> C = circle([0,0,1],[0,0,0], 10)
        >>> allclose(C.tangent([10,0,0]), [0,1,0])
        True
        >>> allclose(C.tangent([0,10,0]), [-1,0,0])
        True
        >>> allclose(C.tangent([-10,0,0]), [0,-1,0])
        True
        >>> allclose(C.tangent([0,-10,0]), [1,0,0])
        True
        """
        p = array(p, float)
        v = (p - self.o)
        v /= norm(v)
        b = self.o + ((cross(v, self.N) - v) / 3)*self.r
        mb = _mirror_(self.o, p, b) 
        mbb = mb - b
        return mbb/norm(mbb)

    def curvature(self, os, rs, p = None):
        """
        Return the curvature of the circle related to the sphere with center os and radius rs.

        >>> C = circle([0,0,1],[0,0,0], 10)
        >>> C.curvature([10,0,0], 10)
        0.10000000000000001

        >>> C = circle([0,0,1],[0,0,0], 5)
        >>> C.curvature([-5,0,0], 5)
        0.20000000000000001

        >>> C = circle([0,0,1],[0,0,0], 10000)
        >>> C.curvature([10000,0,0], 10000)
        0.0001

        """
        oc = self.o
        rc = self.r
        N = self.N
        do = oc - os
        ndo = norm(do)
        kn = do / (ndo*rc)
        import pdb; pdb.set_trace()
        assert(allclose(norm(kn),1./rc))
        if p is None:
            p = cross(N, (1,0,0))
        t = self.tangent(p)
        assert(allclose(norm(N),1.))
        assert(allclose(norm(t),1.))
        skg = sign(dot(kn, cross(N, t)))
        kg = abs(ndo / (rc*rs)) * skg
        return kg

    def __belong__(self, p):
        return p in self.plane() and norm(p-self.o) == self.r

class arc:
    def __init__(self, o, r, a, b, N=None, rtol=1e-05, atol=1e-08):
        """
        Create a primitive graphic object Arc related to a circle with center o and radius r.
        The arc start in point a and end in point b.

        Restrictions:
        Radius of the circle must be positive and not zero.
        Both point must belong to the circle, if not raise assert error (AssertionError).
        """
        self.o = array(o, float)
        self.r = float(r)
        self.a = array(a, float)
        self.b = array(b, float)
        if N is None:
            self.N = cross(self.a - self.o, self.b - self.o)
            if any(isnan(self.N)) or allclose(self.N, zeros(len(self.N))):
                raise RuntimeError('You must specify the plane for this arc')
        else:
            self.N = array(N, float)
        self.N /= norm(self.N)
        if not (r > 0):
                raise RuntimeError('I expect a circle (r > 0), not a point. r = %5.16f' % r)
        if not (allclose(norm(self.a - self.o), r)):
                raise RuntimeError('I expect |a-o| == r, but |a-o| == %5.16f and r = %5.16f' % (norm(self.a - self.o), r))
        if not (allclose(norm(self.b - self.o), r)):
                raise RuntimeError('I expect |b-o| == r, but |b-o| == %5.16f and r = %5.16f' % (norm(self.b - self.o), r))

    def circle(self):
        """
        Return the circle associated to the arc.
        """
        return circle(self.N, self.o, self.r)

    def angle(self):
        """
        Return the angle of the arc over the circle.

        >>> A = arc([0,0,0], 10, [10,0,0], [0,10,0])
        >>> degrees(A.angle())
        90.0

        >>> A = arc([0,0,0], 10, [10,0,0], [10,0,0], N=[0,0,1])
        >>> degrees(A.angle())
        0.0

        >>> A = arc([0,0,0], 10, [10,0,0], [-10,0,0], N=[0,0,1])
        >>> degrees(A.angle())
        180.0
        """
        return arccos(dot((self.a - self.o) / self.r, (self.b - self.o) / self.r))

    def length(self):
        return self.angle()*self.r

    def curvature(self):
        """
        >>> A = arc([0,0,0], 10, [10,0,0], [0,10,0])
        >>> A.curvature()
        0.10000000000000001

        >>> A = arc([0,0,0], 10, [10,0,0], [10,0,0], N=[0,0,1])
        >>> A.curvature()
        0.10000000000000001

        >>> A = arc([0,0,0], 10, [10,0,0], [-10,0,0], N=[0,0,1])
        >>> A.curvature()
        0.10000000000000001

        >>> A = arc([-0.92182400498382733, -0.045187763882833146, 0.99639835812885114], \
                    0.01994653, [-0.91334669824142878, -0.042479022595480906, 1.0142494700075555], \
                    [-0.91458105639989895, -0.03514017512460646, 1.0120332505246663], \
                    N=[-0.90184917, -0.02137358,  0.431522  ])
        >>> A.curvature()
        50.134033338129491
        """
        return self.circle().curvature(self.o, self.r, p = self.a)

    def __str__(self):
        return 'arc(%s, %f, %s, %s, N=%s)' % (list(self.o), self.r, list(self.a), list(self.b), list(self.N))

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4

