# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
# Primitive Geometry Object: Triangle
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

from numpy import array, allclose, dot, ndarray, sqrt, cross
from numpy.linalg import norm
from line import line

class plane:
    def __init__(self, n, o):
        self.n = array(n, float)
        self.o = o

    def normalize(self):
        self.n /= norm(self.n)

    def normalized(self):
        return plane(self.n / norm(self.n), o)

    def d(self, t=None):
        """
        >>> p = plane([2,1,-1],[0,0,1])
        >>> p.d()
        1.0

        >>> p = plane([0,2,0],[0,3./2.,0])
        >>> p.d()
        -3.0

        >>> p = plane([2,3,10],[0,0,-4])
        >>> p.d()
        40.0

        >>> p = plane([2,3,10],[0,0,-4])
        >>> allclose(p.d([0,0,-4]), 0)
        True

        """
        if t is None:
            return -dot(self.n,self.o)
        else:
            return -dot(self.n,self.o-array(t,float))

    def dist(self, p):
        """
        >>> t1 = triangle([0,0,0], [1,0,0], [0,1,0]).plane()
        >>> t1.dist([1,0,0])
        0.0

        >>> t1.dist([0,0,0])
        0.0

        >>> t1.dist([0,1,0])
        0.0

        >>> t1.dist([0.25,0.25,0])
        0.0

        >>> t1.dist([0.25,0.25,0.25])
        0.25

        >>> p = plane([2,1,-1],[0,0,1])
        >>> p.dist([3,1,-2])
        4.0824829046386304

        >>> p = plane([0,2,0],[0,3./2.,0])
        >>> p.dist([3,1,-2])
        0.5

        >>> p = plane([-11, -11, 121],[0,1,2])
        >>> p.dist([-3.87811061, -2.87811061,  1.29488898])
        0.0
        """
        d = self.d()
        u = abs(dot(self.n,p) + d)
        d = sqrt(dot(self.n,self.n))
        return u/d

    def project(self, p):
        """
        Projection of p over the surface

        >>> t1 = triangle([0,0,0], [1,0,0], [0,1,0]).plane()
        >>> allclose(t1.dist(t1.project([1,0,0])), 0)
        True

        >>> t1 = triangle([0,0,0], [1,0,0], [0,1,0]).plane()
        >>> allclose(t1.dist(t1.project([1,0,10])), 0)
        True
        """
        l = line(self.n, p)
        d = dot((self.o - l.o), self.n) / dot(l.d, self.n)
        return l(d) 

    def affine(self, a):
        o + dot(a, b)

    def __str__(self):
        return "<plane n:%s o:%s>" % (self.n, self.o)

def _same_side_(p1, p2, a, b):
        cp1 = cross(b-a, p1-a)
        cp2 = cross(b-a, p2-a)
        return dot(cp1, cp2) >= 0

class triangle:
    def __init__(self, a, b, c):
        self.p = (array(a, float), array(b, float), array(c, float))

    def normal(self):
        """
        Normal of the plane descripted by the face

        >>> t1 = triangle([0,0,0], [1,0,0], [0,1,0])
        >>> allclose(t1.normal(), array([ 0., 0.,  1.]))
        True

        >>> t1 = triangle([0,0,0], [0,0,1], [1,0,0])
        >>> allclose(t1.normal(), array([ 0., 1.,  0.]))
        True

        >>> t1 = triangle([0,0,0], [1,0,0], [0,0,1])
        >>> allclose(t1.normal(), array([ 0., -1.,  0.]))
        True

        >>> t1 = triangle([0,0,0], [0,1,0], [0,0,1])
        >>> allclose(t1.normal(), array([ 1., 0.,  0.]))
        True

        >>> t1 = triangle([1,0,0], [0,1,0], [0,0,1])
        >>> allclose(t1.normal(), array([ 0.57735,  0.57735,  0.57735]))
        True
        """
        N = cross(self.p[1] - self.p[0], self.p[2] - self.p[0])
        return N / norm(N)

    def plane(self):
        """
        Return the plane associated to the triangle

        >>> t = triangle([1,0,0], [0,1,0], [0,0,0])
        >>> p = t.plane()
        >>> allclose(p.n, array([ 0., 0.,  1.]))
        True

        >>> t = triangle([0,0,0], [0,1,0], [0,0,1])
        >>> p = t.plane()
        >>> allclose(p.n, array([ 1., 0., 0.]))
        True

        >>> t = triangle([1,0,0], [0,0,0], [0,0,1])
        >>> p = t.plane()
        >>> allclose(p.n, array([ 0., 1., 0.]))
        True
        """
        return plane(self.normal(), self.p[0])

    def lines(self):
        """
        >>> t1 = triangle([0,0,0], [1,0,0], [0,0,1])
        >>> map(str,list(t1.lines()))
        ['{ x [ 0.  0.  0.] + [ 0.  0.  0.] | x in R }', '{ x [ 1.  0.  0.] + [ 1.  0.  0.] | x in R }', '{ x [ 0.  0.  1.] + [ 0.  0.  1.] | x in R }']
        """
        for i in range(3):
            yield line(self.p[i], self.p[i%3])
        return

    def __call__(self, v):
        """
        Point over the triangle with affine coordinates

        >>> t1 = triangle([0,0,0], [1,0,0], [0,0,1])
        >>> allclose(t1([1,0,0]), array([0, 0, 0]))
        True
        >>> allclose(t1([0,1,0]), array([1, 0, 0]))
        True
        >>> allclose(t1([0,0,1]), array([0, 0, 1]))
        True
        >>> allclose(t1([0.5,0.5,0.5]), array([0.5, 0, 0.5]))
        True
        """
        return dot(array(self.p).T,array(v))

    def __eq__(a, b):
        """
        >>> t1 = triangle([0,0,0], [1,0,0], [0,1,0])
        >>> t1 == triangle([1,0,0], [0,0,0], [0,1,0])
        False

        >>> t1 == triangle([1,0,0], [0,1,0], [0,0,0])
        True
        """
        return allclose(a.p, b.p) or allclose(a.p[1:] + a.p[:1], b.p) or allclose(a.p[2:] + a.p[:2], b.p)

    def __contains__(self, p):
        return (_same_side_(p, self.p[0], self.p[1], self.p[2]) and 
                _same_side_(p, self.p[1], self.p[0], self.p[2]) and 
                _same_side_(p, self.p[2], self.p[0], self.p[1])) 

    def __str__(self):
        return "<triangle a:%s b:%s c:%s>" % self.p

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:

