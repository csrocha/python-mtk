# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

from numpy import array, allclose, dot, zeros, cross, nanmax, seterr
from numpy.linalg import norm

class line:
    def __init__(self, d, o):
        self.d = array(d, dtype=float)
        self.o = array(o, dtype=float)

    def normalize(self):
        """
        Normalize the line

        >>> l = line([1,0,0], [0,0,0])
        >>> d = map(l.dist, [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ])
        >>> l.normalize()
        >>> allclose(d, map(l.dist, [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ]))
        True

        >>> l = line([1,0,0], [1,0,0])
        >>> d = map(l.dist, [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ])
        >>> l.normalize()
        >>> allclose(d, map(l.dist, [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ]))
        True

        > >> l = line([1,-1,0], [1,-1,1])
        > >> d = map(l.dist, [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ])
        > >> l.normalize()
        > >> (d, map(l.dist, [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ]))
        > >> allclose(d, map(l.dist, [ [0,0,0], [1,0,0], [0,1,0], [0,0,1] ]))
        True
        """
        n = norm(self.d)
        self.d /= n
        self.o /= n

    def normalized(self):
        """
        Create a new normalized line
        """
        n = norm(self.d)
        return line(self.d/n, self.o/n)

    def dist(self, p):
        """
        >>> l1 = segment([1,0,0], [0,0,0]).line()
        >>> l1.dist([1,0,0])
        0.0
        >>> l1.dist([0,0,0])
        0.0
        >>> l1.dist([1,1,0])
        1.0
        >>> l1.dist([-1,0,0])
        0.0
        >>> l2 = segment([1, 0, 0],[0, 1, 0]).line()
        >>> l2.dist([ 7.55336799, -6.55336799,  0.        ])
        0.0
        >>> l2.dist([-6.55336799,  7.55336799,  0.        ])
        0.0
        """
        x0 = array(p)
        x1 = self.o
        x2 = self.o + self.d
        return norm(cross(x0-x1, x0-x2))/norm(self.d)

    def __call__(self, v):
        return self.d * v + self.o

    def __eq__(self, l):
        """
        >>> l1 = segment([1,0,0], [0,0,0]).line()
        >>> l2 = segment([2,0,0], [1,0,0]).line()
        >>> l1 == l2
        True
        """
        return sum(map(self.dist, map(l, [-1, 0, 1]))) == 0.0

    def __str__(self):
        return "{ x %s + %s | x in R }" % (str(self.d),str(self.o))

class segment:
    def __init__(self, a, b):
        self.a = array(a, dtype=float)
        self.b = array(b, dtype=float)

    def affine(self, p):
        """
        >>> a = segment([  0.,-10. , 0.],[ 0., 0., 0.])
        >>> a.affine([0.,-5.,0.])
        0.5

        >>> a.affine([0.,5.,0.])
        -0.5

        """
        olderr = seterr(divide='ignore')
        nom = (self.a - self.b)
        if allclose(nom, zeros(nom.shape)): return False
        res = nanmax((p - self.b) / nom)
        seterr(**olderr)
        return res

    def line(self):
        """
        Return the line associated to the segment
        """
        return line(self.a-self.b, self.b)

    def __call__(self, i):
        """
        Return a point of the segment using affine coordinates
        """
        return self.a * i + self.b * (1-i)

    def __str__(self):
        return "[%s-%s]" % (str(self.a), str(self.b))

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:

