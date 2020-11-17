# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
# Primitive Geometry Object: Sphere
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

from numpy import array, allclose, dot
from numpy.linalg import norm

class sphere:
    def __init__(self, o, r):
        self.o = array(o)
        self.r = r

    def dist(self, p):
        return norm(p - self.o) - self.r

    def __eq__(a, b):
        """
        """
        return allclose(a.o, b.o) and self.r==self.r 

    def __str__(self):
        return "<sphere o:%s r:%f>" % (self.o, self.r)

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:

