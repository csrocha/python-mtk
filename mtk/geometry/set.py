# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
# Set of Primitive Geometry Objects
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

from numpy import array, allclose, dot
from numpy.linalg import norm

class set:
    def __init__(self, I):
        self.bag = set(I)

    def intersection(self, gset):
        raise NotImplemented


    def __contains__(self, obj):
        raise NotImplemented


def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:

