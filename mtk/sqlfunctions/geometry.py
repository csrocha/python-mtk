# -*- coding: ISO-8859-1 -*-
# $Id: storage.py 64 2009-05-19 17:29:15Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
""" Plugins de funciones geometricas para la clase Storage.
"""

from mtk.geometry import neighbours
from numpy import array

def sqlf_6_dist3d(Ax,Ay,Az,Bx,By,Bz):
    """
    Calcule 3d distance for SQL queries

    >>> sqlf_6_dist3d(0,0,0,0,0,0)
    0.0
    """
    return neighbours.d(array(map(float,[Ax,Ay,Az])), array(map(float,[Bx,By,Bz])))

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

