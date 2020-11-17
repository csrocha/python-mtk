# -*- coding: ISO-8859-1 -*-
# $Id: lj.py 76 2009-05-30 17:39:56Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Modulo que calcula la funcion de lennard-jones

>>> from pkg_resources import resource_filename
>>> from os.path import join as joinpath
>>> path = resource_filename('mtk', '')
>>> from mtk.storage import Storage
>>> S = Storage()
>>> S.loadpdb(joinpath(path, "data", "test", "obj01.pdb"))
'obj01.pdb'
>>> L = list(S.do("SELECT name FROM molecule GROUP BY name"))
>>> L == list(S.do("SELECT name FROM molecule INTERSECT SELECT symbol FROM potentials"))
True
>>> d = list(S.do("SELECT x,y,z,A,B FROM molecule, potentials WHERE symbol==name"))
>>> len(d)
597
>>> from mtk.geometry import neighbours as N
>>> from numpy import array
>>> x = array(d)
>>> D = N.distance(x)
>>> n = N.nearest(D, 3.0)
>>> p = N.pairs(n)
>>> "%8.4f" % lj(D, p, x[:,3:])
'316.6550'
"""

def lj(dist, net, parm):
    S = 0
    for i, j in net:
        r6 = dist[i][j]**6
        r12 = r6**2
        A = parm[i][0]
        B = parm[j][1]
        S += A/r6 + B/r12
    return S

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

