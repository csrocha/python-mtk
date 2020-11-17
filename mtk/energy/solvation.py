# -*- coding: ISO-8859-1 -*-
# $Id: solvation.py 92 2009-06-11 17:41:32Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Modulo que calcula la funcion de solvatacion
"""
from time import time
from numpy import array
from mtk.geometry.neighbours import d as dist3d

def _solvation_old_(S, radius=6.0, filter=[], pdbid=str()):
    """
    Calcule solvation free energy of a molecule

    >> > from mtk.tests import testfile
    >> > from mtk.storage import Storage
    >> > S = Storage()
    >> > id = S.loadpdb(testfile("obj01.pdb"))
    >> > sol = solvation(S)
    >> > abs((sol - -3429.536) / sol) < 0.00001
    True
    """

    sql = \
        "SELECT SUM(c.e) "\
        "FROM molecule AS ma, molecule AS mb, "\
        "     equivalents AS ea, equivalents AS eb, "\
        "     contacts AS c "\
        "WHERE ma.serial < mb.serial "\
        "  AND mb.pdbid==? "\
        "  AND dist3d(ma.x,ma.y,ma.z, mb.x,mb.y,mb.z) <= %f "\
        "  AND ma.name==ea.atom AND ma.resName == ea.aa "\
        "  AND mb.name==eb.atom AND mb.resName == eb.aa "\
        "  AND ((c.A == ea.class AND c.B == eb.class) OR "\
        "       (c.B == ea.class AND c.A == eb.class)) "\
        % radius
    filter = " AND ".join([sql,] + filter)
    q = S.do(sql, pdbid)
    return q.next()[0]

from time import time

class solvation_vol:
    """
    Generate a Volume with solvation information

    >>> from mtk.tests import testfile
    >>> from mtk.storage import Storage
    >>> from mtk.geometry.vol import create_from_coords
    >>> S = Storage()
    >>> id = S.loadpdb(testfile("obj01.pdb"))
    >>> C = S.get_coords()
    >>> V = create_from_coords(C, resolution=1.7)

    Interminable test. Impossible to solve.
    TODO: Accelerate it.
    > >> V.forallnodes(solvation_vol(S, 6.0))
    """

    def __init__(self, S, radius=6.0, pdbid=str()):
        self.S = S
        self.pdbid = pdbid
        self.sql = """
            SELECT SUM(e)
            FROM sol%(id)i
            WHERE
                dist3d(_ax,_ay,_az, ?,?,?) <= %(radius)f
            AND dist3d(_bx,_by,_bz, ?,?,?) <= %(radius)f
            """ % { 'id': id(self), 'radius': radius }
        sql = ["""
            CREATE TEMP TABLE sol%(id)i
            (_ax REAL, _ay REAL, _az REAL,
             _bx REAL, _by REAL, _bz REAL,
             e REAL)
             """ % { 'id': id(self) }, """
            INSERT INTO sol%(id)i
            SELECT ma.x, ma.y, ma.z,
                mb.x, mb.y, mb.z,
                c.e
            FROM molecule AS ma, molecule AS mb,
                 equivalents AS ea, equivalents AS eb,
                 contacts AS c
            WHERE
                dist3d(ma.x, ma.y, ma.z, mb.x, mb.y, mb.z) <= 2*%(radius)f
            AND mb.pdbid==ma.pdbid
            AND ma.serial < mb.serial
            AND ma.name==ea.atom AND ma.resName == ea.aa
            AND mb.name==eb.atom AND mb.resName == eb.aa
            AND ((c.A == ea.class AND c.B == eb.class) OR
                (c.B == ea.class AND c.A == eb.class))
               """ % { 'id': id(self), 'radius': radius } ]
        sql = ["""
            CREATE TEMP TABLE sol%(id)i
            (_ax REAL, _ay REAL, _az REAL,
             _bx REAL, _by REAL, _bz REAL,
             e REAL)
             """ % { 'id': id(self) }, """
            INSERT INTO sol%(id)i
            SELECT ma.x, ma.y, ma.z,
                mb.x, mb.y, mb.z,
                c.e
            FROM molecule AS ma, molecule AS mb,
                 equivalents AS ea, equivalents AS eb,
                 contacts AS c
            WHERE
                mb.pdbid==ma.pdbid
            AND ma.serial < mb.serial
            AND ma.name==ea.atom AND ma.resName == ea.aa
            AND mb.name==eb.atom AND mb.resName == eb.aa
            AND ((c.A == ea.class AND c.B == eb.class) OR
                (c.B == ea.class AND c.A == eb.class))
            """ % { 'id': id(self), 'radius': radius } ]
        self.S._profile = False
        self.S.do(sql)
        self.count = 0
        self.radius = radius

    def __del__(self):
        self.S.do("DROP TABLE sol%i" % id(self))
        self.S.do("VACUUM")

    def __call__(self, x, v):
        values = list(x[:3]) * 2
        t = time()
        r = self.S.do(self.sql, *values).next()[0]
        #print "X:", time() - t, self.count
        self.count  = self.count + 1
        if r == None: return 0.
        return r

    def get_pairs(self):
        r = self.S.do("""
            SELECT *
            FROM sol%i
                 """ % id(self) )
        return [ x for x in r if dist3d(array(x[:3]),array(x[3:-1])) <=
                self.radius ]

def gettypes(pdb):
    pass

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

