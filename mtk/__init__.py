# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 76 2009-05-30 17:39:56Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""This module is the Molecular Toolkit.
The main objetive of the module is molecular research.

.. moduleauthor:: Cristian S. Rocha <crocha@dc.uba.ar>

"""

import forcefield
import inferator

"""Inicializador del modulo mtk

    *** Testing reading files and geometric calculus ***

>>> from pkg_resources import resource_filename
>>> from os.path import join as joinpath
>>> import io.pdb_ff as pdbff
>>> import geometry.rms as rms
>>> path = resource_filename('mtk', '')
>>> obj01 = pdbff.reader(open(joinpath(path, "data", "test", "obj01.pdb")))
>>> obj02 = pdbff.reader(open(joinpath(path, "data", "test", "obj02.pdb")))
>>> crd01 = obj01.get_coords()
>>> crd02 = obj02.get_coords()
>>> "%08.4f" % rms.rmsd(crd01, crd02)
'001.0794'
>>> "%08.4f" % rms.fit(crd01, crd02)
'001.0785'

>>> obj01 = pdbff.reader(open(joinpath(path, "data", "test", "obj01ca.pdb")))
>>> obj02 = pdbff.reader(open(joinpath(path, "data", "test", "obj02ca.pdb")))
>>> crd01 = obj01.get_coords()
>>> crd02 = obj02.get_coords()
>>> "%08.4f" % rms.rmsd(crd01, crd02)
'000.4633'
>>> "%08.4f" % rms.fit(crd01, crd02)
'000.4633'

   NOTE: Pymol resolve better fit than our fit function.
   RMS is equal than pymol.

    *** Prepare data ***

>>> from storage import Storage
>>> from numpy import array
>>> S = Storage()
>>> S.loadpdb(joinpath(path, "data", "test", "obj01.pdb"))
'obj01.pdb'

    *** Generate relacion between the molecule and potentials ***

>>> q = array(list(S.do("SELECT x,y,z,A,B FROM molecule AS m, potentials AS p "\
               "WHERE m.name==p.symbol AND m.name='CA'")))

    *** Disjoin Coordinates of Potentials ***

>>> C = q[:,:3]
>>> P = q[:,3:]

    *** Calculate distances ***

>>> from geometry import neighbours
>>> D = neighbours.distance(C)

    *** Calculate LJ of neighbours of x distance ***

>>> from energy import lj
>>> N = neighbours.pairs(neighbours.nearest(D, 1.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.0000', 0)
>>> N = neighbours.pairs(neighbours.nearest(D, 3.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.0000', 0)
>>> N = neighbours.pairs(neighbours.nearest(D, 4.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.1022', 164)
>>> N = neighbours.pairs(neighbours.nearest(D, 7.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.1331', 662)
>>> N = neighbours.pairs(neighbours.nearest(D, 10.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.1365', 1376)
>>> N = neighbours.pairs(neighbours.nearest(D, 15.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.1379', 3480)
>>> N = neighbours.pairs(neighbours.nearest(D, 16.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.1379', 3918)
>>> N = neighbours.pairs(neighbours.nearest(D, 17.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.1380', 4354)
>>> N = neighbours.pairs(neighbours.nearest(D, 20.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.1380', 5522)
>>> N = neighbours.pairs(neighbours.nearest(D, 40.0))
>>> "%08.4f" % lj.lj(D, N, P), len(N)
('000.1380', 6806)

    *** Select pairs of atoms with distance less of 7.0 A ***

>>> q = S.do(\
        "SELECT ma.serial, mb.serial  "\
        "FROM molecule AS ma, molecule AS mb "\
        "WHERE ma.serial < mb.serial AND ma.name == 'CA' AND mb.name=='CA'"\
        "  AND dist3d(ma.x,ma.y,ma.z, mb.x,mb.y,mb.z) <= 7.0 "\
        )
>>> L = list(q)
>>> len(L)
331

    *** Select pairs of atoms with distance less of 7.0 A, with equivalents
        table ***

>>> q = S.do(\
        "SELECT ma.serial, mb.serial "\
        "FROM molecule AS ma, molecule AS mb, "\
        "     equivalents AS ea, equivalents AS eb "\
        "WHERE ma.serial < mb.serial AND ma.name == 'CA' AND mb.name=='CA'"\
        "  AND dist3d(ma.x,ma.y,ma.z, mb.x,mb.y,mb.z) <= 7.0 "\
        "  AND ma.name==ea.atom AND ma.resName == ea.aa "\
        "  AND mb.name==eb.atom AND mb.resName == eb.aa "\
        )
>>> M = list(q)
>>> len(M)
331

    *** Select pairs of atoms with distance less of 7.0 A, with equivalents
        table and contact information ***

>>> q = S.do(\
        "SELECT ma.serial, mb.serial, c.e, c.ep "\
        "FROM molecule AS ma, molecule AS mb, "\
        "     equivalents AS ea, equivalents AS eb, "\
        "     contacts AS c "\
        "WHERE ma.serial < mb.serial AND ma.name=='CA' AND mb.name=='CA'"\
        "  AND dist3d(ma.x,ma.y,ma.z, mb.x,mb.y,mb.z) <= 7.0 "\
        "  AND ma.name==ea.atom AND ma.resName == ea.aa "\
        "  AND mb.name==eb.atom AND mb.resName == eb.aa "\
        "  AND ((c.A == ea.class AND c.B == eb.class) OR "\
        "       (c.B == ea.class AND c.A == eb.class)) "\
        )
>>> M = list(q)
>>> len(M)
331

    *** Calculate sum of contacts of the same molecule ***

>>> q = S.do(\
        "SELECT SUM(c.e) "\
        "FROM molecule AS ma, molecule AS mb, "\
        "     equivalents AS ea, equivalents AS eb, "\
        "     contacts AS c "\
        "WHERE ma.serial < mb.serial AND ma.name=='CA' AND mb.name=='CA'"\
        "  AND dist3d(ma.x,ma.y,ma.z, mb.x,mb.y,mb.z) <= 6.0 "\
        "  AND ma.name==ea.atom AND ma.resName == ea.aa "\
        "  AND mb.name==eb.atom AND mb.resName == eb.aa "\
        "  AND ((c.A == ea.class AND c.B == eb.class) OR "\
        "       (c.B == ea.class AND c.A == eb.class)) "\
        )
>>> q.next() == (-175.28100000000052,)
True

"""

def have_gpu():
    """Check if pycuda is installed

    :returns: bool -- Is true if exists pycuda module.
    """
    try:
        import pycuda.autoinit
        return True
    except:
        return False

def have_mayavi():
    """Check if mayavi is installed.

    :returns: bool -- Is true if exists mayavi module.
    """
    try:
        import enthought.mayavi 
        return True
    except:
        return False

def have_vtk():
    """Check if vtk is installed.

    :returns: bool -- Is true if exists vtk module.
    """
    try:
        import enthought.tvtk 
        return True
    except:
        return False

def resource(library, filename):
    """Return the full path to the file in an standard library.

    :param library: repository name where the file is stored.
    :type library: str.
    :param filename: name of the file to get full path.
    :type filename: str.
    :returns: str -- The full path to the file.

    In the following example whe find the 'UnitSphere_500.csv'
    file in the 'csv' library and we check if it exists.
    Here the file exists.

    >>> from os.path import exists
    >>> fn = resource("csv", "UnitSphere_500.csv")
    >>> exists(fn)
    True

    In the next example whe find the 'UnitSphere_500.csv'
    file in the 'csv' library and we check if it exists.
    But in this example it not exists.

    >>> fn = resource("csv", "UnitSphere_500.cvs")
    >>> exists(fn)
    False
    """
    import pkg_resources
    from os.path import join
    return pkg_resources.resource_filename(__name__,
                                           join("data",library,filename))

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

