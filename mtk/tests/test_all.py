# -*- coding: ISO-8859-1 -*-
# $Id: __version__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

import mtk
import mtk.signal
import mtk.inferator
import mtk.forcefield
import mtk.dock
import mtk.energy
import mtk.energy.electrostatic
import mtk.energy.geocomp
import mtk.energy.lj
import mtk.energy.solvation
import mtk.geometry
import mtk.geometry.neighbours
import mtk.geometry.kdtree
import mtk.geometry.rms
import mtk.geometry.vol
import mtk.io
import mtk.io.pdb_ff
import mtk.io.tbl_ff
import mtk.io.prm_ff
import mtk.io.pqr_ff
import mtk.sqlfunctions
import mtk.sqlfunctions.geometry
import mtk.sqlfunctions.set
import mtk.view
import mtk.storage
import mtk.signal

import unittest
import doctest

from mtk import have_gpu, have_mayavi, have_vtk

modules = [ mtk,
            mtk.signal,
            mtk.inferator,
            mtk.forcefield,
            mtk.storage,
            mtk.dock,
            mtk.energy,
            mtk.energy.electrostatic,
            mtk.energy.lj,
            mtk.energy.geocomp,
            mtk.energy.solvation, # Have TODO!
            mtk.geometry,
            mtk.geometry.kdtree,
            mtk.geometry.neighbours,
            mtk.geometry.rms,
            mtk.geometry.transformation,
            mtk.geometry.vol,
            mtk.io,
            mtk.io.pdb_ff, # Download take time
            mtk.io.pqr_ff,
            mtk.io.prm_ff,
            mtk.io.tbl_ff, # Will not use in future
            mtk.sqlfunctions,
            mtk.sqlfunctions.geometry,
            mtk.sqlfunctions.set,
            mtk.view,
          ]

if have_gpu():
    import mtk.geometry.vol_cuda
    modules.append(mtk.geometry.vol_cuda)
else:
    print "Cuda support will not be tested."

if have_mayavi():
    import mtk.view.vol
    import mtk.io.mayavi_ff
    modules.append(mtk.view.vol) # Complex testing
    modules.append(mtk.io.mayavi_ff)
else:
    print "Mayavi support will not be tested."

if have_vtk():
    import mtk.io.vtk_ff
    modules.append(mtk.io.vtk_ff)
else:
    print "Vtk support will not be tested."

def getTestSuite():
    suite = unittest.TestSuite()
    for mod in modules:
        try:
            suite.addTest(mod.test_suite())
        except AttributeError, r:
            print "In module %s:" % mod.__file__
            print r
            break
    return suite

