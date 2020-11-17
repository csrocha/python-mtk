# -*- coding: ISO-8859-1 -*-
# $Id: tbl_ff.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Write Volumens to Mayavi files to visualization"""

from numpy import array
import cPickle
from enthought.tvtk.api import tvtk
from enthought.mayavi.api import Engine
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.mayavi.modules.outline import Outline
from enthought.mayavi.modules.image_plane_widget import ImagePlaneWidget
from enthought.mayavi.modules.iso_surface import IsoSurface

class writer:
    """Write pynum arrays to file as an scene

    TODO: Doesnt work. I think is not really usefull.

    > >> A = array([[[0.0, 1.0, 0.0],[0.0, 1.0, 0.0],[0.0, 1.0, 0.0]]])
    > >> w = writer(A)
    > >> w.write('test.vti')

    """
    def __init__(self, A):
        self.A = A

    def write(self, filename):
        spoints = tvtk.StructuredPoints(origin=(0,0,0), spacing=(1,1,1),
                                    dimensions=self.A.shape)
        spoints.point_data.scalars = self.A.transpose().flatten()
        spoints.point_data.scalars.name = "scalars"

        e = Engine()
        e.start()
        scene = e.new_scene()
        src = VTKDataSource(data = spoints)
        e.add_source(src)
        e.add_module(Outline())
        e.add_module(ImagePlaneWidget())
        e.add_module(IsoSurface())

        fobj = open(filename, 'w')
        cPickle.dump(scene, fobj)
        fobj.close()

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

