# -*- coding: ISO-8859-1 -*-
# $Id: vol.py 89 2009-06-09 16:35:56Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""View tools for volume data

    Testing simple atoms draw


>>> from mtk.geometry.vol import create_from_coords
>>> import sys
>>> sys.path.append("/mnt/data/Projects/protein_docking/trunk/python-mtk/")
>>> from mtk.geometry.vol import Volume
>>> from numpy import array
>>> a = array([[ 0, 0, 0 ], [ -1, 2, -3], [5, -6, 7]])
>>> V = create_from_coords(a, resolution=1.7, inc=3.0)
>>> V.put_ball(a[0], 2.0, lambda r,v: r)
>>> V.put_ball(a[1], 2.0, lambda r,v: r)
>>> V.put_ball(a[2], 2.0, lambda r,v: r)

>>> View = Viewer()
>>> View.append(V)
>>> View.run()

    More complex atoms draw

>>> from mtk import resource
>>> from mtk.io.pdb_ff import reader
>>> p = reader(open(resource("test","obj01.pdb")))
>>> C = p.get_coords()
>>> V = create_from_coords(C, resolution=0.5, inc=3.0)
>>> for c in C: V.put_ball(c, 2.0, lambda r,v: abs(r)<3 and 1 or v)

>>> View = Viewer()
>>> View.append(V)
>>> View.run()

"""
from numpy import array, ogrid, sin, ravel
from enthought.traits.api import HasTraits, Instance, Property
from enthought.traits.ui.api import View, Item, HSplit, VSplit, InstanceEditor
from enthought.tvtk.api import tvtk
from enthought.tvtk.array_handler import get_vtk_array_type
from enthought.tvtk.pyface.scene_model import SceneModel
from enthought.tvtk.pyface.scene_editor import SceneEditor
from enthought.mayavi import mlab
from enthought.mayavi.api import Engine
from enthought.mayavi.plugins.app import Mayavi
from enthought.mayavi.sources.array_source import ArraySource
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene
from enthought.mayavi.core.ui.engine_view import EngineView
from enthought.mayavi.tools.mlab_scene_model import \
                    MlabSceneModel
from enthought.mayavi.modules.api import \
        Outline, Surface, ImagePlaneWidget, IsoSurface
from enthought.mayavi.filters.api import \
        Contour, SetActiveAttribute

from mtk import have_mayavi

if not have_mayavi():
    raise ImportError("Mayavi no está instalado. Instalelo para utilizar este módulo")

class Viewer(HasTraits):
    scene = Instance(MlabSceneModel, ())

    engine_view = Instance(EngineView)

    current_selection = Property

    view = View(HSplit(Item(name='engine_view',
                               style='custom',
                               resizable=True,
                               show_label=False
                               ),
                           Item(name='scene',
                                editor=SceneEditor(),
                                show_label=False,
                                resizable=True,
                                height=500,
                                width=500),
                    ),
            resizable=True,
            scrollable=True
            )

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.engine_view = EngineView(engine=self.scene.engine)
        self.scene.engine.on_trait_change(self._selection_change,
                                          'current_selection')

    def append(self, vol_scalar, vol_vector=None):
        e = self.scene.engine
        s = vol_scalar._data.real.transpose().copy()
        if vol_vector != None:
            V = vol_vector._data.real.transpose().copy()
            V = array([ (v,0,0) for v in V.flat ])
            V = V.reshape(vol_scalar._data.shape + (3,))
            src = ArraySource(origin=vol_scalar.min, spacing=vol_scalar.delta,
                              scalar_data = s, scalar_name='scalar',
                              vector_data = V, vector_name='vector')
        else:
            src = ArraySource(origin=vol_scalar.min, spacing=vol_scalar.delta,
                              scalar_data = s, scalar_name='scalar')
        e.add_source(src)
        #c = Contour()
        #e.add_filter(c)
        #c.add_module(Surface())
        e.add_module(Outline())
        e.add_module(ImagePlaneWidget())
        e.add_module(IsoSurface())
        pass

    def _selection_change(self, old, new):
        self.trait_property_changed('current_selection', old, new)

    def _get_current_selection(self):
        return self.scene.engine.current_selection

    def run(self):
        self.configure_traits()

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

