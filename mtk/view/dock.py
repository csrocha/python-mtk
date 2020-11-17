# -*- coding: ISO-8859-1 -*-
# $Id: vol.py 89 2009-06-09 16:35:56Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""View tools for dock results

"""
from numpy import array, ogrid, sin, ravel, nonzero
from enthought.traits.api import HasTraits, Instance, Property, Range, on_trait_change
from enthought.traits.ui.api import View, Item, VGroup, HSplit, VSplit, InstanceEditor
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
        Outline, Surface, ImagePlaneWidget, IsoSurface, Text

from mtk import have_mayavi

if not have_mayavi():
    raise ImportError("Mayavi no está instalado. Instalelo para utilizar este módulo")

class Viewer(HasTraits):
    """
    >>> from numpy import array, argmax, argmin, all, identity, allclose, zeros
    >>> from mtk.geometry.vol import BasicVolume
    >>> from mtk.dock import dock
    >>> c, s = 9j, 1

    Testing del docking en el eje Z

    >>> Ad = array([ \
                   [ \
                    [ c, s ],\
                    [ 0, 0 ],\
                   ],[ \
                    [ 0, 0 ],\
                    [ 0, 0 ],\
                   ] \
                  ], dtype=complex)
    >>> Bd = array([ \
                   [ \
                    [ s, c ],\
                    [ 0, 0 ],\
                   ],[ \
                    [ 0, 0 ],\
                    [ 0, 0 ],\
                   ] \
                  ], dtype=complex)
    >>> A = zeros((4,4,4), dtype=complex)
    >>> B = zeros((4,4,4), dtype=complex)
    >>> A[1:3,1:3,1:3] = Ad
    >>> B[1:3,1:3,1:3] = Bd
    >>> Va = BasicVolume((0,0,0), (1,1,1), A)
    >>> Vb = BasicVolume((0,0,0), (1,1,1), B)
    >>> fVb = lambda t: Vb
    >>> D = dock(Va, fVb, [ identity(4) ])
    >>> V = Viewer(D, Va, Vb)
    >>> V.run()
    >>> del V

    Testing del docking en el eje Y

    >>> Ad = array([ \
                   [ \
                    [ c, 0 ],\
                    [ s, 0 ],\
                   ],[ \
                    [ 0, 0 ],\
                    [ 0, 0 ],\
                   ] \
                  ], dtype=complex)
    >>> Bd = array([ \
                   [ \
                    [ s, 0 ],\
                    [ c, 0 ],\
                   ],[ \
                    [ 0, 0 ],\
                    [ 0, 0 ],\
                   ] \
                  ], dtype=complex)
    >>> A = zeros((4,4,4), dtype=complex)
    >>> B = zeros((4,4,4), dtype=complex)
    >>> A[1:3,1:3,1:3] = Ad
    >>> B[1:3,1:3,1:3] = Bd
    >>> Va = BasicVolume((0,0,0), (1,1,1), A)
    >>> Vb = BasicVolume((0,0,0), (1,1,1), B)
    >>> fVb = lambda t: Vb
    >>> D = dock(Va, fVb, [ identity(4) ])
    >>> V = Viewer(D, Va, Vb)
    >>> V.run()
    >>> del V

    Testing del docking en el eje Z

    >>> Ad = array([ \
                   [ \
                    [ c, 0 ],\
                    [ 0, 0 ],\
                   ],[ \
                    [ s, 0 ],\
                    [ 0, 0 ],\
                   ] \
                  ], dtype=complex)
    >>> Bd = array([ \
                   [ \
                    [ s, 0 ],\
                    [ 0, 0 ],\
                   ],[ \
                    [ c, 0 ],\
                    [ 0, 0 ],\
                   ] \
                  ], dtype=complex)
    >>> A = zeros((4,4,4), dtype=complex)
    >>> B = zeros((4,4,4), dtype=complex)
    >>> A[1:3,1:3,1:3] = Ad
    >>> B[1:3,1:3,1:3] = Bd
    >>> Va = BasicVolume((0,0,0), (1,1,1), A)
    >>> Vb = BasicVolume((0,0,0), (1,1,1), B)
    >>> fVb = lambda t: Vb
    >>> D = dock(Va, fVb, [ identity(4) ])
    >>> V = Viewer(D, Va, Vb)
    >>> V.run()
    >>> del V

    """
    scene = Instance(MlabSceneModel, ())

    engine_view = Instance(EngineView)

    current_selection = Property

    view = View(VSplit(
                    HSplit(Item(name='engine_view',
                               style='custom',
                               resizable=True,
                               show_label=False
                               ),
                           Item(name='scene',
                                editor=SceneEditor(),
                                show_label=False,
                                resizable=True,
                                height=400,
                                width=500),
                    ),
                    VGroup('_', 'result', 'score'),
            ),
            resizable=True,
            scrollable=True
            )

    def __init__(self, (rR, rT, rS), R, L, S=None, **traits):
        assert(len(rR) == len(rT))
        assert(len(rT) == len(rS))
        rSR = [ float(i.real) for i in rS ]
        self.add_trait('result', Range(0,len(rR)-1,0))
        self.add_trait('score', Range(min(rSR),max(rSR),rSR[0]))
        HasTraits.__init__(self, **traits)
        self.engine_view = EngineView(engine=self.scene.engine)
        self.scene.engine.on_trait_change(self._selection_change,
                                          'current_selection')

        self.results = {'R': rR, 'T':rT, 'S': rSR}
        self.vol = {}
        self.src = {}
        self.moviles = []
        self.reference = 'receptor'
        self.update = True

        self.append(R, 'receptor')
        if S != None: self.append(S, 'solution')
        self.append(L, 'ligand', movil=True)

        self.put_text('')

        self.update_plot()

    def put_text(self, txt):
        self.text = Text(text=txt)
        self.scene.engine.add_module(self.text)
        pass

    def append(self, S, name, movil=False):
        self.vol[name]= S
        self.src[name] = ArraySource()
        self.src[name].scalar_data = S._data.real #.transpose()
        self.src[name].origin = S.min
        self.src[name].spacing = S.delta

        e = self.scene.engine
        e.add_source(self.src[name])
        e.add_module(Outline())
        e.add_module(ImagePlaneWidget())
        e.add_module(IsoSurface())

        if movil : self.moviles.append(name)

    @on_trait_change('score')
    def update_score(self):
        if self.update == True:
            idx = nonzero(array(self.results['S']) >= self.score)
            if len(idx) > 0 and len(idx[0]) != 0:
                self.result =  int(idx[0][0])

    @on_trait_change('result')
    def update_plot(self):
        if self.update == True:
            self.update = False
            d = self.results['T'][self.result]
            self.text.text = "(%f %f %f)" % tuple(d)
            ref = self.vol[self.reference]
            for movil in self.moviles:
                mov = self.vol[movil]
                self.src[movil].origin = mov.min + (ref.center + d - mov.center)
            self.score = float(self.results['S'][self.result])
            self.update = True
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

