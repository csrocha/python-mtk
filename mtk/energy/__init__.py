# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Inicializador del modulo de cálculo de energias.

>>> from mtk.geometry.vol import create_from_coords, Volume, BasicVolume
>>> from mtk.geometry.rms import array_rmsd
>>> import mtk.energy.electrostatic_c as c_elec
>>> c_calc = c_elec.calc
>>> from mtk.energy.electrostatic import calc as py_calc
>>> from numpy import array, zeros, sum, complex, allclose
>>> from time import time


  Calcula el volumen de electrostática usando cálculo en C y en Python
  sobre dos átomos de fantasia (0,0,0) y (1.5,1,-1.5) con cargas 1 y 2
  respectivamente. Los volúmenes (grillas) tienen un delta de 0.5.

>>> A = array([[ 0.0, 0.0, 0.0, 1.0, 1.0 ], [ 1.5, 1.0, -1.5, 2.0, 2.0 ]])
>>> A1 = zeros((5,5,5), dtype=complex)
>>> V1 = BasicVolume((-1,-1,-1), (0.5,0.5,0.5), A1)
>>> A2 = zeros((5,5,5), dtype=complex)
>>> V2 = BasicVolume((-1,-1,-1), (0.5,0.5,0.5), A2)
>>> c_calc(V1,A)
>>> py_calc(V2,A)

> >> _view(V1, V1._data, V2._data)

>>> allclose(array_rmsd(V1._data, V2._data), 0.000)
True

>>> a = array([[ 0.0, 0.0, 0.0, 1.0, 1.0 ]]*2)
>>> V1 = create_from_coords(a[1:4], resolution=1.0, inc=10.0, dtype=complex)
>>> V2 = create_from_coords(a[1:4], resolution=1.0, inc=10.0, dtype=complex)
>>> c_t = time()
>>> c_calc(V1, a)
>>> c_t = time() - c_t
>>> py_t = time()
>>> py_calc(V2, a)
>>> py_t = time() - py_t
>>> all(V1.shape == V2.shape)
True
>>> allclose(array_rmsd(V1._data, V2._data), 0.000)
True
>>> (c_t * 100) / py_t < 0.1
True

>>> atoms = array( [[ 0., 0., 0., 1.],[1., 1., 1., 1.]] )
>>> py_t = time()
>>> Vpy = create_geometry_vol_py( atoms[:,:3], atoms[:,3] )
>>> py_t = time() - py_t
>>> c_t = time()
>>> Vc = create_geometry_vol_c( atoms[:,:3], atoms[:,3] )
>>> c_t = time() - c_t
>>> (c_t * 100) / py_t < 5.2
True

>>> [ sum(Vpy._data.flatten() == i) for i in [ 1., 9j ] ] == [ sum(Vc._data.flatten() == i) for i in [ 1., 9j ] ]
True

"""
from numpy import array
from mtk.geometry.vol import BasicVolume
from mtk.geometry import centre

def create_geometry_vol_py(coords, vdw, size=None, resolution=1.0, solv=0, core=9j, surf=1, r_solv=1):
    """
    Create a volume for geometry complementary docking.

    >>> from numpy import array
    >>> atoms = array( [[ 0., 0., 0., 1.]] )
    >>> V = create_geometry_vol_py( atoms[:,:3], atoms[:,3] )
    >>> [ sum(V._data.flatten() == i) for i in [ 0., 1., 9j ] ]
    [92, 26, 7]

    >>> atoms = array( [[ 0., 0., 0., 1.],[5., 5., 5., 1.]] )
    >>> V = create_geometry_vol_py( atoms[:,:3], atoms[:,3] )
    >>> [ sum(V._data.flatten() == i) for i in [ 1., 9j ] ]
    [52, 14]

    >>> from time import time
    >>> atoms = array( [[ 0., 0., 0., 1.]] * 100 )
    >>> t = time()
    >>> V = create_geometry_vol_py( atoms[:,:3], atoms[:,3] )
    >>> time() - t < 2
    True

    """
    from numpy import ones, float
    from mtk.energy.geocomp import geometry_complementary as sas_f
    from mtk.geometry.vol import create_from_coords

    if coords.shape[0] == 3:
        coords = coords.T
    elif coords.shape[1] == 3:
        None
    else:
        raise RuntimeError("Coordinates must be 3xn or nx3 coordinates. If its 3x3 its take row as points.")

    if size == None:
        inc = max(vdw)+r_solv
        V = create_from_coords(coords, resolution=resolution, inc=inc,
                               init_array=ones, dtype=complex)
    else:
        size = size + 4*(max(vdw)+r_solv)
        newmin = centre(coords) - size/2
        V = BasicVolume(newmin, array([resolution]*3), ones(size/resolution,
                                                            dtype=complex))
    V._data = V._data * solv

    F = lambda r, v, r_core, r_solv: sas_f(r, v,
                                           core=core, surf=surf, solv=solv,
                                           r_core=r_core, r_solv=r_solv)
    for i in xrange(coords.shape[0]):
        x = coords[i]
        radius = vdw[i]
        V.put_ball(x, radius + r_solv, lambda r, v: F(r, v, radius, radius + r_solv))
    return V

def create_geometry_vol_c(coords, vdw, size=None, resolution=1.0, solv=0, core=9j, surf=1, r_solv=1):
    """
    Create a volume for geometry complementary docking.

    >>> from numpy import array
    >>> atoms = array( [[ 0., 0., 0., 1.]] )
    >>> V = create_geometry_vol_c( atoms[:,:3], atoms[:,3] )
    >>> [ sum(V._data.flatten() == i) for i in [ 0., 1., 9j ] ]
    [92, 26, 7]

    >>> atoms = array( [[ 0., 0., 0., 1.],[5., 5., 5., 1.]] )
    >>> V = create_geometry_vol_c( atoms[:,:3], atoms[:,3] )
    >>> [ sum(V._data.flatten() == i) for i in [ 1., 9j ] ]
    [52, 14]

    >>> from time import time
    >>> atoms = array( [[ 0., 0., 0., 1.]] * 100 )
    >>> t = time()
    >>> V = create_geometry_vol_c( atoms[:,:3], atoms[:,3] )
    >>> time() - t < 0.1
    True

    """
    from numpy import ones, float
    from mtk.energy.geocomp import geometry_complementary as sas_f
    from mtk.geometry.vol import create_from_coords
    from mtk.geometry.vol_c import put_ball

    if coords.shape[0] == 3:
        coords = coords.T
    elif coords.shape[1] == 3:
        None
    else:
        raise RuntimeError("Coordinates must be 3xn or nx3 coordinates. If its 3x3 its take row as points.")

    if size == None:
        inc = max(vdw)+r_solv
        V = create_from_coords(coords, resolution=resolution, inc=inc,
                               init_array=ones, dtype=complex)
    else:
        size = size + 4*(max(vdw)+r_solv)
        newmin = centre(coords) - size/2
        V = BasicVolume(newmin, array([resolution]*3), ones(size/resolution,
                                                            dtype=complex))
    V._data = V._data * solv

    order_list = [ solv, surf, core ];
    # TODO: Pasar este código a C
    for i in xrange(coords.shape[0]):
        x = coords[i]
        radius = vdw[i]
        assign_map = { radius: core, radius + r_solv: surf };
        put_ball(V, x, assign_map, order_list)
    return V

def _view(V, A, B):
    """
    >>> from mtk.geometry.rms import array_rmsd
    >>> from numpy import array, allclose
    >>> from mtk.io.pqr_ff import reader
    >>> from mtk import resource
    >>> parms = { 'resolution': 0.5, 'solv': 0, 'core': 9j, 'surf': 1, 'r_solv': 1 }

    >>> atoms = array( [[ 0., 0., 0., 1.]] )
    >>> Vpy = create_geometry_vol_py(atoms[:,:3], atoms[:,3], **parms)
    >>> Vc = create_geometry_vol_c(atoms[:,:3], atoms[:,3], **parms)

    > >> _view(Vpy, Vpy._data, Vc._data)

    >>> Rd = reader(open(resource("test", "1A2K_l_b.pqr")))
    >>> C = array([ (i['x'], i['y'], i['z'], 1., i['occupancy']) for h,i in Rd if h == 'ATOM'])
    >>> Vpy = create_geometry_vol_py(C[:30,:3], C[:30,4], **parms)
    >>> Vc = create_geometry_vol_c(C[:30,:3], C[:30,4], **parms)

    > >> _view(Vpy, Vpy._data, Vc._data)

    >>> allclose(Vpy._data, Vc._data)
    True
    >>> array_rmsd(Vpy._data, Vc._data) == 0
    True

    >>> [ sum(Vpy._data.flatten() == i) for i in [ 0, 1, 9j ] ]
    [13047, 2739, 2034]
    >>> [ sum(Vc._data.flatten() == i) for i in [ 0, 1, 9j ] ]
    [13047, 2739, 2034]
    """
    from mtk.geometry.vol import BasicVolume
    from enthought.mayavi.plugins.app import Mayavi
    from mtk.view.vol import view
    ViewMayavi = Mayavi()
    mayavi = ViewMayavi.script
    view(BasicVolume(V.min, V.delta, A), mayavi)
    view(BasicVolume(V.min, V.delta, B), mayavi)
    ViewMayavi.main()

"""
Declaro como la función principal para realizar geometrias
aquella programada en C que es mucho más rápida.
"""
create_geometry_vol = create_geometry_vol_c

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

