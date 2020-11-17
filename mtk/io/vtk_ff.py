# -*- coding: ISO-8859-1 -*-
# $Id: tbl_ff.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Write Volumens to VTK files to visualization"""

import vtk

from numpy import array, allclose, transpose
from tvtk.api import tvtk
from tvtk.array_handler import get_vtk_array_type

def _cell_to_list(data):
        res = []
        c = 0
        while c < len(data):
            l = data[c]
            res.append([])
            for j in range(0,l):
                c += 1
                res[-1].append(data[c])
            c += 1
        return res

class reader:
    """
    >>> P = array([[[0.0, 1.0, 0.0],[0.0, 1.0, 0.0],[0.0, 1.0, 0.0]]])
    >>> w = writer(points=P)
    >>> w.write('test1')

    >>> r = reader()
    >>> D, I = r.read('test1.vti')
    >>> I.dimensions
    array([1, 3, 3])
    >>> I.spacing
    array([ 1.,  1.,  1.])
    >>> I.origin
    array([ 0.,  0.,  0.])
    >>> allclose(D['scalars'], array([[[0.0, 1.0, 0.0],[0.0, 1.0, 0.0],[0.0, 1.0, 0.0]]]))
    True

    >>> points = array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], 'f')
    >>> triangles = array([[0,1,3], [0,3,2], [1,2,3], [0,2,1]])
    >>> temperature = array([10., 20., 30., 40.])
    >>> w = writer(points=points,polygons=triangles,scalars=temperature)
    >>> w.write('test2')

    >>> r = reader()
    >>> D, I = r.read('test2.vtp')
    >>> allclose(D['points'], points)
    True
    >>> allclose(D['polys'], triangles)
    True

    """
    def __init__(self):
        pass

    def read(self, filename):
        extension = filename[-3:]
        if extension=='vtp':
            return self.readPolyData(filename)
        elif extension=='vti':
            return self.readStructuredPoints(filename)
        return None

    def readPolyData(self, filename):
        reader = tvtk.XMLPolyDataReader(file_name=filename)
        reader.update()
        polydata = reader.output
        data = {'scalars': polydata.point_data.scalars,
                'points': polydata.points,
                'polys': _cell_to_list(polydata.polys.to_array()),
               }
        return data, polydata
         
    def readStructuredPoints(self, filename):
        reader = tvtk.XMLImageDataReader(file_name=filename)
        reader.update()
        image = reader.output
        data = image.point_data.scalars.to_array()
        data = { 'scalars': transpose(data.reshape(image.dimensions),(0,2,1)) }
        return data, image

    def read_as_vol(self, filename):
        D, I = self.read(filename)
        return BasicVolume(I.origin, I.spacing, D['scalars'])

class writer:
    """Write pynum arrays to file as an StructuredPoints

    >>> P = array([[[0.0, 1.0, 0.0],[0.0, 1.0, 0.0],[0.0, 1.0, 0.0]]])
    >>> w = writer(points=P)
    >>> w.write('test1')

    >>> points = array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], 'f')
    >>> triangles = array([[0,1,3], [0,3,2], [1,2,3], [0,2,1]])
    >>> temperature = array([10., 20., 30., 40.])
    >>> w = writer(points=points,polygons=triangles,scalars=temperature)
    >>> w.write('test2')

    """
    def __init__(self, origin=None, spacing=None,
                 points=None,
                 vertices=None, lines=None, polygons=None, strips=None,
                 scalars=None, vectors=None):
        self.points = array(points, float)
        self.lines = lines
        self.vertices = vertices
        self.polygons = polygons
        self.strips = strips
        self.origin = (0,0,0) if origin is None else origin
        self.spacing = (1,1,1) if spacing is None else spacing
        self.scalars = array(scalars, float) if scalars is not None else None
        self.vectors = array(vectors, float) if vectors is not None else None

    def write(self, filename):
        isPolydata = [None, None, None, None] != [ self.vertices, self.lines, self.polygons, self.strips ]

        if self.points is not None and not isPolydata:
            return self.writeStructuredPoints(filename)
        elif self.points is not None and isPolydata:
            return self.writePolyData(filename)
        return False

    def writePolyData(self, filename):
        vertices = array(self.vertices, int) if self.vertices is not None else None
        lines = array(self.lines, int) if self.lines is not None else None
        polys = array(self.polygons, int) if self.polygons is not None else None
        strips = array(self.strips, int) if self.strips is not None else None
        P = tvtk.PolyData(points=array(self.points, float),
                          verts=vertices,
                          lines=lines,
                          polys=polys,
                          strips=strips)
        if not self.scalars is None:
            P.point_data.scalars = self.scalars
            P.point_data.scalars.name = "scalars" 
        if not self.vectors is None:
            P.point_data.vectors = self.vectors
            P.point_data.vectors.name = "vectors" 
        if filename[-4:] != '.vtp':
            filename += '.vtp'
        w = tvtk.XMLPolyDataWriter(input=P, file_name=filename)
        w.write()

    def writeStructuredPoints(self, filename):
        spoints = tvtk.StructuredPoints(origin=self.origin, \
                                        spacing=self.spacing, \
                                        dimensions=self.points.shape)
        spoints.point_data.scalars = self.points.transpose().flatten()
        spoints.point_data.scalars.name = "scalars"
        spoints.scalar_type = get_vtk_array_type(self.points.dtype)
        if filename[-4:] != '.vti':
            filename += '.vti'
        w = tvtk.XMLImageDataWriter(input=spoints, file_name=filename)
        w.write()

def read_vol(filename):
    from mtk.geometry.vol import BasicVolume

    if filename[-4:] != '.vti':
        return False

    R = reader()
    D, I = R.read(filename)
    return BasicVolume(I.origin, I.spacing, D['scalars'])

def write_vol(vol, filename):
    from mtk.geometry.vol import BasicVolume

    if filename[-4:] != '.vti':
        return False

    W = writer(origin = vol.min, spacing = vol.delta, points = vol._data)
    W.write(filename)

def read_pol(filename):
    from mtk.geometry.polygon import polygon

    if filename[-4:] != '.vtp':
        return False

    R = reader()
    D, I = R.read(filename)
    return polygon(D['points'], D['polys']), D['scalars'] if 'scalars' in D else None

def write_pol(pol, filename, scalars=None):
    from mtk.geometry.polygon import polygon

    if filename[-4:] != '.vtp':
        return False

    W = writer(points = pol.v, polygons = pol.f, scalars = scalars)
    W.write(filename)
    return True

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:

