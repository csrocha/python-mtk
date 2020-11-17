# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
# Primitive Geometry Arc Object
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

from math import copysign
from numpy import array, zeros, arange, float, abs, argmax, all, ndenumerate,round, dot, outer, arccos, pi, inner, ndarray
from numpy.linalg import det, norm
from mtk.geometry.vol import BasicVolume
from mtk.geometry.iter import face, line
from mtk.geometry.planecut import cutfaces
from mtk.geometry.triangle import triangle
from mtk.geometry.sphere import sphere
from mtk.geometry.intersection import triangle_sphere
from mtk.geometry.line import segment
from mtk.geometry.arc import arc

from hashlib import sha1
from numpy import array_repr, allclose

class hndarray(ndarray):
    def __new__(cls, values):
        this = ndarray.__new__(cls, shape=values.shape, dtype=values.dtype, buffer=values.data)
        return this

    def __init__(self, values):
        s = array_repr(self).replace(' ','')
        self.__hash = s.__hash__()

    def __eq__(self, other):
        return allclose(self, other)

    def __hash__(self):
        return self.__hash

    def __setitem__(self, key, value):
        raise Exception('hashable arrays are read-only')

class polygon:
    def __init__(self, vertexs, faces):
        self.v = array(vertexs, float)
        self.f = array(faces, float)

        # Map of Vectors to Faces
        self.vtof = dict([(i, set()) for i in range(self.v.shape[0])])

        for f in range(len(self.f)):
            face = self.f[f]
            for i in range(len(face)):
                self.vtof[face[i]].add(f)

        # Calculate neighbourhoods vectors
        self.nb = dict([ (v, self.neighbourhoods(v)) for v in range(self.v.shape[0]) ])

    def neighbourhoods(self, vertex_id):
        faces = self.vtof[vertex_id]
        pvertexs = set([ v for f in faces for v in self.f[f] ])

        n = set()
        for v in pvertexs:
            if len(faces.intersection(self.vtof[v])) == 2:
                   n.add(v)

        return n

    def _search_vertexs_near_to_sphere(self, vertex_id, radix):
        """
        Search minimal vertexs set outside the sphere
        
        >>> v = [[0,0,0],[10,0,0],[10,10,0],[0,10,0],[-10,10,0],[-10,0,0],[-10,-10,0],[0,-10,0]]
        >>> f = [[0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6],[0,6,7],[0,7,1]]
        >>> P = polygon(v, f)
        >>> P._search_vertexs_near_to_sphere(0, 5)
        set([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        """
        S = sphere(self.v[vertex_id], radix)

        # Search minimal vertex set outside the sphere
        i_depth = set([vertex_id])
        i_visited = set()
        i_result = set()

        while i_depth:
            i = i_depth.pop()
            i_visited.add(i)
            for k in self.neighbourhoods(i):
                if k not in i_visited and S.dist(self.v[k]) < 0:
                    i_depth.add(k)
                elif k not in i_visited:
                    i_result.add(k)

        return i_result

    def area(self, vertex_id, radix, debug=1):
        """
        >>> v = [[0,0,0],[10,0,0],[10,10,0],[0,10,0],[-10,10,0],[-10,0,0],[-10,-10,0],[0,-10,0],[10,-10,0]]
        >>> f = [[0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6],[0,6,7],[0,7,8],[0,8,1]]

        >>> P = polygon(v, f)
        >>> P.area(0, 5)

        >>> from arc import circle
        >>> C = circle([0,0,1], [0,0,0], 5)
        >>> C.area()

        >>> P = polygon(v, f)
        >>> P.area(0, 3)
        >>> from arc import circle
        >>> C = circle([0,0,1], [0,0,0], 3)
        >>> C.area()

        """
        S = sphere(self.v[vertex_id], radix)

        faces = self._search_vertexs_near_to_sphere(vertex_id,radix)

        # Search faces witch intersept the sphere
        fs = set([ f for v in faces for f in self.vtof[v] ])
        AT = 0 # Termino del angulo de las tangentes
        KT = 0 # Termino de curvatura
        ta = {}
        tb = {}
        for f in fs:
            T = triangle(*(map(lambda i: self.v[i], self.f[f])))
            I = triangle_sphere(T, S)
            for alpha,oc,rc,a,b in I:
                A = arc(oc, rc, a, b, T.normal())
                # Calculo de la torsion de la tangente (Fase I)
                ta[hndarray(a)] = A.circle().tangent(a)
                tb[hndarray(b)] = A.circle().tangent(b)
                # Curvatura del arco
                kg = A.curvature()
                KT += A.length()*kg

        AT = sum([ arccos(min(1., dot(ta[k], tb[k]))) for k in ta.keys() ])

        if debug:
            from mtk.io.vtk_ff import writer
            w = writer(points=self.v, polygons=self.f, scalars=range(len(self.v)))
            w.write('test_poly')

            l = len(ta.keys())
            keys = ta.keys()

            ps = array([ k.T for k in keys ])

            w = writer(points=ps, lines=[(i,(i+1) % l) for i in range(l)], vectors=[ta[v] for v in keys])
            w.write('test_tang_a')

            w = writer(points=ps, lines=[(i,(i+1) % l) for i in range(l)], vectors=[tb[v] for v in keys])
            w.write('test_tang_b')

            w = writer(points=ps, lines=[(i,(i+1) % l) for i in range(l)], scalars=[arccos(min(1.,dot(ta[v],tb[v]))) for v in keys])
            w.write('test_tang_d')

        return ( 2*pi - AT - KT ) * S.r**2

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())


# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:
