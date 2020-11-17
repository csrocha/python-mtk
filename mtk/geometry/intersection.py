# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
# Primitive Geometry Object Intersections
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

from numpy import empty_like, dot, array, ndarray, seterr, allclose, sqrt, arccos, pi
from numpy.linalg import norm
from line import segment, line
from sphere import sphere
from triangle import triangle, plane
from bitarray import bitarray

def _perp( a ) :
    b = empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

def line_line(l1, l2):
    """
    Line intersection check

    >>> l1 = segment( [0.0, 0.0], [1.0, 0.0] ).line()
    >>> l2 = segment( [4.0, -5.0], [4.0, 2.0] ).line()
    >>> line_line(l1, l2)
    array([ 4.,  0.])
    
    >>> l1 = segment( [2.0, 2.0,1.0], [4.0, 3.0,1.0] ).line()
    >>> l2 = segment( [6.0, 0.0,1.0], [6.0, 3.0,1.0] ).line()
    >>> line_line(l1,l2)
    array([ 6.,  4.,  1.])

    >>> l1 = segment( [1.0, 0.0,0.0], [0.0, 0.0,0.0] ).line()
    >>> l2 = segment( [1.0, 1.0,0.0], [0.0, 1.0,0.0] ).line()
    >>> line_line(l1,l2)
    False
    """
    olderr = seterr(all='ignore')

    da = l1.d
    db = l2.d
    dp = l1.o-l2.o
    dap = _perp(da)
    denom = dot( dap, db)
    if denom == 0.0: return False
    num = dot( dap, dp )

    seterr(**olderr)

    return (num / denom)*db + l2.o

def segment_segment(s1, s2):
    """
    Segment intersection.

    >>> s1 = segment( [1.0, 0.0,0.0], [0.0, 0.0,0.0] )
    >>> s2 = segment( [1.0, 1.0,0.0], [0.0, 1.0,0.0] )
    >>> segment_segment(s1,s2)
    False

    >>> s1 = segment( [2.0, 2.0,1.0], [4.0, 3.0,1.0] )
    >>> s2 = segment( [6.0, 0.0,1.0], [6.0, 3.0,1.0] )
    >>> segment_segment(s1,s2)
    False

    >>> s1 = segment( [2.0, 2.0,1.0], [8.0, 5.0,1.0] )
    >>> s2 = segment( [6.0, 0.0,1.0], [6.0, 6.0,1.0] )
    >>> segment_segment(s1,s2)
    array([ 6.,  4.,  1.])
    """
    l1=s1.line()
    l2=s2.line()
    i = line_line(l1, l2)
    if isinstance(i, bool): return False
    k = s1.affine(i)
    return k >= 0 and k <= 1 and i

def line_sphere(l, s):
    """
    Solve two points of intersection of a line and circle

    >>> s = sphere([0, 0, 0], 10.0)
    >>> a = segment([0, 0, 0], [1, 0, 0]).line()
    >>> f, g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    >>> allclose([s.dist(f), s.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([1, 0, 0], [0, 0, 0]).line()
    >>> f, g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    >>> allclose([s.dist(f), s.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([1, 0, 0],[0, 1, 0]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    >>> allclose([s.dist(f), s.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([1, 1, 0],[0, 1, 0]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    >>> allclose([s.dist(f), s.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([0, 0, 0],[1, 0, 0]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([0, 0, 0],[0, 1, 0]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([0, 0, 0],[0, 0, 1]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
    
    >>> a = segment([1, 0, 0],[0, 0, 0]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([0, 1, 0],[0, 0, 0]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([0, 0, 1],[0, 0, 0]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True
 
    >>> a = segment([ 0., 0., 0.],[ 10.,-10.,  0.]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([ 10.,-10.,  0.],[  0.,-10. , 0.]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True

    >>> a = segment([  0.,-10. , 0.],[ 0., 0., 0.]).line()
    >>> f,g = line_sphere(a, s)
    >>> allclose([a.dist(f), a.dist(g)], [0.0, 0.0])
    True


    """
    o = l.o
    l = l.d
    c = s.o
    r = s.r
    LO = dot(l,o)
    LC = dot(l,c)
    OC = dot(o,c)
    A = LO - LC
    AA = A*A

    LL = dot(l,l)
    OO = dot(o,o)
    CC = dot(c,c)
    RR = r*r

    B = OO + CC - RR - 2*OC

    C = LL

    tsqr = AA - C*B

    if tsqr < 0:
        return tuple()

    tsqr = sqrt(tsqr)
    k1 = (-A + tsqr)/LL
    k2 = (-A - tsqr)/LL

    return (l*k1+o, l*k2+o)

def line_plane(l, p):
    """
    Intersection of line and a plane

    >>> l = segment([0,0,0],[0,1,0]).line()
    >>> p = triangle([1,0,0],[0,0,0],[0,0,1]).plane()
    >>> allclose(line_plane(l, p), [0,0,0])
    True

    >>> l = segment([1,0,0],[0,0,0]).line()
    >>> p = triangle([0,1,0],[0,0,0],[0,0,1]).plane()
    >>> allclose(line_plane(l, p), [0,0,0])
    True

    >>> l = segment([1,0,0],[-1,0,0]).line()
    >>> p = triangle([0,1,0],[0,0,0],[0,0,1]).plane()
    >>> allclose(line_plane(l, p), [0,0,0])
    True

    """
    d = dot((p.o - l.o), p.n) / dot(l.d, p.n)
    return l(d)

def segment_sphere(seg, sph):
    """
    Return the intersection points between segment and sphere.

    >>> s = sphere([0, 0, 0], 10.0)
    >>> a = segment([1, 1, 0],[0, 1, 0])
    >>> segment_sphere(a, s)
    []

    >>> a = segment([15, 1, 0],[-15, 1, 0])
    >>> ss = segment_sphere(a, s)
    >>> allclose(map(a.line().dist, ss), [0.,0.])
    True
    >>> allclose(map(s.dist, ss), [0.,0.])
    True

    >>> s = sphere([1, 1, 1], 10.0)
    >>> a = segment([15, 1, 1],[-15, 1, 0])
    >>> ss = segment_sphere(a, s)
    >>> allclose(map(a.line().dist, ss), [0.,0.])
    True
    >>> allclose(map(s.dist, ss), [0.,0.])
    True
 
    >>> s = sphere([0, 0, 0], 5.0)
    >>> a = segment([ 0., 0., 0.],[ 10.,-10.,  0.])
    >>> segment_sphere(a, s)
    [array([ 3.53553391, -3.53553391,  0.        ])]

    >>> a = segment([ 10.,-10.,  0.],[  0.,-10. , 0.])
    >>> segment_sphere(a, s)
    []

    >>> a = segment([  0.,-10. , 0.],[ 0., 0., 0.])
    >>> segment_sphere(a, s)
    [array([ 0., -5.,  0.])]

    """
    ints = line_sphere(seg.line(), sph)
    if ints:
        return [ a for a,i in zip(ints, map(seg.affine, ints)) if i >= 0 and i <= 1 ]
    return [] 

def _circle_in_sphere(p, c, r, s):
    ret = True
    for v in [ [1,0,0], [0,1,0], [0,0,1] ]:
        a = p.project(v)-c
        na = norm(a)
        if na == 0: continue
        an = a / na
        a = c + an*r
        ret &= allclose(p.dist(a), 0.)
        if not ret: import pdb; pdb.set_trace()
        ret &= allclose(s.dist(a), 0.)
        if not ret: import pdb; pdb.set_trace()
    return ret

def plane_sphere(p, s):
    """

    >>> p = triangle([0,0,0],[1,0,0],[0,0,1]).plane()
    >>> s = sphere([0,0,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([0,0,0],[1,0,0],[0,0,1]).plane()
    >>> s = sphere([0,1,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([0,0,0],[1,0,0],[0,0,1]).plane()
    >>> s = sphere([1,0,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([0,1,0],[1,1,0],[0,1,1]).plane()
    >>> s = sphere([0,0,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([0,-1,0],[1,-1,0],[0,-1,1]).plane()
    >>> s = sphere([0,0,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([0,0,0],[1,0,0],[0,1,0]).plane()
    >>> s = sphere([0,0,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([0,0,1],[1,0,1],[0,1,1]).plane()
    >>> s = sphere([0,0,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([1,0,0],[0,1,0],[0,0,1]).plane()
    >>> s = sphere([0,0,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([1,0,0],[0,1,0],[0,0,1]).plane()
    >>> s = sphere([0,0,0], 10)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    >>> p = triangle([-0.908749,-0.045065,1.024429],[-0.90862,-0.055086,1.024449],[-0.908307,-0.045029,1.038119]).plane()
    >>> s = sphere([-0.923142, -0.045219,  0.997029], 0.02)
    >>> c,r = plane_sphere(p, s)
    >>> _circle_in_sphere(p, c, r, s)
    True

    """

    p.normalize()

    d = dot(s.o-p.o, p.n)

    if d > s.r:
        return False
    else:
        return (s.o - d*p.n, sqrt(s.r*s.r - d*d))

def _arc_(a, b, c, r):
    return (arccos(dot((a - c) / r, (b - c) / r)),c,r,a,b)

_tim_ = {
    # Caso Primero
    '110': [(0,0,1,0)],
    '011': [(1,0,2,0)],
    '101': [(2,0,0,0)],
    # Caso Segundo
    '200': [(0,0,0,1)],
    '020': [(1,0,1,1)],
    '002': [(2,0,2,1)],
    # Caso Tercero
    '211': [(0,1,1,0),(2,0,0,0)],
    '121': [(1,1,2,0),(0,0,1,0)],
    '112': [(2,1,0,0),(1,0,2,0)],
    # Caso Cuarto
    '220': [(0,1,1,0),(1,1,0,0)],
    '022': [(1,1,2,0),(2,1,1,0)],
    '202': [(2,1,0,0),(0,1,2,0)],
    # Caso Quinto
    '222': [(0,1,1,0),(1,1,2,0),(2,1,0,0)],
}

def triangle_sphere(t, s):
    """
    Return the intersection points between segment and sphere.

    Primer Caso Nulo. Sin intersecciones. Vertices todos adentro.

    >>> s = sphere([0, 0, 0], 10.0)
    >>> t = triangle([0, 0, 0],[1, 0, 0],[0, 1, 0])
    >>> triangle_sphere(t, s)
    []

    Segundo Caso Nulo. Sin intersecciones. Vertices todos afuera.

    >>> s = sphere([0, 0, 0], 10.0)
    >>> t = triangle([-20, -11, 0],[20, -11, 0],[0, 20, 0])
    >>> triangle_sphere(t, s)
    [(6.2831853071795862, array([ 0.,  0.,  0.]), 10.0, None, None)]

    Tercer Caso Nulo. Sin intersecciones. Vertices todos afuera.

    >>> s = sphere([0, 0, 0], 10.0)
    >>> t = triangle([20, 0, 0],[30, 0, 0],[25, 10, 0])
    >>> triangle_sphere(t, s)
    []

    Primer caso. Dos intersecciones. Dos vertices adentro y uno afuera.

    >>> t = triangle([0, 0, 0],[11, 0, 0],[0, 1, 0])
    >>> triangle_sphere(t, s)
    [(0.009094794229760993, array([ 0.,  0.,  0.]), 10.0, array([ 10.,   0.,   0.]), array([ 9.99958643,  0.09094669,  0.        ]))]

    Segundo caso. Dos intersecciones. Un vertice adentro y dos afuera.

    >>> t = triangle([0, 0, 0],[11, 0, 0],[11, 1, 0])
    >>> triangle_sphere(t, s)
    [(0.090659887200743763, array([ 0.,  0.,  0.]), 10.0, array([ 9.95893206,  0.90535746,  0.        ]), array([ 10.,   0.,   0.]))]

    Tercer caso. Cuatro intersecciones. Un vertice adentro y dos afuera.

    >>> t = triangle([-11, 0, 0],[11, 0, 0],[0, 1, 0])
    >>> triangle_sphere(t, s)
    [(0.009094794229760993, array([ 0.,  0.,  0.]), 10.0, array([ 10.,   0.,   0.]), array([ 9.99958643,  0.09094669,  0.        ])), (0.009094794229760993, array([ 0.,  0.,  0.]), 10.0, array([-9.99958643,  0.09094669,  0.        ]), array([-10.,   0.,   0.]))]

    Cuarto caso. Cuatro intersecciones. Tres vertices afuera.

    >>> t = triangle([-15, -11, 0],[15, -11, 0],[0, 11, 0])
    >>> triangle_sphere(t, s)
    [(0.13980786035649304, array([ 0.,  0.,  0.]), 10.0, array([ 0.69847012,  9.97557715,  0.        ]), array([-0.69847012,  9.97557715,  0.        ])), (2.5334834342706425, array([ 0.,  0.,  0.]), 10.0, array([-9.54130421, -2.99391283,  0.        ]), array([ 9.54130421, -2.99391283,  0.        ]))]

    Quinto caso. Seis intersecciones. Tres vertices afuera.

    >>> t = triangle([-11, 0, 0],[11, 0, 0],[0, 11, 0])
    >>> triangle_sphere(t, s)
    [(0.10578747987904831, array([ 0.,  0.,  0.]), 10.0, array([ 10.,   0.,   0.]), array([ 9.94409721,  1.05590279,  0.        ])), (0.21157495975809609, array([ 0.,  0.,  0.]), 10.0, array([ 1.05590279,  9.94409721,  0.        ]), array([-1.05590279,  9.94409721,  0.        ])), (0.10578747987904831, array([ 0.,  0.,  0.]), 10.0, array([-9.94409721,  1.05590279,  0.        ]), array([-10.,   0.,   0.]))]
    """
    olderr = seterr(all='ignore')

    # Calcula si el plano del triangulo cruza la esfera, sino devuelve una lista vacia.
    C = plane_sphere(t.plane(),s)
    if not C: return []

    # Calcula el mapa de los vertices del triangulo que están dentro de la esfera.
    # Distancias negativas están dentro de la esfera.
    ds = bitarray(map(lambda p: s.dist(p) <= 0, t.p))

    # Calcula los puntos de intersección de cada uno de los segmentos del triángulo
    # a la esfera.
    segs = [ segment(t.p[i], t.p[(i+1)%3]) for i in range(3) ]
    ssis = [ segment_sphere(seg, s) for seg in segs ]
    ssif = [ p for pts in ssis for p in pts ]

    # Calcula el número de intersecciones que ocurrieron.
    c = len(ssif)

    # Split small circle in origin and radious.
    o,r = C

    seterr(**olderr)

    # Devuelve los arcos según la disposición de los vértives y la
    # cantidad de cortes.

    if c == 0 and ds in [bitarray('111')]:
        return []

    elif c == 0 and ds in [bitarray('000')] and o not in t:
        return []

    elif c == 0 and ds in [bitarray('000')] and o in t:
        return [(2*pi,)+C+(None,None)]
       
    elif c in [2, 4, 6]:
        T = _tim_['%i%i%i' % tuple(map(len, ssis))]
        return [ _arc_(ssis[aa][bb], ssis[cc][dd], C[0], C[1]) for aa,bb,cc,dd in T ]

    raise RuntimeError

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4

