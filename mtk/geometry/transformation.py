#!/usr/bin/python
# -*- coding: latin-1 -*-

from numpy import sin, cos, radians, \
        array, zeros, identity, \
        dot, cross, arange, \
        pi, allclose
from math import acos
from numpy.linalg import norm
from sys import stderr

def normalize(v):
    """
    Normalize a vector

    >>> a = array([1,0,0])
    >>> allclose(normalize(a), a)
    True
    >>> a = array([1,1,1])
    >>> allclose(normalize(a), [0.57735]*3)
    True
    """
    return v / norm(v)

def angleOf(a, b):
    """
    Calcule the angle between a and b

    >>> a = array([1,0,0])
    >>> b = array([0,1,0])
    >>> allclose(angleOf(a, b), radians(90))
    True
    """
    a, b = normalize(a), normalize(b)
    return acos(dot(a, b))

def applyTv(R, v):
    """
    Apply a transformation R to a afine vector v
    """
    try:
        v = dot(R, v)
    except ValueError:
        print >> stderr, R.shape, v.shape
        raise
    return v[:3] / v[3]

def applyTm(R, M):
    """
    Apply a transformation R to a matrix M with vectors
    """
    try:
        v = dot(R, M)
    except ValueError:
        print >> stderr, R.shape, M.shape
        raise
    return v[:3,:] / v[3,:]

def transform(T, x):
    """
    Apply a transformation T to the vector x
    """
    v = dot(T, array( list(x)[:3] + [1.0] ))
    return v/v[3]

def translation(v):
    """
    generate a matrix for translation in v direction.
    """
    T = identity(4)
    T[0:3,3] = v
    return T

def rotation(A, B, G, o = (0,0,0)):
    """
    generate a matrix for rotation in angles A,B,G in axis X,Y,Z with o origin.
    """
    T = identity(4)
    cosA, cosB, cosG = map(cos, (A, B, G))
    sinA, sinB, sinG = map(sin, (A, B, G))
    T = array( [
       [  cosB*cosG,  sinA*sinB*cosG + cosA*sinG, -cosA*sinB*cosG+sinA*sinG, o[0] ],
       [ -cosB*sinG, -sinA*sinB*sinG + cosA*cosG,  cosA*sinB*sinG+sinA*cosG, o[1] ],
       [       sinB,                  -sinA*cosB,                 cosA*cosB, o[2] ],
       [ 0, 0, 0, 1 ]
    ] )
    return T

def axisRotation(w, phi):
    """
    Genera una matriz rotacion a partir de una rotación phi sobre el eje n.
    http://mathworld.wolfram.com/RodriguesRotationFormula.html

    >>> w = array([1,0,0])
    >>> axisRotation(w, 0)
    array([[ 1.,  0.,  0.,  0.],
           [ 0.,  1.,  0.,  0.],
           [ 0.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  1.]])
    >>> D = normalize(array([1,1,1]))
    >>> R = normalize(array([1,0,0]))
    >>> a = -acos(dot(D, R))
    >>> d = normalize(cross(D, R))
    >>> AX = axisRotation(d, a)
    >>> allclose( AX[0:3,0], D)
    True
    >>> allclose( applyTv(AX, [1,0,0,1]), D )
    True
    >>> allclose( applyTv(AX, array([[1,0,0,1]]).T), D )
    True
    """
    W = array([[    0,-w[2], w[1], 0],
               [ w[2],    0,-w[0], 0], 
               [-w[1], w[0],    0, 0],
               [    0,    0,    0, 0]], dtype=float)
    R = identity(4) + W*sin(phi) + dot(W, W)*(1-cos(phi))

    return R

def isEqual(A, B, tolerance=0.2, v=(1,0,0)):
    return dot(A, v) == dot(B, v)

def genRotationsByAngle(step_alpha, step_beta, step_gamma, o=(0,0,0)):
    c = 0
    R = []
    for alpha in arange(0, radians(180), step_alpha):
        for beta in arange(0, radians(360), step_beta):
            for gamma in arange(0, radians(360), step_gamma):
                R.append(rotation(alpha, beta, gamma, o))
    return R

def genRotationsByDirections(directions, step_gamma, reference=array([1.,0.,0.])):
    """
    Generate rotations based in normalized directions and a step angle.

    Checking directions of a cube.

    >>> d = [ [1., 1., 1.], [1., 1., -1.], [1., -1., 1.], [1., -1., -1.], \
              [-1., 1., 1.], [-1., 1., 1.], [-1., -1., 1.], [-1., -1., -1.], ]
    >>> R = genRotationsByDirections(d, radians(360))
    >>> v = array([1., 0., 0., 1.])
    >>> allclose(map(abs, array([ list(applyTv(r, v).flat) for r in R ]).flat), 0.577350269)
    True

    Generating from a list of points sphere

    >>> from mtk import resource
    >>> import csv
    >>> R = csv.reader(open(resource("csv", "UnitSphere_500.csv")))
    >>> directions = array([ map(float, line) for line in R ])
    >>> allclose(directions[0,:], array([ 0.80635011, -0.58584762, -0.08112963]))
    True
    >>> R = genRotationsByDirections(directions, radians(15))
    >>> all([ allclose(applyTv(R[0], [1,0,0,1]), directions[0,:] ) for i in xrange(len(R)) ])
    True
    >>> len(R)
    6048
    >>> v = array([1., 0., 0., 1.])
    >>> from numpy.linalg import norm
    >>> all([ allclose(norm(applyTv(r, v)), 1.) for r in R ])
    True
    """
    from numpy.linalg import norm
    directions = array(directions)
    reference = reference / norm(reference)
    R = []
    for i in xrange(directions.shape[0]):
        direction = directions[i,:] / norm(directions[i,:])
        if not allclose(direction, reference):
            angle = angleOf(direction, reference)
            rotation_axis = -cross(direction, reference)
            rotation_axis = rotation_axis / norm(rotation_axis)
            B = axisRotation(rotation_axis, angle)
        else:
            B = identity(4,float)
        for phi in arange(0, radians(360), step_gamma):
            R.append(dot(axisRotation(direction, phi), B))
    return R

def objRotations(R, fd):
    """
    >>> from mtk import resource
    >>> import sys, csv
    >>> R = csv.reader(open(resource("csv","UnitSphere_500.csv")))
    >>> directions = array([ map(float, line) for line in R ])
    >>> R = genRotationsByDirections(directions, radians(15))
    >>> objRotations(R, open('/home/crocha/testing.obj', 'w'))
    """
    vf = "v %0.4f %0.4f %0.4f\n"
    ff = "f %s\n"
    V = array([
        [19.000000, -1.000000, 0.000000, 1.0],
        [19.000000, -1.000000, -1.000000, 1.0],
        [19.000000, 1.000000, 1.000000, 1.0],
        [19.000000, 1.000000, 0.000000, 1.0],
        [21.000000, -1.000000, 0.000000, 1.0],
        [21.000000, 1.000000,  0.000000, 1.0],
    ]).T
    F = [
        [3, 6, 4,],
        [1, 5, 2,],
        [1, 4, 6, 5,],
    ]
    b = 0
    for r in R:
        VR = applyTm(r, V)
        for i in xrange(0, VR.shape[1]): 
            fd.write(vf % tuple(VR[:,i].flat))
        for fr in F:
            fd.write(ff % ' '.join(map(str, map(lambda x: x + b, fr))))
        b = b + VR.shape[1]

def drawRotations(R = genRotationsByAngle(radians(45),radians(45),radians(45)), t="vectors"):
    """
    Dibuja las rotaciones en mlan

    > >> from mtk import resource
    > >> from enthought.mayavi import mlab
    > >> import csv
    > >> R = csv.reader(open(resource("csv","UnitSphere_500.csv")))
    > >> directions = array([ map(float, line) for line in R ])
    > >> R = genRotationsByDirections(directions, radians(15))
    > >> r = drawRotations(R=R)
    > >> len(R)
    6048
    > >> mlab.show()

    > >> from enthought.mayavi import mlab
    > >> R = genRotationsByAngle(radians(15), radians(15), radians(15))
    > >> r = drawRotations(R=R)
    > >> mlab.show()
    """
    from sys import stderr
    from numpy import vstack, hstack

    if t=="points":
        from enthought.mayavi.mlab import points3d
        # Draw points
        V = array([10,0,0,1])
        points3d(0,0,0,color=(1.0,0.5,0.5))
        P = []
        for r in R:
            P.append(applyTv(r, V))
        P = hstack(P)
        points3d(P[0,:], P[1,:], P[2,:])
    if t=="vectors":
        from enthought.tvtk.tools import visual

        V = array([[10,0,0,1],
                   [11,0,0,1],
                   [11,.5,0,1],
                   [11,.5,.5,1]]).T

        for r in R:
            P = applyTm(r, V)
            visual.Curve(points=P[0:3,:])
        pass

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

