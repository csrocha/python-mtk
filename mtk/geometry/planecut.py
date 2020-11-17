from numpy import *
from numpy import linalg

def cutpoint(plane, line):
    """
    Return the intersection point between line and plane.

    >>> face = (array((0.,0.,0.)), array((1.,0.,0.)), array((0.,1.,0.)))
    >>> plane = (face[0], cross(face[0] - face[1], face[0] - face[2]))
    >>> line = (array((-.5,-.5,-.5)), array((.5,.5,.5)))
    >>> cutpoint(plane, line)
    array([ 0.,  0.,  0.])
    >>> line = (array((1.,1.,-.5)), array((1.,1.,.5)))
    >>> cutpoint(plane, line)
    array([ 1.,  1.,  0.])
    """
    oplane, nplane = plane
    oline, vline = (line[0], line[1] - line[0])
    f = dot(vline, nplane)
    if f != 0:
        d = dot((oplane - oline), nplane) / f
        p = d * vline + oline
        n = linalg.norm(line[0] - line[1])
        if allclose(n, linalg.norm(p - line[0]) + linalg.norm(p - line[1])):
            return p
    return None

def cutline(plane, face):
    """
    Return start and end points for the intersection line between face and plane.

    >>> face = (array((0.,0.,0.)), array((1.,0.,0.)), array((0.,1.,0.)))
    >>> plane = (face[0], cross(face[0] - face[1], face[0] - face[2]))
    >>> face = (array((.5,.5,-.5)), array((1.,1.,.5)), array((0.,0.,0.5)))
    >>> cutline(plane, face)
    (array([ 0.75,  0.75,  0.  ]), array([ 0.25,  0.25,  0.  ]))
    """
    face = tuple(face)
    faceO = face
    faceD = face[1:] + ( face[0], )
    p = [ cutpoint(plane, L) for L in zip(faceO, faceD) ]
    p = [ x for x in p if x != None ]
    return tuple(p)

def cutfaces(plane, faces, vertexs):
    """
    Return a list of closed polygons generated from the cut of the faces
    in a plane.

    >>> vertexs = [ 1, 1, 1, \
                    1, 1,-1, \
                    1,-1, 1, \
                    1,-1,-1, \
                   -1, 1, 1, \
                   -1, 1,-1, \
                   -1,-1, 1, \
                   -1,-1,-1, ]
    >>> vertexs = array(vertexs).reshape((8,3))
    >>> faces = [ (0,1,5), (0,5,4), \
                  (4,5,7), (4,7,6), \
                  (6,7,3), (6,3,2), \
                  (3,0,2), (3,1,0) ]
    >>> plane = (array((0.,0.,0.)), array((0.,0.,1.)))
    >>> K = cutfaces(plane, faces, vertexs)
    >>> faces = [ (0,1,5,4), \
                  (4,5,7,6), \
                  (6,7,3,2), \
                  (2,3,1,0) ]
    >>> K = cutfaces(plane, faces, vertexs)
    >>> faces = [ (0,1,5,4), \
                  (4,5,7,6), \
                  (6,7,3,2), \
                  (2,3,1,0), \
                  (0,2,4,6), \
                  (1,3,7,5), \
                ]
    """
    lines = [ cutline(plane, map(lambda i: vertexs[i], face))
             for face in faces ]
    lines = filter(lambda l: len(l) == 2, lines)
    return lines

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:
