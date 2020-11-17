from numpy import array, argmax, argmin, float, mgrid, dot, allclose, arange, round
from math import copysign
from numpy.linalg import norm

intround = lambda x: int(round(x))

def line(a, b, debug=False):
    """
    >>> list(line((0,0), (3,3)))
    [array([ 0.,  0.]), array([ 1.,  1.]), array([ 2.,  2.]), array([ 3.,  3.])]
    >>> list(line((0,0), (-3,-3)))
    [array([ 0.,  0.]), array([-1., -1.]), array([-2., -2.]), array([-3., -3.])]
    >>> list(line((0,0), (3,-3)))
    [array([ 0.,  0.]), array([ 1., -1.]), array([ 2., -2.]), array([ 3., -3.])]
    >>> list(line((0,0), (-3,3)))
    [array([ 0.,  0.]), array([-1.,  1.]), array([-2.,  2.]), array([-3.,  3.])]
    >>> list(line((3,3), (0,0)))
    [array([ 3.,  3.]), array([ 2.,  2.]), array([ 1.,  1.]), array([ 0.,  0.])]
    >>> list(line((1,1), (-1,-1)))
    [array([ 1.,  1.]), array([ 0.,  0.]), array([-1., -1.])]
    >>> list(line((1,1), (1,-1)))
    [array([ 1.,  1.]), array([ 1.,  0.]), array([ 1., -1.])]
    >>> list(line((1,1), (-1,1)))
    [array([ 1.,  1.]), array([ 0.,  1.]), array([-1.,  1.])]
    >>> list(line((0,0), (1,2)))
    [array([ 0.,  0.]), array([ 0.5,  1. ]), array([ 1.,  2.])]
    >>> list(line((0, 0, 0), (0, 0, 0)))
    [array([0, 0, 0])]
    """
    if debug: import pdb; pdb.set_trace()
    a, b = array(map(intround,a)), array(map(intround,b))
    delta = b - a
    if all(delta == 0):
        yield a
    else:
        axis = argmax(abs(delta))
        factors = delta / float(delta[axis])
        xtop = delta[axis]
        dire = int(copysign(1, xtop))
        for x in xrange(0, xtop + dire, dire):
            yield a + x*factors

def face(coords, debug=False):
    assert(len(coords) <= 4)
    if len(coords) == 3:
        for x in face3(*coords):
            yield x
        coords = [ coords[0], coords[2], coords[1] ]
        for x in face3(*coords):
            yield x
    else:
        coords.append(coords[0])
        for x in face3b(*coords[:3]):
            yield x
        for x in face3b(*coords[2:]):
            yield x


def face3(a, b, c, debug=False):
    """
    > >> list(face3((0,0,0), (2,0,0), (0,2,0)))
    [array([0, 0, 0]), array([ 1.,  0.,  0.]), array([ 0.,  1.,  0.]), array([ 2.,  0.,  0.]), array([ 1.,  1.,  0.]), array([ 0.,  2.,  0.])]
    
    >>> P = [array([0., 0., 0.]), array([ 1.,  2.,  0.]), array([ 1.,  1.,  0.]), array([ 1.,  0.,  0.]), array([ 1., -1.,  0.]), array([ 1., -2.,  0.])]
    >>> L = list(face3((0,0,0), (1,2,0), (1,-2,0)))
    
    > >> [ any([(o == p).all() for p in P]) for o in L ] == [True, False, True, True, True, True, True, True, True]
    True
    
    >>> L = list(face3((1,2,0), (0,0,0), (1,-2,0)))
    >>> all([ any([(o == p).all() for p in P]) for o in L ])
    True
    >>> L = list(face3((1,2,0), (1,-2,0), (0,0,0)))
    >>> all([ any([(o == p).all() for p in P]) for o in L ])
    True
    >>> L = list(face3((1,-2,0), (1,2,0), (0,0,0)))
    >>> all([ any([(o == p).all() for p in P]) for o in L ])
    True
    >>> L1 = list(face3((0,0,0), (0,10,0), (0,0,10)))
    >>> L2 = list(face3((0,0,0), (0,0,10), (0,10,0)))
    >>> all([ any([(o == p).all() for p in L1]) for o in L2 ])
    True
    """
    a, b, c = map(array, (a,b,c))
    ab = norm(a - b)
    ac = norm(a - c)
    bc = norm(b - c)
    if ab > ac and ab > bc:
        A = line(a, b)
        partner = [ line(a, c), line(c, b) ]
    elif ac > ab and ac > bc:
        A = line(a, c)
        partner = [ line(a, b), line(b, c) ]
    elif bc > ab and bc > ac:
        A = line(b, c)
        partner = [ line(b, a), line(a, c) ]
    else:
        print a, b, c
        print ab, bc, ac
    for ax in A:
        try:
            bx = partner[0].next()
        except:
            partner.pop(0)
            bx = partner[0].next()
        for x in line(ax, bx):
            yield x

def face3b(a, b, c, debug=False):
    """
    >>> list(face3((0,0,0), (1,2,0), (1,-2,0)))
    >>> list(face3b((0,0,0), (1,2,0), (1,-2,0)))
    []
    """
    F = array([a, b, c], dtype=float)

    done = []
    o = dot(array([1./3,1./3,1./3]), F)
    d = 1/128. # 1/(max([norm(a - o), norm(b - o), norm(c - o)]))
    l = 1. + d

    for x in arange(0, l, d):
        for y in arange(0, l-x, d):
            z = 1 - (x+y)
            if z >= 0:
                t = array(round(dot(array([x, y, z]), F)), dtype=int)
                yield t

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

