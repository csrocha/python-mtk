# -*- coding: ISO-8859-1 -*-
# $Id: neighbours.py 94 2009-06-17 11:52:27Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Neighbours groups calculator"""

from numpy import *
import numpy.linalg as LA

def distance(vectors):
    """Calcule distance Matrix of vectors
    
    >>> a = array([[0,0,0], [0,1,0], [1,0,0]])
    >>> distance(a)
    array([[ 0.        ,  1.        ,  1.        ],
           [ 1.        ,  0.        ,  1.41421356],
           [ 1.        ,  1.41421356,  0.        ]])

    >>> a = array([[0,0,0,1], [0,1,0,1], [1,0,0,1]])
    >>> distance(a)
    array([[ 0.        ,  1.        ,  1.        ],
           [ 1.        ,  0.        ,  1.41421356],
           [ 1.        ,  1.41421356,  0.        ]])
    """
    vectors = array(vectors)[:,:3]
    n,m = vectors.shape
    delta = zeros((n,n),'d')
    for d in xrange(m):
        data = vectors[:,d]
        delta += (data - data[:,newaxis])**2
    return sqrt(delta)

def distance_to(x, vectors):
    """Calculate distance from x to each item of vectors
    
    >>> a = [ [0,0,0], [0,1,0], [2,0,0], [1,1,1] ]
    >>> alltrue(distance_to([0,0,0], a) == array([ 0., 1.,  2., sqrt(3.)]))
    True
    """
    x = array(x)[:3]
    v = array(vectors)[:,:3]
    r = array(map(LA.norm, x - v))
    return r

def d(a, b):
    """
    Calcule distance between two vectors

    >>> a = array([0, 0, 0, 1])
    >>> b = array([1, 0, 1, 1])
    >>> d(a,b) == distance(array([a,b]))[0,1]
    True
    """
    d = a[:3] - b[:3]
    d2 = dot(d, d)
    return sqrt(sum(d2))

def nearest(D, d):
    """Groups data in neighbors of distance d

    >>> a = array([[0,0,0], [0,1,0], [1,0,0]])
    >>> D = distance(a)
    >>> nearest(D, 1.0)
    [[1, 2], [0], [0]]
    """
    o = [[] for i in xrange(D.shape[0])]
    for i in xrange(D.shape[0]):
        for j in xrange(i+1, D.shape[0]):
            if D[i][j] <= d:
                o[i].append(j)
                o[j].append(i)
    return o

def pairs(N):
    """Generate a list of pairs from a groups of nearest
    
    >>> a = array([[0,0,0], [0,1,0], [1,0,0]])
    >>> D = distance(a)
    >>> N = nearest(D, 1.0)
    >>> pairs(N)
    [(0, 1), (0, 2), (1, 0), (2, 0)]
    """
    P = []
    for i in xrange(len(N)): 
        P += [ (i,j) for j in N[i] ]
    return P

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

