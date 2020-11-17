# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
#
# Source: http://www.cs.mun.ca/~rod/2500/notes/numpy-arrays/numpy-arrays.html
#
# line segment intersection using vectors
# see Computer Graphics by F.S. Hill
#
from numpy import *

def perp( a ) :
    b = empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return 
def seg_intersect(a1,a2, b1,b2) :
    """
    segment intersection check

    >>> p1 = array( [0.0, 0.0] )
    >>> p2 = array( [1.0, 0.0] )
    >>> p3 = array( [4.0, -5.0] )
    >>> p4 = array( [4.0, 2.0] )

    >>> seg_intersect( p1,p2, p3,p4)
    array([ 4.,  0.])

    >>> p1 = array( [2.0, 2.0] )
    >>> p2 = array( [4.0, 3.0] )
    >>> p3 = array( [6.0, 0.0] )
    >>> p4 = array( [6.0, 3.0] )

    >>> seg_intersect( p1,p2, p3,p4)
    array([ 6.,  4.])

    >>> p1 = array( [2.0, 2.0,1.0] )
    >>> p2 = array( [4.0, 3.0,1.0] )
    >>> p3 = array( [6.0, 0.0,1.0] )
    >>> p4 = array( [6.0, 3.0,1.0] )

    >>> seg_intersect( p1,p2, p3,p4)
    array([ 6.,  4.,  1.])


    """
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = dot( dap, db)
    num = dot( dap, dp )
    return (num / denom)*db + b1

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

