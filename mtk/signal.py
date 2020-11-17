# -*- coding: ISO-8859-1 -*-
#
# $Id: vol.py 85 2009-06-07 23:46:26Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

from numpy import array, zeros, ndenumerate
from numpy.fft import fftn, ifftn, rfftn, irfftn

def naive_xcorr(A, B):
    """
    >>> A = [ [ [ 0+1j ] ] ]
    >>> B = [ [ [ 1+0j ] ] ]
    >>> naive_xcorr(A, B)
    array([[[ 0.-1.j]]])
    >>> A = [ [ [ 0+1j, 1 ] ] ]
    >>> B = [ [ [ 1+0j, 1 ] ] ]
    >>> naive_xcorr(A, B)
    array([[[ 1.-1.j,  1.-1.j]]])
    """
    A, B = array(A), array(B)
    assert(A.shape == B.shape)
    assert(A.dtype == B.dtype)
    shape = A.shape
    M = zeros(shape, dtype=A.dtype)
    for i, v in ndenumerate(M):
        for j, a in ndenumerate(A):
            M[i] += a.conjugate() * B[ tuple((array(i) + array(j)) % array(shape)) ]
    return M

def rxcorr(A, B):
    """Correlation between two matrices

    >>> from numpy import allclose

    >>> A = [ 0, 1 ]
    >>> B = [ 1, 0 ]
    >>> allclose(rxcorr(A, B), array([ 0*1 + 1*0, 1*1 + 0*0 ]))
    True
    >>> allclose(rxcorr(A, B), naive_xcorr(A,B))
    True
    >>> A = [[ 0, 1 ],[ 0, 1 ]]
    >>> B = [[ 1, 0 ],[ 0, 1 ]]
    >>> allclose(rxcorr(A, B), array( \
               [[ 0*1 + 1*0 + 0*0 + 1*1, 1*0 + 0*0 + 1*1 + 0*1 ], \
                [ 0*1 + 1*0 + 0*0 + 1*1, 1*1 + 0*0 + 0*1 + 1*0 ]]))
    True
    >>> allclose(rxcorr(A, B), naive_xcorr(A,B))
    True
    """
    Fa = rfftn(array(A))
    Fb = rfftn(array(B))
    Fab = Fa * Fb.conjugate()
    return irfftn(Fab)

def xcorr(A, B):
    """Correlation between two matrices

    >>> from numpy import allclose, ones

    >>> A = [ [ [ 0+1j ] ] ]
    >>> B = [ [ [ 1+0j ] ] ]
    >>> allclose(xcorr(A, B), array([[[ 0.-1.j]]]))
    True
    >>> allclose(xcorr(A, B), naive_xcorr(A, B))
    True
    >>> A = [ 0+1j, 1+0j ]
    >>> B = [ 1+0j, 0+1j ]
    >>> allclose(xcorr(A, B), naive_xcorr(A, B))
    True
    >>> A = [ 0, 1 ]
    >>> B = [ 1, 0 ]
    >>> allclose(xcorr(A, B), array([ 0*1 + 1*0, 1*1 + 0*0 ]))
    True
    >>> A = [[ 0, 1 ],[ 0, 1 ]]
    >>> B = [[ 1, 0 ],[ 0, 1 ]]
    >>> allclose(xcorr(A, B), array([[ 0*1 + 1*0 + 0*0 + 1*1, 1*0 + 0*0 + 1*1 + 0*1 ], \
                                      [ 0*1 + 1*0 + 0*0 + 1*1, 1*1 + 0*0 + 0*1 + 1*0 ]]))
    True
    >>> A = [[[ 0, 1, 1j ],[ 0, 1, 3j ]]]
    >>> B = [[[ 1, 0, 2j ],[ 0, 1, 4j ]]]
    >>> allclose(naive_xcorr(A, B), xcorr(A, B))
    True
    >>> A = [[[ 0, 10+1j, 1+100j, 10+1j, 0 ]]]
    >>> B = [[[ 0, 10+1j, 1+100j, 10+1j, 0 ]]]
    >>> A = array(A).conjugate()
    >>> allclose(naive_xcorr(A, B), xcorr(A, B))
    True

    >>> O = ones((4,4,4),dtype=complex)
    >>> O[1:3,1:3,1:3]=9j
    >>> A = zeros((12,12,12), dtype=complex)
    >>> A[4:8,4:8,4:8] = O
    >>> B = zeros((12,12,12), dtype=complex)
    >>> B[4:8,4:8,4:8] = O
    >>> C = xcorr(A.conjugate(), B)
    >>> allclose([C.real.min()],[ -592 ] )
    True
    >>> allclose([C.real.max()],[ 24 ])
    True

    """
    Fa = fftn(array(A)).conjugate()
    Fb = fftn(array(B))
    Fab = Fa * Fb
    return ifftn(Fab)

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

