# -*- coding: ISO-8859-1 -*-
# $Id: storage.py 64 2009-05-19 17:29:15Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
""" Plugins de funciones de conjunto para la clase Storage.
"""

def sqlf_2_eqitems(A, B):
    """
    Devuelve el número de items equivalentes en dos conjuntos.

    >>> sqlf_2_eqitems("1,2,3", "2,3")
    2
    >>> sqlf_2_eqitems("A,AB,C", "A,C")
    2
    >>> sqlf_2_eqitems("3,AB,C", "A,C")
    1
    """
    A = A.split(',')
    B = B.split(',')
    if len(A) > len(B):
        A, B = B, A
    c = 0
    for i in A:
        if i in B:
            c = c + 1
    return c

def sqlf_1_sort_group(A):
    """
    Ordena una lista de elementos separados por comas

    >>> sqlf_1_sort_group("B,C,D,A")
    'A,B,C,D'
    >>> sqlf_1_sort_group("B,D,C,A")
    'A,B,C,D'
    """
    return ','.join(sorted(A.split(',')))

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

