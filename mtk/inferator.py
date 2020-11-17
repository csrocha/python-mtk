# -*- coding: ISO-8859-1 -*-
#
# $Id:$
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
""" Módulo para inferir átomos que no pueden ser reconocidos de forma 
    directa por la base de datos.
"""

def serialize_atoms(atoms):
    """
    Generate a list of sequences with residues with atoms

    >>> A = [ ('N', 'PHE', 1, 1.0),('CA', 'PHE', 1, 1.0),('C', 'PHE', 1, 1.0),('O', 'PHE', 1, 1.0) ]
    >>> serialize_atoms(A) == {1: ('PHE', [('N', 1.0), ('CA', 1.0), ('C', 1.0), ('O', 1.0)])}
    True
    >>> A += [ ('N', 'PHE', 2, 1.0),('CA', 'PHE', 2, 1.0),('C', 'PHE', 2, 1.0),('O', 'PHE', 2, 1.0) ]
    >>> serialize_atoms(A) == {1: ('PHE', [('N', 1.0), ('CA', 1.0), ('C', 1.0), ('O', 1.0)]), \
    2: ('PHE', [('N', 1.0), ('CA', 1.0), ('C', 1.0), ('O', 1.0)])}
    True

    """
    D = {}
    for name, resname, resid, property in atoms:
        if not resid in D:
            D[resid] = (resname, [])
        D[resid][1].append( (name, property) )
    return D

def build_tree(serial):
    """
    Build a tree of atoms

    TODO: No implemented. It's must help to search the similiar residue
    estructure in the database.

    >>> A = [ ('N', 'PHE', 1, 1.0),('CA', 'PHE', 1, 1.0),('C', 'PHE', 1, 1.0),('O', 'PHE', 1, 1.0) ]
    >>> build_tree(serialize_atoms(A))
    Traceback (most recent call last):
        ...
    NotImplementedError

    """
    raise NotImplementedError
    T = {}
    for i in serial:
        resname, atoms = serial[i]
        if not resname in T:
            T[resname] = { }

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

