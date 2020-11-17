#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
# $Id: drawmol.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
#
#  Code based in Modelling Protein Docking using Shape
#                Complementarity, Electrostatics and
#                Biochemical Information.
#                J Mol Biol, Vol. 272, No. 1. (1997), pp. 106-20.
#                by H. A. Gabb, R. M. Jackson, M. J. Sternberg

from mtk.geometry.neighbours import distance, distance_to
from numpy import array, sum, where, logical_and, logical_not, piecewise

def epsilon(r):
    if r <= 6.0:
        return 4
    elif r <= 8.0:
        return 38*r - 224
    else:
        return 80

def charge_1(dists, charges):
    """
    Calcule the charge in an coordinate influed by atoms in dists distance with charges charge.
    Inputs dists and charges are array 1D of same range of floats.

    >>> from numpy import allclose
    >>> dists = array(range(5)) + 0.1
    >>> charges = array(range(5))
    >>> vdws = array([ 1 ] * 5)
    >>> c = charge_1(dists, charges)
    >>> allclose(c, 0.95120) 
    True

    >>> dists = array(range(10)) + 0.1
    >>> charges = array(range(10))
    >>> vdws = array([ 1 ] * 10)
    >>> c = charge_1(dists, charges)
    >>> allclose(c, 1.36864) 
    True
    
    """
    charge = charges / ( map(epsilon, dists) * dists )
    return sum(charge)

def charge_2(dists, charges):
    """
    Calcule the charge in an coordinate influed by atoms in dists distance with charges charge.
    Inputs dists and charges are array 1D of same range of floats.

    >>> from numpy import allclose
    >>> dists = array(range(5)) + 0.1
    >>> charges = array(range(5))
    >>> vdws = array([ 1 ] * 5)
    >>> c = charge_2(dists, charges)
    >>> allclose(c, 0.95120)
    True

    >>> dists = array(range(10)) + 0.1
    >>> charges = array(range(10))
    >>> vdws = array([ 1 ] * 10)
    >>> c = charge_2(dists, charges)
    >>> allclose(c, 1.36864)
    True

    """
    d6 = dists <= 6.0
    d8 = dists <= 8.0
    d6_8 = logical_and(logical_not(d6), d8)
    epsilons = (d6*4.0) + \
            d6_8*(38.0*dists-224.0) + \
            logical_not(d8)*80.0
    charge = (charges / ( epsilons * dists ))
    return sum(charge)

class electrostatic_vol:
    def __init__(self, atoms):
        self.atoms = atoms

    def __call__(self, x, v):
        rs = distance_to(x, self.atoms[:,0:3]) # TODO: Reducir tiempo!!!
        if all(rs > self.atoms[:,-2]):
            #return charge_2(rs, self.atoms[:,-1])
            return charge_1(rs, self.atoms[:,-1])
        else:
            return 0
        return r

def calc(V, atoms):
    V.forallnodes(electrostatic_vol(atoms))

def test_suite():
    import doctest
    return doctest.DocTestSuite()

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
