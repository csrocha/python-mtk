#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
# $Id: drawmol.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

def geometry_complementary(r, v, core=1+10j, surf=10+1j, solv=0, r_core=2.0, r_solv=3.0):
    """
    Function to complete complemetary geometry.

    >>> geometry_complementary(0, 1+10j) == 1+10j
    True
    """
    if v!=core and r<=r_core:
        return core
    if v==solv and r<=r_solv:
        return surf
    return v

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

