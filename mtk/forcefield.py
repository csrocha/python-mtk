#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
# $Id:$
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""
.. module:: forcefield
   :platform: Unix, Windows
   :synopsis: This module help to resolve values of forcefield of a molecule.

.. moduleauthor:: Cristian S. Rocha <crocha@dc.uba.ar>

"""

import os.path
from mtk.storage import Storage
from numpy import zeros

def editDistance(A, B):
    """
    Return the edit distance between two sequences.

    >>> A = [u'C', u'CA', u'HA', u'HN', u'N', u'O']
    >>> B = [u'C', u'CA', u'CD', u'HA', u'HD', u'HN', u'N', u'O']
    >>> editDistance(A, B) == 2
    True
    """
    if len(A) == 0 or len(B) == 0:
        return max(len(A), len(B))

    def D(a,b):
        if a == b: return 0
        else: return 1

    M = zeros((len(A), len(B)))
    for i in range(len(A)):
        M[i,0] = D(A[i], B[0])
    for i in range(len(B)):
        M[0,i] = D(B[i], A[0])
    for i in range(1, len(A)):
        for j in range(1, len(B)):
            M[i,j] = min( [
                M[i-1,j-1] + D(A[i],B[j]),
                M[i,j-1]   + 1,
                M[i-1,j]   + 1,
                ])
    return M[-1,-1]

def relativeEditDistance(A, B):
    """
    Return the edit distance of sequences relatives to their lenght.

    >>> A = [u'C', u'CA', u'HA', u'HN', u'N', u'O']
    >>> B = [u'C', u'CA', u'CD', u'HA', u'HD', u'HN', u'N', u'O']
    >>> relativeEditDistance(A, B) == 2.0/14.0
    True
    >>> A = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
    >>> B = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE']
    >>> relativeEditDistance(A, B) == 5.0/len(A+B)
    True
    """
    return editDistance(A, B) / (len(A) + len(B))

def aminoacids(S):
    """
    Return list of aminoacids in the database

    >>> S = Storage()
    >>> A = aminoacids(S)
    >>> len(A) > 0
    True
    """
    d = S.do("SELECT * FROM aminoacids")

    A = {}
    for resname, triname, letter in d:
        if not triname in A:
            A[triname] = []
        A[triname].append(resname)
    return A


def residuesTable(S):
    """
    Return residues on the database

    >>> S = Storage()
    >>> T = residuesTable(S)
    >>> len(T) > 0
    True
    """
    forcefield = os.path.join(S.path, 'test', 'amber99.prm')
    S.loadprm(forcefield)
    d = S.do("SELECT * FROM ffbiotype")
    table = {}
    for i in d:
        id, atomname, resname, atomid = i
        if not resname in table:
            table[resname] = []
        table[resname].append(atomname)
    for i in table:
        table[i].sort()
    return table

def relateAminoacids(S, pdbid):
    """
    Return the most related aminoacid of pdbid
    >>> from mtk.tests import opentestfile
    >>> S = Storage()
    >>> id = S.loadpdb(opentestfile("1A2K_l_b.pdb"))
    >>> R = relateAminoacids(S, id)
    >>> len(R) > 0
    True
    """
    T = residuesTable(S)
    A = aminoacids(S)
    d = S.do("SELECT * FROM molecule WHERE pdbid = ?", pdbid)
    # Generate list of atoms for each residue from the pdb
    residues = {}
    for i in d:
        pdbid, record, serial, name, altloc, resname, domain, resid = i[:8]
        if not resid in residues:
            residues[resid] = (resname, [], None)
        residues[resid][1].append(name)
    for i in residues:
        residues[i][1].sort()
    # Determine with residue forcefield is related to each residue from the pdb
    relation = {}
    for i in residues:
        M = []
        for j in A[residues[i][0]]:
            M.append((j, relativeEditDistance(residues[i][1], T[j])))
        relation[i] = min( M, key=lambda (a,b): b)[0]

    return relation

def moleculeForceField(S, pdbid):
    """
    Generate forcefield information for the molecule

    TODO: Not Implemented. I must inferer the values from some molecules. A
    complex theme.

    >>> from mtk.tests import opentestfile
    >>> S = Storage()
    >>> id = S.loadpdb(opentestfile("1A2K_l_b.pdb"))
    >>> id = S.loadpdb(opentestfile("1ACB_r_b.pdb"))
    >>> R = moleculeForceField(S, id)

    Test not terminated
    """
    R = relateAminoacids(S, pdbid)
    query = ["""
        CREATE TEMPORARY TABLE IF NOT EXISTS Tresidues
        (resSerial INTEGER, TaaName TEXT)
        """]
    #S.do(

    return R

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

