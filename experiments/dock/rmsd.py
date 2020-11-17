# -*- coding: ISO-8859-1 -*-
#!/usr/bin/python

import csv
from numpy import load, array, radians, sum, identity
from mtk.io.pqr_ff import reader
from mtk import resource
from mtk.energy import create_geometry_vol
from mtk.geometry.transformation import applyTm, genRotationsByDirections
from mtk.geometry.vol import BasicVolume
from mtk.geometry.rms import rmsd
from mtk.dock import load as loaddock
from structures import read, read_rotations

def calcule_rmsd(T, x, Rcentre, Lcentre, LM):
    t =  (Lcentre - Rcentre) + x
    T[:3,3] = array(t)
    Cl = applyTm(T, LM[:,:4].T)
    v_rmsd = rmsd(LM[:,:3], Cl.T[:,:3])
    return v_rmsd

R = read_rotations()
results = loaddock("dock_results-0000.csv")
Mr, Ml, resolution = read()

Rcentre = sum(Mr[:,:3], axis=0)/Mr.shape[0]
Lcentre = sum(Ml[:,:3], axis=0)/Ml.shape[0]

for i in xrange(len(results[0])):
    r, x, s = results[0][i], results[1][i], results[2][i]
    T = R[r]
    RMSD = calcule_rmsd(T, x, Rcentre, Lcentre, Ml)
    if RMSD < 100:
        print s.real, RMSD

