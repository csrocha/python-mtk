# -*- coding: ISO-8859-1 -*-
#!/usr/bin/python

import csv
from numpy import load, array, radians
from mtk.io.pqr_ff import reader
from mtk import resource
from mtk.energy import create_geometry_vol
from mtk.geometry.transformation import applyTm, genRotationsByDirections
from mtk.geometry.vol import BasicVolume
from mtk.dock import load as loaddock
from structures import read, read_rotations

R = read_rotations()
results = loaddock("dock_results-0000.csv")
Mr, Ml, resolution = read()

print len(Mr), len(Ml)

parms = { 'resolution': resolution, 'solv': 0, 'core': 1+10j, 'surf': 10+1j,
         'r_solv': resolution*2 }

from mtk.geometry import size as geosize

l_minimal = array([ geosize(Ml).max() ] *3)
r_minimal = array([ geosize(Mr).max() ] *3)
size = l_minimal*2 + r_minimal

T = R[0]

SASr = create_geometry_vol(Mr[:,:3], Mr[:,4], size=size, **parms)
SASl = create_geometry_vol(applyTm(T, Ml[:,:4].T).T[:,:3], Ml[:,4], size=size, **parms)
#S = BasicVolume(SASr.min, SASr.delta, load(open('dock_fft.npy')))
S = None

from numpy import all, abs, load

from mtk.view.dock import Viewer
DockViewer = Viewer(results, SASr, SASl, S=S)

#import vtk
#
#AtomPosition = vtk.vtkPoints()
#AtomType = vtk.vtkFloatArray()
#for row in Mr:
#    AtomPosition.InsertNextPoint(row[:3])
#    AtomType.InsertNextValue(row[4])
#
#Molecule=vtk.vtkPolyData()
#Molecule.SetPoints(AtomPosition)
#Molecule.GetPointData().SetScalars(AtomType)
#
DockViewer.run()

