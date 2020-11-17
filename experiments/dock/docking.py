# -*- coding: ISO-8859-1 -*-
#!/usr/bin/python

from numpy import array, identity, radians
from mtk.dock import dock, save as savedock
from mtk.energy import create_geometry_vol
from mtk.geometry.transformation import applyTm
from time import time
from numpy import save
from structures import read, read_rotations
import os.path 

stime = time()

print "t %04.4f: Leyendo las rotaciones." % (time()-stime)

R = read_rotations()

print "t %04.4f: Leyendo las proteínas." % (time()-stime)

Mr, Ml, resolution = read()
parms = { 'resolution': resolution, 'solv': 0, 'core': 9j, 'surf': 1,
         'r_solv': resolution*2 }

print "t %04.4f: Creando geometría del receptor." % (time()-stime)

from mtk.geometry import size as geosize

l_minimal = array([ geosize(Ml).max() ] *3)
r_minimal = array([ geosize(Mr).max() ] *3)
size = l_minimal*2 + r_minimal

Va = create_geometry_vol(Mr[:,:3], Mr[:,4], size=size, **parms)

fVb = lambda T: create_geometry_vol(applyTm(T, Ml[:,:4].T)[:3], Ml[:,4],
                                    size=size, **parms)

if False:
    from mtk.view.vol import Viewer
    VolViewer = Viewer()
    VolViewer.append(Va)
    VolViewer.append(fVb(R[:1]))
    VolViewer.run()

print "t %04.4f: Filtrando resultados existentes." % (time()-stime)
start = 0
for i in xrange(len(R)):
    if os.path.exists("dock_results-%04i.csv" % i):
        start = i+1
R = R[start:]

print "Empezando en %i" % start

print "t %04.4f: Dockeando." % (time()-stime)

def progress_function(state, data):
    if state in ["receptor copy", "dock"]:
        print "Shape:", array(data.shape)
    if state in ["results selection"]:
        r, FFT = data
        save(open("dock_fft-%04i.npy" % (start + r), "w"), FFT)
        pass
    if state in ["docked"]:
        t, r, D = data
        savedock(D, "dock_results-%04i.csv" % (start + r))
        print "t %04.4f: Stored %i results." % (time()-stime, len(D[0]))
    print "t %04.4f: Docking in state %s." % (time()-stime, state)

D = dock(Va, fVb, R, filter_function=lambda v: v > 400, progress_function=progress_function)

print "t %04.4f: Finish." % (time()-stime)

