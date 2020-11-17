#!/usr/bin/python

from numpy import *

F = array([[ 0, 0, 0 ], [ 3, 0, 0 ], [ 0, 10, 10 ]], dtype=float).T
F = array([[ 0, 0, 0 ], [ 3, 0, 0 ], [ 0, 10, 0 ]], dtype=float)

R = arange(0,1.1,.1)

X, Y, Z = mgrid[0:1.1:0.1, 0:1.1:0.1, 0:1.1:0.1]

for x,y,z in zip(X.flat,Y.flat,Z.flat):
    print x, y, z
    #dot(F, array([x, y, z]))
                #print array(dot(F, array([x, y, z])),dtype=int)

