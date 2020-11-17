# -*- coding: ISO-8859-1 -*-
from numpy import array, radians, identity
from mtk.io.pqr_ff import reader
from mtk import resource
import csv
from mtk.geometry.transformation import applyTm, genRotationsByDirections

def read_rotations_from_Sphere():
    Dfd = csv.reader(open(resource("csv", "UnitSphere_500.csv")))
    D = array([ map(float, line) for line in Dfd ])
    R = genRotationsByDirections(D, radians(15.))
    return R

def read_identity():
    return [ identity(4) ]

read_rotations = read_identity

def load_1A2K_b():
    Rl = reader(open(resource("test", "1A2K_l_b.pqr")))
    Rr = reader(open(resource("test", "1A2K_r_b.pqr")))
    Ml = array([ (i['x'], i['y'], i['z'], 1., i['occupancy'])
                for h,i in Rl if h == 'ATOM'])
    Mr = array([ (i['x'], i['y'], i['z'], 1., i['occupancy'])
                for h,i in Rr if h == 'ATOM'])
    resolution = 1.2
    return Mr, Ml, resolution

def load_papain():
    Rl = reader(open(resource("test", "papain.pqr")))
    Rr = reader(open(resource("test", "papain.pqr")))
    Ml = array([ (i['x'], i['y'], i['z'], 1., i['occupancy'])
                for h,i in Rl if h == 'ATOM'])
    Mr = array([ (i['x'], i['y'], i['z'], 1., i['occupancy'])
                for h,i in Rr if h == 'ATOM'])
    resolution = 0.25
    return Mr, Ml, resolution

def load_1A2K_b():
    Rl = reader(open(resource("test", "1A2K_l_b.pqr")))
    Rr = reader(open(resource("test", "1A2K_r_b.pqr")))
    Ml = array([ (i['x'], i['y'], i['z'], 1., i['occupancy'])
                for h,i in Rl if h == 'ATOM'])
    Mr = array([ (i['x'], i['y'], i['z'], 1., i['occupancy'])
                for h,i in Rr if h == 'ATOM'])
    resolution = 0.25
    return Mr, Ml, resolution

def read_piramid():
    Mr = array([ \
                     [ .0, .0, .0, 1., 1. ],\
                    \
                     [ 1., 1., 1., 1., 1. ],\
                     [ 1.,-1., 1., 1., 1. ],\
                     [-1., 1., 1., 1., 1. ],\
                     [-1.,-1., 1., 1., 1. ],\
                    \
                     [ 2., 2., 2., 1., 1. ],\
                     [ 2., 0., 2., 1., 1. ],\
                     [ 2.,-2., 2., 1., 1. ],\
                     [ 0., 2., 2., 1., 1. ],\
                     [ 0.,-2., 2., 1., 1. ],\
                     [-2., 2., 2., 1., 1. ],\
                     [-2., 0., 2., 1., 1. ],\
                     [-2.,-2., 2., 1., 1. ],\
                   ])
    Ml = array([ [ .0, .0, .0, 1., 1. ] ])
    resolution = 1.4
    return Mr, Ml, resolution

def read_half_cube():
    r = 1
    c = [ -r*2.5, 0, r*2.5 ]
    Result = array([ [ c[i], c[j], c[k], 1., r ] for i in range(len(c)) for j in
                    range(len(c)) for k in range(len(c)) ])
    Mr = Result[:len(c)**3/2,:]
    Ml = Result[len(c)**3/2:,:]
    resolution = 0.5
    return Mr, Ml, resolution

def read_dot_cube():
    r = 1
    c = [ -r*2.4, 0, r*2.4 ]
    Result = array([ [ c[i], c[j], c[k], 1., r ]
                    for i in range(len(c))
                    for j in range(len(c))
                    for k in range(len(c))
                    if not (i == 1 and j == 1 and k ==1)
                   ])
    Mr = Result
    Ml = array([ [ 0, 0., 0., 1., r ] ])
    resolution = 0.25
    return Mr, Ml, resolution

read = load_1A2K_b
#read=read_half_cube
