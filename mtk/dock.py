# -*- coding: ISO-8859-1 -*-
#
# $Id: vol.py 85 2009-06-07 23:46:26Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""

Docking testing.

Este test trata de resolver el problema de docking, pero la verdad es que no
es la mejor forma de testearlo. Hay que revisar los números finales.

>>> from numpy import array, argmax, argmin, all, identity, allclose
>>> from mtk.geometry.vol import BasicVolume
>>> c, s = 9j, 1

Testing del docking en el eje Z

>>> Ad = array([ \
               [ \
                [ c, s ],\
                [ 0, 0 ],\
               ],[ \
                [ 0, 0 ],\
                [ 0, 0 ],\
               ] \
              ], dtype=complex)
>>> Bd = array([ \
               [ \
                [ s, c ],\
                [ 0, 0 ],\
               ],[ \
                [ 0, 0 ],\
                [ 0, 0 ],\
               ] \
              ], dtype=complex)
>>> A = zeros((4,4,4), dtype=complex)
>>> B = zeros((4,4,4), dtype=complex)
>>> A[1:3,1:3,1:3] = Ad
>>> B[1:3,1:3,1:3] = Bd
>>> Va = BasicVolume((0,0,0), (1,1,1), A)
>>> Vb = BasicVolume((0,0,0), (1,1,1), B)
>>> fVb = lambda t: Vb
>>> C = Collector()
>>> dock(Va, fVb, [ identity(4) ], progress_function=C)
True
>>> Dr = C.score.real
>>> allclose((min(Dr), max(Dr)), (c*c*1, s*s*1))
True
>>> C.translation[argmax(Dr)]
array([ 0.,  0.,  1.])
>>> C.translation[argmin(Dr)]
array([ 0.,  0., -1.])

Testing del docking en el eje Y

>>> Ad = array([ \
               [ \
                [ c, 0 ],\
                [ s, 0 ],\
               ],[ \
                [ 0, 0 ],\
                [ 0, 0 ],\
               ] \
              ], dtype=complex)
>>> Bd = array([ \
               [ \
                [ s, 0 ],\
                [ c, 0 ],\
               ],[ \
                [ 0, 0 ],\
                [ 0, 0 ],\
               ] \
              ], dtype=complex)
>>> A = zeros((4,4,4), dtype=complex)
>>> B = zeros((4,4,4), dtype=complex)
>>> A[1:3,1:3,1:3] = Ad
>>> B[1:3,1:3,1:3] = Bd
>>> Va = BasicVolume((0,0,0), (1,1,1), A)
>>> Vb = BasicVolume((0,0,0), (1,1,1), B)
>>> fVb = lambda t: Vb
>>> C = Collector()
>>> dock(Va, fVb, [ identity(4) ], progress_function=C)
True
>>> Dr = C.score.real
>>> allclose((min(Dr), max(Dr)), (c*c*1, s*s*1))
True
>>> C.translation[argmax(Dr)]
array([ 0.,  1.,  0.])
>>> C.translation[argmin(Dr)]
array([ 0., -1.,  0.])

Testing del docking en el eje X

>>> Ad = array([ \
               [ \
                [ c, 0 ],\
                [ 0, 0 ],\
               ],[ \
                [ s, 0 ],\
                [ 0, 0 ],\
               ] \
              ], dtype=complex)
>>> Bd = array([ \
               [ \
                [ s, 0 ],\
                [ 0, 0 ],\
               ],[ \
                [ c, 0 ],\
                [ 0, 0 ],\
               ] \
              ], dtype=complex)
>>> A = zeros((4,4,4), dtype=complex)
>>> B = zeros((4,4,4), dtype=complex)
>>> A[1:3,1:3,1:3] = Ad
>>> B[1:3,1:3,1:3] = Bd
>>> Va = BasicVolume((0,0,0), (1,1,1), A)
>>> Vb = BasicVolume((0,0,0), (1,1,1), B)
>>> fVb = lambda t: Vb
>>> C = Collector()
>>> dock(Va, fVb, [ identity(4) ], progress_function=C)
True
>>> Dr = C.score.real
>>> allclose((min(Dr), max(Dr)), (c*c*1, s*s*1))
True
>>> C.translation[argmax(Dr)]
array([ 1.,  0.,  0.])
>>> C.translation[argmin(Dr)]
array([-1.,  0.,  0.])


"""

from numpy import array, ones, zeros, ndenumerate, nonzero, transpose, append, vstack
from numpy.fft import fftn, ifftn, rfftn, irfftn, fftshift
from time import time

def dock(Va, fVb, Ts, filter_function=lambda v: ones(v.shape, dtype=bool), progress_function=None):
    """
    Realiza el docking.

    Devuelve:

        dR: id de la matriz de rotacion.
        dT: traslacion con respecto al centro del receptor.
        dS: score resultante.
    """

    if progress_function: progress_function("start", time)

    if progress_function: progress_function("receptor copy", Va)

    # DISCUSION:
        #  Que significa que la asignacion de conjugados tal como lo plantie
        #  aqui?
    A = Va._data #.conjugate()
    N = Va.shape[0] * Va.shape[1] * Va.shape[2]

    if progress_function: progress_function("receptor fft", A)

    Fa = fftn(A) #.conjugate()

    if progress_function: progress_function("rotation loop", Fa)

    # Inicia la busqueda en la rotación
    for ir in xrange(len(Ts)):
        if progress_function: progress_function("rotation", ir)

        T = Ts[ir]
        Vb = fVb(T)
        B = Vb._data.conjugate()

        if progress_function: progress_function("dock", B)

        R = ifftn(Fa * fftn(B).conjugate())

        if progress_function: progress_function("results selection", (ir, R))

        base = Va.center - Va.min

        # Identificamos los elementos que queremos
        condition = filter_function(R)
        # Tomamos los scores
        dS = R[condition]
        if dS.shape[0] > 0: # chequeamos que hay resultados.
            # Conseguimos sus coordenadas
            idx = transpose(nonzero(condition))
            recoord = (idx > (Va.shape-1)/2)
            dT = (recoord * (idx - Va.shape) + (1-recoord) * idx)*Va.delta
            # Indicamos la rotacion
            dR = ones(len(dS)) * ir
        else:
            dT = zeros((0,3))
            dR = zeros((0))

        if progress_function: progress_function("docked", (time(), ir, (dR, dT, dS)))
        del dR
        del dT
        del dS

    if progress_function: progress_function("finish", time)

    return True

class Collector:
    def __init__(self):
        self._R = zeros(0)
        self._T = zeros((0,3))
        self._S = zeros(0)

    def __call__(self, state, data):
        if state == "docked":
            self._R = append(self._R, data[2][0])
            self._T = vstack((self._T, data[2][1]))
            self._S = append(self._S, data[2][2])

    @property
    def score(self):
        return self._S

    @property
    def rotation(self):
        return self._R

    @property
    def translation(self):
        return self._T

def save((R, T, S), filename):
    import csv

    f = open(filename, 'w')

    c = csv.writer(f)

    c.writerow(("Rotation id", "X", "Y", "Z", "Score"))
    for i in xrange(len(R)):
        c.writerow((R[i],)+ tuple(T[i]) + (S[i],))
    f.close()

def load(filename):
    import csv

    f = open(filename)

    c = csv.reader(f)

    R, T, S = [], [], []

    c.next() # Skip header
    for r,x,y,z,s in c:
        #print r, x, y, z, s
        R.append(int(float(r)))
        T.append((float(x),float(y),float(z)))
        S.append(complex(s))
    f.close()

    return R,T,S

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:CompilerSet makeprg=python\ %<.py

