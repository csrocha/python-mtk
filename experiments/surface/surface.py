# -*- coding: ISO-8859-1 -*-
# $Id: vol.py 89 2009-06-09 16:35:56Z cristian_docking $

from structures import *
from numpy import array, ones, zeros, save, sum
from numpy.fft import fftn, ifftn, rfftn, irfftn, fftshift
from mtk.geometry.vol import BasicVolume
from mtk.view.vol import Viewer
from mtk.geometry.vol_c import put_ball
from mtk.geometry import origin, centre
from mtk.signal import xcorr

shape = (250,250,250)
mincoord = array(shape)/-2
delta = (0.5,0.5,0.5)

class Objeto:
    def filled(self):
        pass

    def surface(self):
        pass

class Ball(Objeto):
    def filled(self, A):
        put_ball(A, array((0.,0.,0.)), {10:1}, [0,1]);

    def surface(self, A):
        put_ball(A, array((0.,0.,0.)), {10:-1, 11:1}, [0,1,-1]);


class TwoBall(Objeto):
    def filled(self, A):
        put_ball(A, array((-6.,0.,0.)), {10:1}, [0,1]);
        put_ball(A, array((6.,0.,0.)), {10:1}, [0,1]);

    def surface(self, A):
        put_ball(A, array((-6.,0.,0.)), {10:-1, 11:1}, [0,1,-1]);
        put_ball(A, array((6.,0.,0.)), {10:-1, 11:1}, [0,1,-1]);

class Cube(Objeto):
    def filled(self, A):
        A._data[40:60,40:60,40:60] = ones((20,20,20))

    def surface(self, A):
        A._data[40:60,40:60,40:60] = ones((20,20,20))
        A._data[41:59,41:59,41:59] = -ones((18,18,18))

class Molecule(Objeto):
    def __init__(self):
        print "Iniciando molecula"
        atoms, L, r = load_1A2K_b()
        n = 1000
        zoom = 2.
        self.coords = atoms[:n,:4]*zoom
        self.radios = atoms[:n,4]*zoom
        self.centre = centre(self.coords[:,:3])

    def filled(self, A):
        print "Generando volumen (%i)" % self.coords.shape[0]
        coords = self.coords
        radios = self.radios
        order_list = [0., 1.]
        for i in xrange(coords.shape[0]):
            x = coords[i][:3] - self.centre
            r = radios[i]
            assign_map = { r: 1 }
            if (i % 100 == 0): print i, x
            put_ball(A, x, assign_map, order_list)

    def surface(self, A):
        print "Generando superficie"
        coords = self.coords
        radios = self.radios
        order_list = [0,1,-1]
        for i in xrange(coords.shape[0]):
            x = coords[i][:3] - self.centre
            r = radios[i]
            assign_map = { r-1: -1, r: 1 }
            if (i % 100 == 0): print i, x
            put_ball(A, x, assign_map, order_list)

class SoftMolecule(Objeto):
    def __init__(self):
        print "Iniciando molecula"
        atoms, L, r = load_1A2K_b()
        n = 1000
        zoom = 2.
        self.coords = atoms[:n,:4]*zoom
        self.radios = atoms[:n,4]*zoom
        self.centre = centre(self.coords[:,:3])

    def filled(self, A):
        print "Generando volumen (%i)" % self.coords.shape[0]
        coords = self.coords
        radios = self.radios
        order_list = [0., 1.]
        for i in xrange(coords.shape[0]):
            x = coords[i][:3] - self.centre
            r = radios[i]
            assign_map = { r: 1 }
            if (i % 100 == 0): print i, x
            put_ball(A, x, assign_map, order_list)

        B = BasicVolume(mincoord, delta, zeros(A.shape))
        put_ball(B, (0,0,0), { 2: 1 }, order_list)

        C = xcorr(A._data, B._data.conjugate())
        View.append(C)


    def surface(self, A):
        print "Generando superficie"
        coords = self.coords
        radios = self.radios
        order_list = [0,1,-1]
        for i in xrange(coords.shape[0]):
            x = coords[i][:3] - self.centre
            r = radios[i]
            assign_map = { r-1: -1, r: 1 }
            if (i % 100 == 0): print i, x
            put_ball(A, x, assign_map, order_list)


Objs = {
#    'ball.npy': Ball(),
#    'twoball.npy': TwoBall(),
#    'cube.npy': Cube(),
    'molecule.npy': Molecule(),
#    'softmolecule.npy': SoftMolecule(),
}


View = Viewer()
for filename, Obj in Objs.items():

    # El objeto
    A = BasicVolume(mincoord, delta, zeros(shape))
    Obj.filled(A)

    # El calculador del área
    B = BasicVolume(mincoord, delta, zeros(shape))
    put_ball(B, array((0.,0.,0.)), {0.70:1}, [0,1]);
    nnodes = sum(B._data.flatten)

    # El filtro
    C = BasicVolume(mincoord, delta, zeros(shape))
    Obj.surface(C)
    C._data = array(C._data > 0, dtype=complex)

    Ad = A._data
    Bd = B._data.conjugate()
    Rd = fftshift(ifftn(fftn(Ad) * fftn(Bd).conjugate())) #/nnodes

    R = BasicVolume(mincoord, delta, Rd)

    RC = BasicVolume(mincoord, delta, R._data * C._data)

#    save(filename, RC._data)

    View.append(A)
    View.append(B)
    View.append(C)
    View.append(R)
    View.append(RC)

View.run()

