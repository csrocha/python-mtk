# -*- coding: ISO-8859-1 -*-
# $Id: vol.py 88 2009-06-09 15:40:30Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Create and manage volume of proteins

    >>> from numpy import array
    >>> a = array([[ 0, 0, 0 ], [ -1, 2, -3], [5, -6, 7]])
    >>> V = create_from_coords(a, resolution=1.7)
    >>> (V.shape == array([5, 6, 7])).all()
    True
    >>> ("%08.4f "*3) % tuple(V.t(array([-1,-6,-3])))
    '000.0000 000.0000 000.0000 '
    >>> ("%08.4f "*3) % tuple(V.t(array([5,2,7])))
    '003.5294 004.7059 005.8824 '
    >>> V = create_from_coords(a, resolution=1.7, inc=2.0)
    >>> (V.shape == array([7, 9, 10])).all()
    True
    >>> ("%08.4f "*3) % tuple(V.t(array([-1,-6,-3])))
    '001.1765 001.1765 001.1765 '
    >>> ("%08.4f "*3) % tuple(V.t(array([5,2,7])))
    '004.7059 005.8824 007.0588 '
    >>> sum(V.s(array([0,0,0])) == V.min)
    3
    >>> sum(round(V.s(array(V.shape))) >= V.max)
    3
"""

from numpy import array, zeros, round, identity, dot, ceil, floor, ndenumerate
from numpy import trunc, int, allclose
from numpy import any as Any, all as All
from numpy.linalg import norm, inv
from neighbours import d
from mtk.geometry.transformation import applyTv, transform
from sys import stderr
from itertools import starmap

def GridIterator(_min, _max, step=array([1,1,1])):
    """ Iterate in a box of the grid

    >>> list(GridIterator([0,0,0], [1,1,1])) ==\
        [ (0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),\
          (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1) ]
    True
    >>> list(GridIterator(array([1,1,1]), array([2,2,2]))) ==\
        [ (1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2),\
          (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2) ]
    True
    >>> list(GridIterator([1,4,2], [2,5,4])) ==\
        [ (1, 4, 2), (1, 4, 3), (1, 4, 4), (1, 5, 2),\
          (1, 5, 3), (1, 5, 4), (2, 4, 2), (2, 4, 3),\
          (2, 4, 4), (2, 5, 2), (2, 5, 3), (2, 5, 4) ]
    True
    >>> list(GridIterator([0,0,0], [1,1,1], [.5,.5,.5])) ==\
        [(0.0, 0.0, 0.0), (0.0, 0.0, 0.5), (0.0, 0.0, 1.0), (0.0, 0.5, 0.0),\
         (0.0, 0.5, 0.5), (0.0, 0.5, 1.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.5),\
         (0.0, 1.0, 1.0), (0.5, 0.0, 0.0), (0.5, 0.0, 0.5), (0.5, 0.0, 1.0),\
         (0.5, 0.5, 0.0), (0.5, 0.5, 0.5), (0.5, 0.5, 1.0), (0.5, 1.0, 0.0),\
         (0.5, 1.0, 0.5), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0), (1.0, 0.0, 0.5),\
         (1.0, 0.0, 1.0), (1.0, 0.5, 0.0), (1.0, 0.5, 0.5), (1.0, 0.5, 1.0),\
         (1.0, 1.0, 0.0), (1.0, 1.0, 0.5), (1.0, 1.0, 1.0)]
    True
    >>> list(GridIterator([0,0,0], [-1,-1,-1])) == [ ]
    True
    >>> list(GridIterator([0,0,0], [0,0,0])) == [(0, 0, 0)]
    True
    """
    step = array(step)
    _min = array(_min, dtype=step.dtype)
    _max = array(_max, dtype=step.dtype)
    actual = _min.copy()
    step_x = array([step[0], 0, 0])
    step_y = array([0, step[1], 0])
    step_z = array([0, 0, step[2]])

    while all(actual <= _max):
        yield tuple(actual)
        actual += step_z
        if actual[2] > _max[2]:
            actual[2] = _min[2]
            actual += step_y
            if actual[1] > _max[1]:
                actual[1] = _min[1]
                actual += step_x
                if actual[0] > _max[0]:
                    raise StopIteration

def create_from_coords(coords, resolution=None, delta=None, shape=None, inc=0, init_array=zeros, dtype=float):
    """
    Return a volume calculated from the given
    coordinates and resolution.

    To be deprecated.
    """
    try:
        assert(len(inc)==3)
    except:
        inc = (inc,inc,inc)

    from mtk.geometry import origin, extreme

    _min = origin(coords) - inc
    _max = extreme(coords) + inc

    return Volume(_min, _max, resolution=resolution, delta=delta, shape=shape, init_array=init_array, dtype=dtype)

class wrapperVolume:
    def __init__(self, V):
        self._V = V

    def __getitem__(self, key):
        return self._V.__getitem__(key)

    def __setitem__(self, key, v):
        self._V.__setitem__(key, v)

    def get_value(self, key):
        return self._V.get_value(key)

class noLimitsVolume(wrapperVolume):
    """
    Inteface without limits for get values out of the valume.

    >>> A = array([[[ 0.0,0.0 ],[ 0.0,0.0 ]],[[ 0.0,0.0 ],[ 0.0,1.0 ]]], dtype=float)
    >>> V = Volume( (0,0,0), (1,1,1), resolution = 1, def_array = A )
    >>> nLV = noLimitsVolume(V, 256)
    >>> nLV[(0,0,0)]
    0.0
    >>> nLV[(1,1,1)]
    1.0
    >>> nLV[(2,2,2)]
    256
    """
    def __init__(self, V, default=0):
        wrapperVolume.__init__(self, V)
        self.default = default

    def __getitem__(self, key):
        try:
            return self._V[key]
        except IndexError:
            return self.default

class rotationLimitsVolume(wrapperVolume):
    """
    Inteface without limits for get values out of the valume as a thorus.

    >>> A = array([[[ 0.0,0.0 ],[ 0.0,0.0 ]],[[ 0.0,0.0 ],[ 0.0,1.0 ]]], dtype=float)
    >>> V = Volume( (0,0,0), (1,1,1), resolution = 1, def_array = A )
    >>> nLV = rotationLimitsVolume(V)
    >>> nLV[(0,0,0)]
    0.0
    >>> nLV[(1,1,1)]
    1.0
    >>> nLV[(2,2,2)]
    0.0
    >>> nLV[(3,3,3)]
    1.0
    """
    def __init__(self, V):
        wrapperVolume.__init__(self, V)

    def __getitem__(self, key):
        key %= self._V.shape
        return self._V[key]

class BasicVolume:
    def __init__(self, _min, _delta, _data):
        """
        Create a basic volume
        """
        self._min = _min
        self._delta = _delta
        self._data = _data
        self.t = lambda x: (x[:3] - _min) / _delta
        self.s = lambda x: _min + ( array(x) * _delta )

    @property
    def shape(self):
        """Shape of the volume"""
        return array(self._data.shape)

    @property
    def dtype(self):
        """Datatype of the volume"""
        return self._data.dtype

    @property
    def min(self):
        """Minimal coordinate of the grid"""
        return array(self._min)

    @property
    def max(self):
        """Maximal coordinate of the grid"""
        return array(self.min + (self.shape-1)*self.delta)

    @property
    def delta(self):
        """Space beetween points of the grid"""
        return self._delta

    @property
    def center(self):
        return (self.max + self.min)/2.0

    @property
    def size(self):
        return self.shape[0]*self.shape[1]*self.shape[2]

    @property
    def metadata(self):
        MD = {
                'min': self.min,
                'max': self.max,
                'delta': self.delta,
        }
        return MD

    def __getitem__(self, key):
        return self._data[tuple(key)]

    def __setitem__(self, key, value):
        self._data[tuple(key)] = value

    def get_value(self, X, wrapper=noLimitsVolume):
        """
        Trilinear interpolation value calculation
        Ref: http://en.wikipedia.org/wiki/Trilinear_interpolation

        >>> from numpy import allclose
        >>> A = array([[[ 0.0,0.0 ],[ 0.0,0.0 ]],[[ 0.0,0.0 ],[ 0.0,1.0 ]]], dtype=float)
        >>> V = Volume( (0,0,0), (1,1,1), resolution = 1, def_array = A )
        >>> V.get_value((0,0,0))
        0.0
        >>> V.get_value((1,1,1))
        1.0
        >>> allclose(V.get_value((0.5,0.5,0.5)), 0.125)
        True
        >>> allclose(V.get_value((0.999,0.999,0.999)), 0.997)
        True
        >>> V.get_value((1.,1.,0.))
        0.0
        >>> V.get_value((1.,1.,1.))
        1.0
        >>> V.get_value((2.,2.,2.))
        Traceback (most recent call last):
            ...
        IndexError: index (2) out of range (0<=index<2) in dimension 0
        >>> A = array([[[1,0,0],[0,0,0],[0,0,0]],\
                       [[0,0,0],[0,1,0],[0,0,0]],\
                       [[0,0,0],[0,0,0],[0,0,1]]])
        >>> V1 = Volume((-1,-1,-1),(1,1,1), resolution=1, def_array=A)
        >>> all([ V1._data[i] == v for i,v in ndenumerate(V1._data) ])
        True
        >>> all([ V1.get_value(V1.s(i)) == v for i,v in ndenumerate(V1._data) ])
        True
        >>> V1 = Volume((0,0,0),(1,1,1), resolution=0.5, def_array=A)
        >>> all([ V1.get_value(V1.s(i)) == v for i,v in ndenumerate(V1._data) ])
        True

        >>> from mtk.io.vtk_ff import write_vol, writer
        >>> from numpy import allclose, arange
        >>> A = array(range(3*3*3),dtype=float).reshape(3,3,3)
        >>> V = Volume( (0,0,0), (2,2,2), resolution = 1, def_array = A )
        >>> B = [ ((x,y,z),V.get_value((x,y,z))) for z in arange(0, 2.5, 0.5) for y in arange(0, 2.5, 0.5) for x in arange(0, 2.5, 0.5) ]
        >>> A, B
        >>> write_vol(V, 'test.vti')
        >>> Ba, Bb = zip(*B)
        >>> W = writer(points=Ba, vertices=[(1,i) for i in range(len(Ba))], scalars=Bb)
        >>> W.write('test.vtp')
        """
        X = array([X[0],X[1],X[2]], float)
        Xo = (X[:3] - self.min) / self.delta
        Xt = array(trunc(Xo), dtype=int)
        Xd = Xo - trunc(Xo)
        if all(Xd ==(0.0, 0.0, 0.0)):
            return self[tuple(Xo)]
        zd = Xd[2]; ldz = 1-zd
        W = wrapper(self)
        i1 = W[Xt+(0,0,0)]*ldz + W[Xt+(0,0,1)]*zd
        i2 = W[Xt+(0,1,0)]*ldz + W[Xt+(0,1,1)]*zd
        j1 = W[Xt+(1,0,0)]*ldz + W[Xt+(1,0,1)]*zd
        j2 = W[Xt+(1,1,0)]*ldz + W[Xt+(1,1,1)]*zd
        yd = Xd[1]
        w1 = i1*(1 - yd) + i2*yd
        w2 = j1*(1 - yd) + j2*yd
        xd = Xd[0]
        return w1*(1 - xd) + w2*xd

    def subvolume(self, min, max):
        """
        Return a subvolume of the volume. Not check min max values.
        """
        _min = self.t(min)
        _max = self.t(max) + 1
        return self._data[_min[0]:_max[0], _min[1]:_max[1], _min[2]:_max[2]]

    def reminmax(self, n_min, n_max):
        """
        Reset the limits of the Volume.

        >>> A = array([[[ 0.0,0.0,0.0 ],\
                        [ 0.0,0.0,0.0 ],\
                        [ 0.0,0.0,0.0 ]],\
                       [[ 0.0,0.0,0.0 ],\
                        [ 0.0,0.0,0.0 ],\
                        [ 0.0,0.0,0.0 ]],\
                       [[ 0.0,0.0,0.0 ],\
                        [ 0.0,0.0,0.0 ],\
                        [ 0.0,0.0,1.0 ]]], dtype=float)
        >>> V = Volume( (0,0,0), (1,1,1), resolution = 0.5, def_array = A )
        >>> K = [ V.get_value(i) for i in GridIterator(V.min, V.max, V.delta) ]
        >>> V.min, V.max, V.delta, V.shape
        (array([ 0.,  0.,  0.]), array([ 1.,  1.,  1.]), array([ 0.5,  0.5,  0.5]), array([3, 3, 3]))
        >>> V.get_value((1,1,1))
        1.0
        >>> V.shape
        array([3, 3, 3])
        >>> omin, omax = V.min, V.max
        >>> V.reminmax((-1,-1,-1),(2,2,2))
        >>> V.shape
        array([7, 7, 7])
        >>> Q = [ V.get_value(i) for i in GridIterator(omin, omax, V.delta) ]
        >>> V.get_value((2,2,2))
        0.0
        >>> K == Q
        True
        """
        p_min = floor(self.t(n_min))
        p_max = ceil(self.t(n_max))
        n_min = self.s(p_min)
        n_max = self.s(p_max)
        d = floor((n_max - n_min) / self.delta) + 1
        NA = zeros(d, dtype=self.dtype)
        c = (d - self.shape) / 2
        shape = tuple(self.shape)
        NA[c[0]:c[0]+shape[0],
           c[1]:c[1]+shape[1],
           c[2]:c[2]+shape[2]] = self._data
        self._data = NA
        self._min = n_min

    def transfer(self, source, T = identity(4)):
        """
        Copy the content of Volume source with a
        Transformation T.

        >>> A = array([[[0,0,0],[0,0,0],[0,0,0]],\
                       [[0,1,0],[0,2,0],[0,3,0]],\
                       [[0,0,0],[0,0,0],[0,0,0]]])
        >>> V1 = BasicVolume((-1,-1,-1),(1,1,1),A)

        >>> from mtk.geometry.transformation import rotation, pi
        >>> from numpy import radians
        >>> V2 = BasicVolume((-1,-1,-1),(0.5,0.5,0.5),zeros((5,5,5)))
        >>> T = rotation(radians(90), 0, 0)
        >>> V2.transfer(V1, T=T)
        >>> all([ allclose(V2.get_value(applyTv(T.T, V1.s(i).tolist() +\
                [1.])), V1.get_value(V1.s(i))) for i, v in ndenumerate(V1._data) ])
        True
        """
        from mtk.geometry.vol_c import transfer
        transfer(source, self, T)
        #v = array([0,0,0,1],dtype=float)
        #for i,a in ndenumerate(self._data):
        #    v[:3] = self.s(i)
        #    self._data[i] = source.get_value(applyTv(T, v))

    def boxiterator(self, bmin, bmax):
        """
        Iterate over a box in nodes between bmin to bmax.
        For iteration return the coordinate in space (x), the
        index in the array (i) and the value in the array (i)

        >>> from numpy import all
        >>> V = BasicVolume((-5,-5,-5),(0.5,0.5,0.5),zeros((20,20,20)))
        >>> S = [ x for x, i, v in V.boxiterator((-3,-3,-3),(3,3,3)) ]
        >>> len(S) == 13*13*13
        True
        """
        #import pdb; pdb.set_trace()
        imin = self.t(array(bmin))
        imax = self.t(array(bmax)) + 1
        slicing = map(lambda (a,b):
                      slice(a,b),zip(imin,imax))
        A = self._data[slicing]
        s = self.s
        for i,v in ndenumerate(A):
            yield (bmin+array(i)*self._delta, imin + i, v)

    def fornodes(self, F, _min, _max):
        slicing = map(lambda (a,b): slice(a,b),zip(_min,_max))
        A = self._data[slicing]
        _min = array(_min)
        s = self.s
        for i,v in ndenumerate(A):
            self._data[slicing][i] = F(s(_min+i), v)

    def _fornodes(self, F, _min, _max):
        """
        Repeat a function on a list of values in a cubic range
        """
        for i in GridIterator(array(_min), array(_max)):
            self._data[i] = F(self.s(i), self._data[i])

    def forallnodes(self, F):
        """
        Repeat a function on a list of values in all the volume

        >>> from numpy import all
        >>> V = Volume((0,0,0), (10,10,10), resolution=1)
        >>> V.forallnodes(lambda x, y: 1)
        >>> all(sum(V._data) == 11.)
        True
        """
        self.fornodes(F, (0,0,0), array(self.shape))

    def dump(self, filename):
        """
        Almacena el volumen en un archivo vti
        """
        from mtk.io.vtk_ff import writer

        if filename[-4:] != '.vti':
            return False

        W = writer(origin = self.min, spacing = self.delta, points = self._data)
        W.write(filename)

        return True

    def old_dump(self, VOLfilename, x_dict={}):
        """
        Almacena el volumen en un archivo.
        """
        import tarfile
        import pickle
        import numpy
        from StringIO import StringIO
        from tempfile import TemporaryFile

        volfile = tarfile.open(VOLfilename, 'w:gz')

        # Store metadata
        mdIO = StringIO()
        MD = self.metadata
        MD.update(x_dict)
        pickle.dump(MD, mdIO)
        mdIO.seek(0)
        info = tarfile.TarInfo('md')
        info.size = len(mdIO.buf)
        volfile.addfile(info, mdIO)

        # Store Array
        volIO = TemporaryFile()
        numpy.save(volIO, self._data)
        info = tarfile.TarInfo('grid')
        info.size = volIO.tell()
        volIO.seek(0)
        volfile.addfile(info, volIO)

        mdIO.close()
        volIO.close()
        volfile.close()

class Volume(BasicVolume):
    def __init__(self, _min, _max, resolution=None, delta=None, shape=None, init_array=zeros, def_array=None, dtype=float):
        """ Return a volume calculated from the given
            coordinates and resolution.

        >>> V = Volume((-1,-1,-1),(1,1,1),resolution=1.)
        >>> V.t((-1,-1,-1))
        array([ 0.,  0.,  0.])
        >>> V.t((1,1,1))
        array([ 2.,  2.,  2.])
        >>> V.shape
        array([3, 3, 3])
        >>> V._data.shape
        (3, 3, 3)
        >>> V._data[2,2,2]
        0.0
        >>> V = Volume( (0,0,0), (1,1,1), resolution = 0.5,\
                       def_array = zeros((2,2,2)) )
        Traceback (most recent call last):
            ...
        RuntimeError: min, max and resolution not related to array shape
        """
        size = array(_max, dtype=float) - array(_min, dtype=float)
        if delta != None and len(delta) != 3:
            raise RuntimeEror("Delta must be an array of 3 floats")
        if delta != None:
            shape = array(size / delta, dtype=int) + 1
        elif resolution != None:
            delta = array([resolution]*3)
            shape = array(ceil(size / delta), dtype=int) + 1
        elif shape != None:
            delta = array(size / shape, dtype=float) + 1
        else:
            raise RuntimeError("Not define one of delta, resolution or shape values. Can't create the grid")

        if def_array != None and all(shape != def_array.shape):
            raise RuntimeError("min, max and resolution not related to array shape")

        _min = array(_min[:3], dtype=float)
        _max = array(_max[:3], dtype=float)
        _delta = array(delta, dtype=float)
        if def_array != None:
            def_array = array(def_array)
            _data = def_array
        else:
            _data = init_array(shape, dtype)
        BasicVolume.__init__(self, _min, _delta, _data)

#    def get_nodes(self):
#        """
#        return an iterator over the volume
#        """
#        return ndindex(self.shape)
#
#    def get_coords(self):
#        """
#        return an iterator over the coordinates of the nodes in the volume
#        """
#        return starmap(self.s, self.get_nodes())
#
    def put_ball(self, x, r, f):
        """ Put a ball with center x and r radius on the Volume using the f writter function.
            f ( r, v ) : float
            r: distance to the center of the atom
            v: value of the Volume in this position before writting.

        >>> from numpy import array
        >>> a = array([[ 0, 0, 0 ], [ -1, 2, -3], [5, -6, 7]])
        >>> V = create_from_coords(a, resolution=1.7, inc=2.0)
        >>> V.put_ball(a[0], 2.0, lambda r,v: r)
        >>> map(lambda x: "%5.4f" % x, V._data[0:5,3:8,1:6].flatten())
        ['0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '2.3854', '1.7720', '2.5239', '0.0000', '0.0000', '2.1213', '1.3964', '2.2760', '0.0000', '0.0000', '3.0150', '2.5573', '3.1257', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '2.0396', '1.2689', '2.2000', '0.0000', '0.0000', '1.7234', '0.6481', '1.9105', '0.0000', '0.0000', '2.7495', '2.2383', '2.8705', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '2.9000', '2.4207', '3.0150', '0.0000', '0.0000', '2.6870', '2.1610', '2.8107', '0.0000', '0.0000', '3.4366', '3.0430', '3.5341', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000']
        """
        r = array([r]*3)
        _min = round(self.t(x[:3] - r))
        _max = round(self.t(x[:3] + r))

        F = lambda o, v: f(norm(o - x[:3]), v)

        self._fornodes(F, _min, _max)

def load(filename):
    from mtk.io.vtk_ff import reader

    if filename[-4:] != '.vti':
        return False

    R = reader()
    D, I = R.read(filename)
    return BasicVolume(I.origin, I.spacing, D['scalars'])

def old_load(VOLfilename):
    """
    Lee el volumen de un archivo.
    """
    import tarfile
    import pickle
    import numpy
    volfile = tarfile.open(VOLfilename, 'r')
    MD = pickle.load(volfile.extractfile('md'))
    A = numpy.load(volfile.extractfile('grid'))
    return BasicVolume(MD['min'], MD['delta'], A), MD

def cloneVolume(V, init_array=None):
    """
    Copy a Volume
    """
    if init_array != None:
        return Volume(V.min, V.res, init_array = init_array)
    else:
        return Volume(V.min, V.delta, V._data.copy())

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4:

