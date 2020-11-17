# -*- coding: ISO-8859-1 -*-
# $Id$
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Create and manage volume of protein using cuda

>>> 1
1
"""

from vol import Volume

def have_gpu():
    try:
        import pycuda.autoinit
        return True
    except:
        return False

if have_gpu():
    import pycuda.gpuarray as gpuarray

    import pycuda.driver as drv
    from pycuda.compiler import SourceModule

    import pycuda.autoinit
    assert isinstance(pycuda.autoinit.device.name(), str)
    assert isinstance(pycuda.autoinit.device.compute_capability(), tuple)
    assert isinstance(pycuda.autoinit.device.get_attributes(), dict)
else:
    raise RuntimeError("No cuda support!")

from jinja2 import Template
from numpy import array, zeros, round, identity, dot, ceil
from numpy.linalg import norm, inv

class VolumeCUDA(Volume):
    def __init__(self, min, max, resolution=None, delta=None, shape=None, init_array=zeros, def_array=None, dtype=float):
        """
        """
        Volume.__init__(self, min, max, resolution=resolution, delta=delta, shape=shape, init_array=init_array, def_array=def_array, dtype=dtype)
        V.A = to_gpu(V.A)
        V.t = Template("""
            __global__ void apply_ball( {{ type_name }} *dest, {{ type_name }} *A)
            {
                const int i = threadIdx.x;
                const int j = threadIdx.y;
                const int k = threadIdx.z;
                const {{type_name}} r = sqrt(x*x + y*y + z*z);

                dest = {{ function }}(r, i, j, k) 
            }
            """)

    def fornodes(self, F, min, max):
        """
        """
        assert isinstance(F, SourceModule)
        dest = numpy.zeros_like(V.A)

        F(drv.Out(drv.Out(dest), drv.In(V.A)))

        V.A = dest

    def put_ball(self, x, r, f):
        """
        """

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

