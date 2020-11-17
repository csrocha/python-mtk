# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
# Unit Tests for intersections
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

import sys, os
import unittest

sys.path.insert(0, os.path.dirname(__file__))

import mtk.geometry.intersection as intersection

from mtk.geometry.sphere import sphere
from mtk.geometry.triangle import triangle
from numpy import allclose, pi

class TestGeometryIntersection_TriangleSphere(unittest.TestCase):
        def setUp(self):
                pass

        def test_simplesurface(self):
                s = sphere([0,0,0], 5.0)
                T = [ triangle([0, 0, 0],[10,-10,0],[0,-10, 0]),
                      triangle([0, 0, 0],[0,-10,0],[-10,-10,0]),
                      triangle([0, 0, 0],[-10,-10,0],[-10,10,0]),
                      triangle([0, 0, 0],[-10,10,0],[0,10,0]),
                      triangle([0, 0, 0],[0,10,0],[10,10,0]),
                      triangle([0, 0, 0],[10,10,0],[10,-10,0]),

                      triangle([0, 1, 2],[10,-11,1],[0,-10, 1]),
                      triangle([0, 1, 2],[0,-10,1],[-11,-10,0]),
                      triangle([0, 1, 2],[-11,-10,0],[-10,11,0]),
                      triangle([0, 1, 2],[-10,11,0],[0,8,0]),
                      triangle([0, 1, 2],[0,8,0],[11,10,-1]),
                      triangle([0, 1, 2],[11,10,-1],[10,-11,0]),
 
                      triangle([0, 0, 0],[10,-10,1],[0,-10, 0]),
                      triangle([0, 0, 0],[0,-10,0],[-10,-10,0]),
                      triangle([0, 0, 0],[-10,-10,0],[-10,10,0]),
                      triangle([0, 0, 0],[-10,10,0],[0,10,0]),
                      triangle([0, 0, 0],[0,10,0],[10,10,1]),
                      triangle([0, 0, 0],[10,10,1],[10,-10,1]),
                    ]
                for t in T:
                        for a,o,r,pa,pb in intersection.triangle_sphere(t, s):
                                # Distance of intersection points is cero to the sphere surface
                                self.assertTrue(allclose( s.dist(pa), 0.),
                                                '%s.dist(%s)=%f and must be zero' %
                                                (t, pa, s.dist(pa)))
                                self.assertTrue(allclose( s.dist(pb), 0.),
                                                '%s.dist(%s)=%f and must be zero' %
                                                (t, pb, s.dist(pb)))
                                # Distance of intersection points is cero to the plane associated to the triangle
                                self.assertTrue(allclose(t.plane().dist(pa), 0.),
                                                '%s.dist(%s)=%f and must be zero' %
                                                (t, pa, t.plane().dist(pa)))
                                self.assertTrue(allclose(t.plane().dist(pb), 0.),
                                                '%s.dist(%s)=%f and must be zero' %
                                                (t, pb, t.plane().dist(pb)))
                                # Origin of small circle is in 0,0,0
                                #self.assertTrue( allclose( o, [0., 0., 0.] ),
                                #                'Incorrect center (%s)' % o)
                                # Value of a is a logic range
                                self.assertTrue(0 <= a and a <= 2*pi,
                                                'Angle out of range (%f)' % a)
 
        def test_complexsurface(self):
                s = sphere([-0.923142, -0.045219,  0.997029], 0.02)
                T = [ triangle([-0.908749, -0.045065,  1.024429], [-0.90862  ,-0.055086 , 1.024449], [-0.908307 ,-0.045029  ,1.038119]),
                    ]
                for t in T:
                        for a,o,r,pa,pb in intersection.triangle_sphere(t, s):
                                # Value of a is a logic range
                                self.assertTrue(0 <= a and a < 2*pi,
                                                'Angle out of range (%f)' % a)
                                # Distance of intersection points is cero to the sphere surface
                                self.assertTrue(allclose( s.dist(pa), 0.),
                                                '%s.dist(%s)=%f and must be zero' %
                                                (t, pa, s.dist(pa)))
                                self.assertTrue(allclose( s.dist(pb), 0.),
                                                '%s.dist(%s)=%f and must be zero' %
                                                (t, pb, s.dist(pb)))
                                # Distance of intersection points is cero to the plane associated to the triangle
                                self.assertTrue(allclose(t.plane().dist(pa), 0.),
                                                '%s.dist(%s)=%f and must be zero' %
                                                (t, pa, t.plane().dist(pa)))
                                self.assertTrue(allclose(t.plane().dist(pb), 0.),
                                                '%s.dist(%s)=%f and must be zero' %
                                                (t, pb, t.plane().dist(pb)))
 
if __name__ == '__main__':
        unittest.main()

# vim:expandtab:smartindent:tabstop=4:softtabstop=4:shiftwidth=4
