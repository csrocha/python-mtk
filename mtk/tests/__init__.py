# -*- coding: ISO-8859-1 -*-
# $Id: __version__.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.

import mtk
from os.path import join
testpath = mtk.__path__ + [ "data", "test" ]

# Check if its the best way
# path = resource_filename('mtk', '')
#
#

def testfile(filename):
    return join(*(testpath + [filename]))

def opentestfile(filename):
    return open(testfile(filename))

