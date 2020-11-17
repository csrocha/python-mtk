# -*- coding: ISO-8859-1 -*-
# $Id: storage.py 109 2009-07-28 23:18:52Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
""" Modulo de almacenamiento de mtk

>>> S = Storage()
"""

import logging

logging.basicConfig(format='%(levelname)s:%(asctime)-15s:%(process)s:%(module)s:%(message)s', filename='/tmp/python-mtk.log')
log = logging.getLogger()
log.setLevel(logging.DEBUG)

