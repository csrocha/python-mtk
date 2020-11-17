# -*- coding: ISO-8859-1 -*-
# $Id: __init__.py 64 2009-05-19 17:29:15Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
""" Plugins for Storage queries
"""

def find_functions(path=None):
    """
    Find modules of plugins
    """
    import sys, os, fnmatch, string

    #LazyModule = sys.modules['.']

    # Walk the current directory to identify all management modules,
    # and store a reference to them within the returned dictionary.

    if path == None:
        path = __path__[0]

    modfiles = filter(lambda x,m=fnmatch.fnmatch:m(x,'*.py'),
                      os.listdir(path))

    functions = []
    for curmodule in modfiles:
        modname = os.path.splitext(curmodule)[0]
        if modname[0] == '_':
            continue
        else:
            M = __import__(modname, globals(), locals())
            functions +=[ (f[5:], vars(M)[f]) for f in dir(M) if f[:5] == 'sqlf_' ]
    return dict(functions)

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

