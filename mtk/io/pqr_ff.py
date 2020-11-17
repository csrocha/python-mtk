# -*- coding: ISO-8859-1 -*-
# $Id: pdb_ff.py 81 2009-06-02 18:55:20Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""PQR file format processor module.

Reading a pqr fileformat from a string as a stream. The iterator of
PQR entries are stored in p as a dict.

Read standard input format

>>> import StringIO
>>> atomline  = "REMARK FILENAME=\\\"somewhere\\\"\\n"
>>> atomline += "ATOM   2989  CB  THR   197      -8.881 -18.844  28.383  0.3654 1.9080\\n\\n"
>>> atomfile = StringIO.StringIO(atomline)
>>> p = reader(atomfile)
>>> p.next() == ('REMARK', {'comment': 'FILENAME="somewhere"',\
                            'record': 'REMARK', 'pdbid': ''})
True
>>> p.next() == ('ATOM', {'resName': 'THR', 'chainID': '', \
                          'name': 'CB', 'record': 'ATOM', 'altLoc': '', \
                          'occupancy': 1.9079999999999999, 'resSeq': 197, \
                          'charge': 0.3654, 'iCode': '', 'pdbid': '', \
                          'y': -18.844000000000001, 'x': -8.8810000000000002, \
                          'serial': 2989, 'z': 28.382999999999999})
True

Read standard input alternative

>>> atomline  = "ATOM   2989  CB  THR   197      -8.881 -18.844  28.383  0.3654 1.9080\\n"
>>> atomline += "ATOM   2990  OG1 THR   197      -9.195 -20.211  28.098 -0.6761 1.7210\\n"
>>> atomline += "ATOM   2991  CG2 THR   197      -7.462 -18.833  28.918 -0.2438 1.9080\\n"
>>> atomline += "ATOM   2992  HA  THR   197      -9.812 -17.243  29.353  0.1007 1.3870\\n"
>>> atomline += "ATOM   2993 HG22 THR   197      -6.949 -19.611  28.545  0.0642 1.4870\\n"
>>> atomline += "ATOM   2994 HG21 THR   197      -7.002 -17.981  28.652  0.0642 1.4870\\n"
>>> atomfile = StringIO.StringIO(atomline)
>>> p = reader(atomfile)
>>> p.next() == ('ATOM', {'resName': 'THR', 'chainID': '',\
                          'name': 'CB', 'record': 'ATOM', \
                          'altLoc': '', 'occupancy': 1.9079999999999999, \
                          'resSeq': 197, 'charge': 0.3654, 'iCode': '', \
                          'pdbid': '', 'y': -18.844000000000001, \
                          'x': -8.8810000000000002, 'serial': 2989, \
                          'z': 28.382999999999999})
True

Write the input to the standard output without any modification

>>> import sys
>>> atomfile.seek(0)                   # Whe need restart the file
>>> sys.stdout.writelines(writter(p))
ATOM   2989  CB  THR   197      -8.881 -18.844  28.383  0.3654 1.9080
ATOM   2990  OG1 THR   197      -9.195 -20.211  28.098 -0.6761 1.7210
ATOM   2991  CG2 THR   197      -7.462 -18.833  28.918 -0.2438 1.9080
ATOM   2992  HA  THR   197      -9.812 -17.243  29.353  0.1007 1.3870
ATOM   2993 HG22 THR   197      -6.949 -19.611  28.545  0.0642 1.4870
ATOM   2994 HG21 THR   197      -7.002 -17.981  28.652  0.0642 1.4870

Get the coords as an array of affine vectors

>>> atomfile.seek(0)                   # Whe need restart the file
>>> c = p.get_coords()
>>> c
array([[ -8.881, -18.844,  28.383,   1.   ],
       [ -9.195, -20.211,  28.098,   1.   ],
       [ -7.462, -18.833,  28.918,   1.   ],
       [ -9.812, -17.243,  29.353,   1.   ],
       [ -6.949, -19.611,  28.545,   1.   ],
       [ -7.002, -17.981,  28.652,   1.   ]])

Write modified protein coordinates to the standard output.

>>> from numpy import ones
>>> i = ones(c.shape)
>>> atomfile.seek(0)                   # Whe need restart the file
>>> q = update_coords(p, i)
>>> sys.stdout.writelines(writter(q))
ATOM   2989  CB  THR   197       1.000   1.000   1.000  0.3654 1.9080
ATOM   2990  OG1 THR   197       1.000   1.000   1.000 -0.6761 1.7210
ATOM   2991  CG2 THR   197       1.000   1.000   1.000 -0.2438 1.9080
ATOM   2992  HA  THR   197       1.000   1.000   1.000  0.1007 1.3870
ATOM   2993 HG22 THR   197       1.000   1.000   1.000  0.0642 1.4870
ATOM   2994 HG21 THR   197       1.000   1.000   1.000  0.0642 1.4870

>>> atomline =  "ATOM   4052  HE3 LYS   307      -7.982  85.475 -165.252  0.1135 1.1000"
>>> atomline += "ATOM   4053  HG2 LYS   307     -10.970  85.337 -165.256  0.0103 1.4870"
>>> atomline += "ATOM   4054  HG3 LYS   307      -9.736  84.396 -164.747  0.0103 1.4870"
>>> atomline += "ATOM   4055  HZ1 LYS   307      -6.613  85.109 -166.951  0.3400 0.6000"
>>> atomline += "ATOM   4056  HZ3 LYS   307      -7.559  85.924 -168.005  0.3400 0.6000"
>>> atomline += "ATOM   4057  HZ2 LYS   307      -6.624  86.743 -166.945  0.3400 0.6000"
>>> atomfile = StringIO.StringIO(atomline)
>>> p = reader(atomfile)
>>> p.next() == ('ATOM', {'resName': 'LYS', 'chainID': '', 'name': 'HE3',\
                          'record': 'ATOM', 'altLoc': '',\
                          'occupancy': 1.1000000000000001, 'resSeq': 307,\
                          'charge': 0.1135, 'iCode': '', 'pdbid': '', \
                          'y': 85.474999999999994, 'x': -7.9820000000000002,\
                          'serial': 4052, 'z': -165.25200000000001})
True
"""

import re
from numpy import array

records = {
   "REMARK": {
        "re": re.compile(
            r"^"
            r"(?P<record>REMARK) "
            r"(?P<comment>.*)"
            r"\s*"
            r"$"
            ),
        "cn": { "record": lambda x: str(x).strip(),
            "comment": lambda x: str(x).strip(),
            },
        "fm": "REMARK "
            "%(comment)s"
            "\n"
        },
   "ATOM": {
        "re": re.compile(
            r"^"
            r"(?P<record>ATOM)  "
            r"(?P<serial>[\d ]{5})"
            r".{1}"
            r"(?P<name>[\w ]{4})"
            r"(?P<altLoc>[\w ]{1})"
            r"(?P<resName>[\w ]{3})"
            r".{1}"
            r"(?P<chainID>[\w ]{1})"
            r"(?P<resSeq>.{4})"
            r"(?P<iCode>.{1})"
            r"\s+"
            r"(?P<x>[\-\+\d\.]+)"
            r"\s+"
            r"(?P<y>[\-\+\d\.]+)"
            r"\s+"
            r"(?P<z>[\-\+\d\.]+)"
            r"\s+"
            r"(?P<charge>[\-\+\d\.]+)"
            r"\s+"
            r"(?P<occupancy>[\-\+\d\.]+)"
            r".*$"
            ),
        "cn": { "record": lambda x: str(x).strip(),
            "serial": lambda x: int(str(x).strip()),
            "name": lambda x: str(x).strip(),
            "altLoc": lambda x: str(x).strip(),
            "resName": lambda x: str(x).strip(),
            "chainID": lambda x: str(x).strip(),
            "resSeq": lambda x: int(str(x).strip()),
            "iCode": lambda x: str(x).strip(),
            "x": lambda x: float(str(x).strip()),
            "y": lambda x: float(str(x).strip()),
            "z": lambda x: float(str(x).strip()),
            "charge": lambda x: float(str(x).strip()),
            "occupancy": lambda x: float(str(x).strip()),
            },
        "fm": "ATOM  "
            "%(serial)5i"
            " "
            "%(name)-4s"
            "%(altLoc)1s"
            "%(resName)-3s"
            " "
            "%(chainID)1s"
            "%(resSeq)4i"
            "%(iCode)1s"
            "   "
            "%(x)8.3f"
            "%(y)8.3f"
            "%(z)8.3f"
            "%(charge)8.4f"
            "%(occupancy)7.4f"
            "\n"
        },
    }

class reader:
    """Class for iterate over a PQR input

    Note: Dont recognize records out of the scope.
    """
    def __init__(self, input, pdbid=''):
        self.input = iter(input)
        self.pdbid = pdbid

    def __iter__(self):
        return self

    def next(self):
        while 1:
            l = self.input.next()
            for r in records:
                m = records[r]['re'].match(l)
                if m:
                    m = m.groupdict()
                    for k,v in m.items():
                        m[k] = records[r]['cn'][k](v)
                    m['pdbid'] = self.pdbid
                    return (r, m)

    def filter(self, records):
        return map(lambda (r,i): i, filter(lambda(r,i):r in records, self))

    def get_coords(self):
        """coordinates(pdb iterator) -> array

        >>> import StringIO
        >>> atomline  = "ATOM      1  N   ASP A   1      11.860  13.207  12.724  1.00 21.64           N \\n"
        >>> atomline += "ATOM      2  CA  ASP A   1      11.669  12.413  13.949  1.00 25.07           C \\n"
        >>> atomline += "ATOM      3  C   ASP A   1      11.824  13.276  15.220  1.00 21.99           C \\n"
        >>> atomline += "ATOM      4  O   ASP A   1      12.513  14.275  15.207  1.00 25.33           O \\n"
        >>> atomline += "ATOM      5  CB  ASP A   1      12.830  11.408  14.093  1.00 38.96           C \\n"
        >>> atomline += "ATOM   2227  N   VAL I  62      16.426  33.763  20.155  1.00 11.95          2384\\n"
        >>> atomfile = StringIO.StringIO(atomline)
        >>> p = reader(atomfile)
        >>> p.get_coords()
        array([[ 11.86 ,  13.207,  12.724,   1.   ],
               [ 11.669,  12.413,  13.949,   1.   ],
               [ 11.824,  13.276,  15.22 ,   1.   ],
               [ 12.513,  14.275,  15.207,   1.   ],
               [ 12.83 ,  11.408,  14.093,   1.   ],
               [ 16.426,  33.763,  20.155,   1.   ]])
        """
        l = []
        for i in self.filter('ATOM'):
            l.append([ i['x'], i['y'], i['z'], 1 ])
        return array(l)

class update_coords:
    """Class for update PQR data"""
    def __init__(self, pdb, coords):
        self.pdb = iter(pdb)
        self.coords = iter(coords)

    def __iter__(self):
        return self

    def next(self):
        r, a = self.pdb.next()
        if r == "ATOM":
            c = self.coords.next()
            a['x'], a['y'], a['z'] = c[0:3]
        return r, a

def _name_adaptor(name):
    if not name[-1].isdigit() and len(name) < 4:
        name = name + " "
    append_size = 4 - len(name)
    name = " "*append_size + name
    return name

class writter:
    """Class for iterate PDQ data"""
    def __init__(self, input):
        self.input = iter(input)

    def __iter__(self):
        return self

    def next(self):
        while 1:
            r, a = self.input.next()
            if r in records:
                a["name"] = _name_adaptor(a["name"])
                return records[r]['fm'] % a

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

