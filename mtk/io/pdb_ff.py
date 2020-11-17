# -*- coding: ISO-8859-1 -*-
# $Id: pdb_ff.py 81 2009-06-02 18:55:20Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""PDB file format processor module.

Reading a pdb fileformat from a string as a stream. The iterator of
PDB entries are stored in p as a dict.

Read standard input format

>>> import StringIO
>>> atomline  = "REMARK FILENAME=\\"/usr/people/nonella/xplor/benchmark1/ALANIN.PDB\\"\\n"
>>> atomline += "ATOM      1  CA  ACE     1      -2.184   0.591   0.910  1.00  7.00      \\n"
>>> atomfile = StringIO.StringIO(atomline)
>>> p = reader(atomfile)
>>> p.next() == ('REMARK', {'comment': 'FILENAME="/usr/people/nonella/xplor/benchmark1/ALANIN.PDB"',\
                            'record': 'REMARK', 'pdbid': ''})
True
>>> p.next() == ('ATOM', {'comment': 'None', 'tempFactor': 7.0, \
                          'resName': 'ACE', 'chainID': '', 'name': 'CA', 'record': 'ATOM',\
                          'altLoc': '', 'occupancy': 1.0, 'element': 'None', 'resSeq': 1,\
                          'charge': '', 'iCode': '', 'pdbid': '', 'y': 0.59099999999999997,\
                          'x': -2.1840000000000002, 'serial': 1, 'z': 0.91000000000000003})
True

Read standard input alternative

>>> atomline  = "ATOM      1  N   ASP A   1      11.860  13.207  12.724  1.00 21.64           N \\n"
>>> atomline += "ATOM      2  CA  ASP A   1      11.669  12.413  13.949  1.00 25.07           C \\n"
>>> atomline += "ATOM      3  C   ASP A   1      11.824  13.276  15.220  1.00 21.99           C \\n"
>>> atomline += "ATOM      4  O   ASP A   1      12.513  14.275  15.207  1.00 25.33           O \\n"
>>> atomline += "ATOM      5  CB  ASP A   1      12.830  11.408  14.093  1.00 38.96           C \\n"
>>> atomline += "ATOM   2227  N   VAL I  62      16.426  33.763  20.155  1.00 11.95      1ACB2384\\n"
>>> atomfile = StringIO.StringIO(atomline)
>>> p = reader(atomfile)
>>> p.next() == ('ATOM', {'comment': '', 'tempFactor': 21.640000000000001, \
    'resName': 'ASP', 'chainID': 'A', 'name': 'N', 'altLoc': '', 'occupancy': \
    1.0, 'y': 13.207000000000001, 'resSeq': 1, 'charge': '', 'iCode': '', \
    'x': 11.859999999999999, 'serial': 1, 'z': 12.724, 'element': 'N', \
    'pdbid': '', 'record': 'ATOM'})
True

Write the input to the standard output without any modification

>>> import sys
>>> atomfile.seek(0)                   # Whe need restart the file
>>> sys.stdout.writelines(writter(p))
ATOM      1  N   ASP A   1      11.860  13.207  12.724  1.00 21.64           N  
ATOM      2  CA  ASP A   1      11.669  12.413  13.949  1.00 25.07           C  
ATOM      3  C   ASP A   1      11.824  13.276  15.220  1.00 21.99           C  
ATOM      4  O   ASP A   1      12.513  14.275  15.207  1.00 25.33           O  
ATOM      5  CB  ASP A   1      12.830  11.408  14.093  1.00 38.96           C  
ATOM   2227  N   VAL I  62      16.426  33.763  20.155  1.00 11.95          2384

Get the coords as an array of affine vectors

>>> atomfile.seek(0)                   # Whe need restart the file
>>> c = p.get_coords()
>>> c
array([[ 11.86 ,  13.207,  12.724,   1.   ],
       [ 11.669,  12.413,  13.949,   1.   ],
       [ 11.824,  13.276,  15.22 ,   1.   ],
       [ 12.513,  14.275,  15.207,   1.   ],
       [ 12.83 ,  11.408,  14.093,   1.   ],
       [ 16.426,  33.763,  20.155,   1.   ]])

Write modified protein coordinates to the standard output.

>>> from numpy import ones
>>> i = ones(c.shape)
>>> atomfile.seek(0)                   # Whe need restart the file
>>> q = update_coords(p, i)
>>> sys.stdout.writelines(writter(q))
ATOM      1  N   ASP A   1       1.000   1.000   1.000  1.00 21.64           N  
ATOM      2  CA  ASP A   1       1.000   1.000   1.000  1.00 25.07           C  
ATOM      3  C   ASP A   1       1.000   1.000   1.000  1.00 21.99           C  
ATOM      4  O   ASP A   1       1.000   1.000   1.000  1.00 25.33           O  
ATOM      5  CB  ASP A   1       1.000   1.000   1.000  1.00 38.96           C  
ATOM   2227  N   VAL I  62       1.000   1.000   1.000  1.00 11.95          2384
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
            r".{3}"
            r"(?P<x>[\-\+\d\. ]{8})"
            r"(?P<y>[\-\+\d\. ]{8})"
            r"(?P<z>[\-\+\d\. ]{8})"
            r"(?P<occupancy>[\-\+\d\. ]{6})"
            r"(?P<tempFactor>[\-\+\d\. ]{6})"
            r"(?:.{10}"
            r"(?P<element>[\w ]{2})"
            r"(?P<charge>[\w ]{2})?"
            r"(?P<comment>.*)|\s*\w*\s*)"
            r"$"
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
            "occupancy": lambda x: float(str(x).strip()),
            "tempFactor": lambda x: float(str(x).strip()),
            "element": lambda x: str(x).strip(),
            "charge": lambda x: x != None and str(x).strip() or "",
            "comment": lambda x: str(x).strip(),
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
            "%(occupancy)6.2f"
            "%(tempFactor)6.2f"
            "          "
            "%(element)2s"
            "%(charge)2s"
            "%(comment)s"
            "\n"
        },
    "HETATM": {
        "re": re.compile(
            r"^"
            r"(?P<record>HETATM)"
            r"(?P<serial>[\d ]{5})"
            r".{1}"
            r"(?P<name>[\w ]{4})"
            r"(?P<altLoc>[\w ]{1})"
            r"(?P<resName>[\w ]{3})"
            r".{1}"
            r"(?P<chainID>[\w ]{1})"
            r"(?P<resSeq>.{4})"
            r"(?P<iCode>.{1})"
            r".{3}"
            r"(?P<x>[\-\+\d\. ]{8})"
            r"(?P<y>[\-\+\d\. ]{8})"
            r"(?P<z>[\-\+\d\. ]{8})"
            r"(?P<occupancy>[\-\+\d\. ]{6})"
            r"(?P<tempFactor>[\-\+\d\. ]{6})"
            r"(?:.{10}"
            r"(?P<element>[\w ]{2})"
            r"(?P<charge>[\w ]{2})?"
            r"(?P<comment>.*)|\s*\w*\s*)"
            r"$"
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
            "occupancy": lambda x: float(str(x).strip()),
            "tempFactor": lambda x: float(str(x).strip()),
            "element": lambda x: str(x).strip(),
            "charge": lambda x: x != None and str(x).strip() or "",
            "comment": lambda x: str(x).strip(),
            },
        "fm": "HETATM"
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
            "%(occupancy)6.2f"
            "%(tempFactor)6.2f"
            "          "
            "%(element)2s"
            "%(charge)2s"
            "%(comment)s"
            "\n"
        },
    "CONECT": {
        "re": re.compile(
            r"^"
            r"(?P<record>CONECT)"
            r"(?P<serial>[\d ]{5})"
            r"(?P<serialcon0>[\d ]{5})"
            r"(?P<serialcon1>[\d ]{5})?"
            r"(?P<serialcon2>[\d ]{5})?"
            r"(?P<serialcon3>[\d ]{5})?"
            r".*"
            r"$"
            ),
        "cn": { "record": lambda x: str(x).strip(),
            "serial": lambda x: int(str(x).strip()),
            "serialcon0": lambda x: int(str(x).strip()),
            "serialcon1": lambda x: int(str(x).strip()) if not x == ' '*5  else None,
            "serialcon2": lambda x: int(str(x).strip()) if not x == ' '*5  else None,
            "serialcon3": lambda x: int(str(x).strip()) if not x == ' '*5  else None,
            },
        "fm": "%(record)6s"
            "%(serial)5i"
            "%(serialcon0)5i"
            "%(serialcon1)5i"
            "%(serialcon2)5i"
            "%(serialcon3)5i"
            "\n"
        },
    }

class reader:
    """Class for iterate over a PDB input

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
    """Class for update PDB data"""
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

class writter:
    """Class for iterate PDB data"""
    def __init__(self, input):
        self.input = iter(input)

    def __iter__(self):
        return self

    def next(self):
        while 1:
            r, a = self.input.next()
            if r in records:
                a["name"] = len(a["name"]) != 3 and " "+a["name"] or a["name"] 
                return records[r]['fm'] % a

def download(pdbid, pdburl="http://www.rcsb.org/pdb/files/%s.pdb"):
    """Download protein from pdb.org using the id

    >>> X = list( download('1p5f') )
    >>> X = [ i for i in X if i[0] == 'ATOM' ]
    >>> X[0] == ('ATOM', {'comment': '', 'tempFactor': 37.43, 'resName': \
         'SER', 'chainID': 'A', 'name': 'N', 'altLoc': '', 'occupancy': \
            1.0, 'element': 'N', 'resSeq': 3, 'charge': '', 'iCode': '', \
            'y': 31.914999999999999, 'x': 9.1549999999999994, 'serial': 1,\
            'z': 16.690000000000001, 'record': 'ATOM', 'pdbid': '1p5f'})
    True
    """

    from urllib import urlopen
    X = urlopen(pdburl % pdbid)
    return reader(X, pdbid=pdbid)

class connectorIterator:
    """Return a list of connectors iterators

    >>> import StringIO
    >>> atomline  = "ATOM      1  N   ASP A   1      11.860  13.207  12.724  1.00 21.64           N \\n"
    >>> atomline += "ATOM      2  CA  ASP A   1      11.669  12.413  13.949  1.00 25.07           C \\n"
    >>> atomline += "ATOM      3  C   ASP A   1      11.824  13.276  15.220  1.00 21.99           C \\n"
    >>> atomline += "ATOM      4  O   ASP A   1      12.513  14.275  15.207  1.00 25.33           O \\n"
    >>> atomline += "ATOM      5  CB  ASP A   1      12.830  11.408  14.093  1.00 38.96           C \\n"
    >>> atomline += "ATOM   2227  N   VAL I  62      16.426  33.763  20.155  1.00 11.95          2384\\n"
    >>> atomline += "CONECT    1    2    6   22   33                                                \\n"
    >>> atomline += "CONECT    2    1    3    0    0                                                \\n"
    >>> atomline += "CONECT    3    2    4    8   34                                                \\n"
    >>> atomline += "CONECT    4    3    5    6   36                                                \\n"
    >>> atomline += "CONECT    5    4   40    0    0                                                \\n"
    >>> atomfile = StringIO.StringIO(atomline)
    >>> p = connectorIterator(reader(atomfile))
    >>> list(p)[0:3] == [{'pdbid': '', 'serialcon': 2, 'serial': 1}, {'pdbid': \
        '', 'serialcon': 6, 'serial': 1}, {'pdbid': '', 'serialcon': 22, \
        'serial': 1}]
    True
    """
    def __init__(self, pdbiter):
        self.pdbid = pdbiter.pdbid
        self.iter = iter(filter(lambda(r,i):r == 'CONECT', pdbiter))
        self.a = None
        self.C = []

    def __iter__(self):
        return self

    def next(self):
        if self.C == []:
            r, d = self.iter.next()
            self.a, self.C = int(d['serial']), \
                    [ int(d['serialcon%i' % i]) for i in xrange(0,4)
                      if not d['serialcon%i' % i] == None ]
        c = self.C[0]
        self.C = self.C[1:]
        return { "pdbid": self.pdbid, "serial": self.a, "serialcon": c }

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

