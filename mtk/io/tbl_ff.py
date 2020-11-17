# -*- coding: ISO-8859-1 -*-
# $Id: tbl_ff.py 62 2009-05-19 09:15:03Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Load tables into a dict

    >>> input = "a:int    b:int c:int\\n" \
                "# Que lindo dia\\n" \
                "1 2    3\\n" \
                "4    5 6\\n"
    >>> from StringIO import StringIO
    >>> r = reader(StringIO(input))
    >>> for i in r: print i
    {'a': 1, 'c': 3, 'b': 2}
    {'a': 4, 'c': 6, 'b': 5}
    >>> r = reader(StringIO(input))
    >>> index(list(r), 'a')
    {1: [0], 4: [1]}
    >>> info = [ { 'a': 1, 'b': 0.1, 'q': 9 }, \
                 { 'a': 4, 'b': 0.2, 'q': 8 }, \
                 { 'a': 1, 'b': 0.1, 'q': 7 }, \
                 { 'a': 4, 'b': 0.4, 'q': 6 }, \
                 ]
    >>> for i in Join(reader(StringIO(input)), 'a', info, 'a'): print i
    {'a': 1, 'q': 9, 'c': 3, 'b': 0.10000000000000001}
    {'a': 1, 'q': 9, 'c': 3, 'b': 0.10000000000000001}
    {'a': 4, 'q': 8, 'c': 6, 'b': 0.20000000000000001}
    {'a': 4, 'q': 8, 'c': 6, 'b': 0.20000000000000001}
    >>> for i in Join(info,'a',reader(StringIO(input)), 'a'): print i
    {'a': 1, 'q': 9, 'c': 3, 'b': 2}
    {'a': 4, 'q': 8, 'c': 6, 'b': 5}
    {'a': 1, 'q': 7, 'c': 3, 'b': 2}
    {'a': 4, 'q': 6, 'c': 6, 'b': 5}
    >>> for i in NaturalJoin(info, reader(StringIO(input))): print i
    {'a': 1, 'q': 9, 'c': 3, 'b': 2}
    {'a': 4, 'q': 8, 'c': 6, 'b': 5}
    {'a': 1, 'q': 7, 'c': 3, 'b': 2}
    {'a': 4, 'q': 6, 'c': 6, 'b': 5}
    >>> for i in NaturalJoin(reader(StringIO(input)), info): print i
    {'a': 1, 'q': 9, 'c': 3, 'b': 0.10000000000000001}
    {'a': 1, 'q': 9, 'c': 3, 'b': 0.10000000000000001}
    {'a': 4, 'q': 8, 'c': 6, 'b': 0.20000000000000001}
    {'a': 4, 'q': 8, 'c': 6, 'b': 0.20000000000000001}


    >>> IN = "pais, prefijo\\n" \
             "Argentina, 54\\n" \
             "EEUU, 1\\n"
    >>> import csv
    >>> import StringIO
    >>> INf = StringIO.StringIO(IN)
    >>> r = csv.DictReader(INf)
    >>> r.next()
    {'pais': 'Argentina', ' prefijo': ' 54'}
    >>> r.next()
    {'pais': 'EEUU', ' prefijo': ' 1'}
    >>> r.next()
    Traceback (most recent call last):
        ...
    StopIteration
"""
import re
import sqlite3

spaces = re.compile('\s*')
comments = re.compile('#.*$')

class reader:
    def __init__(self, input):
        self.input = iter(input)
        keys = [ k.split(':') for k in spaces.split(input.next()) ]
        self.keys = [ i[0] for i in keys if len(i)==2 ]
        self.types = [ eval(i[1]) for i in keys if len(i)==2 ]

    def __iter__(self):
        return self

    def next(self):
        while True:
            l = self.input.next()
            L = [ i for i in spaces.split(comments.sub('', l)) if i != '' ]
            if L != []:
                i = dict( [ (k, t(v)) for k,t,v in zip(self.keys, self.types,  L) ] )
                return i

def index(iter, key):
    """Index a table by a key"""
    idx = {}
    for dict in iter:
        kv = dict[key]
        if not kv in idx:
            idx[kv] = []
        idx[kv].append(iter.index(dict))
    return idx

class Join:
    """Join two tables"""

    def __init__(self, data_iter, data_key, tbl_list, tbl_key):
        self.data_iter = iter(data_iter)
        self.data_key = data_key
        self.tbl_list = list(tbl_list)
        self.tbl_idx = index(self.tbl_list, tbl_key)
        self.S = []

    def __iter__(self):
        return self

    def next(self):
        if len(self.S) == 0:
            item = self.data_iter.next()
            key = item[self.data_key]
            data_to_append = [ self.tbl_list[i].items() for i in  self.tbl_idx[key] ] 
            self.S = []
            for d in data_to_append:
                self.S.append( dict(item.items() + d) )
        _r = self.S[0]
        self.S = self.S[:-1]
        return dict(_r)

def NaturalJoin(tableA, tableB):
    tableA = list(tableA)
    tableB = list(tableB)
    keysA = tableA[0].keys()
    keysB = tableB[0].keys()
    key = None
    for k in keysA:
        if k in keysB:
            key = k
            break
    if key == None:
        raise RuntimeError
    return Join(tableA, key, tableB, key)

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

