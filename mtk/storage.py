# -*- coding: ISO-8859-1 -*-
# $Id: storage.py 109 2009-07-28 23:18:52Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
""" Modulo de almacenamiento de mtk

>>> S = Storage()
"""

import sqlite3
import csv
import os.path
import sys
import sqlfunctions
from pkg_resources import resource_filename
from time import time
from sys import stderr

class AssignError (RuntimeError):
    def __init__(self, message, notassigned):
        self.message = message
        self.notassigned = notassigned

    def str():
        return self.message

def normalizeTuple(tuple):
    """
    Make a tuple unique searcheable in the storage

    >>> normalizeTuple([ 1, 2 ])
    [1, 2]
    >>> normalizeTuple([ 2, 1 ])
    [1, 2]
    >>> normalizeTuple([ 1, 2, 3 ])
    [1, 2, 3]
    >>> normalizeTuple([ 3, 2, 1 ])
    [1, 2, 3]
    >>> normalizeTuple([ 1, 2, 3, 4 ])
    [1, 2, 3, 4]
    >>> normalizeTuple([ 4, 3, 2, 1 ])
    [1, 2, 3, 4]
    >>> normalizeTuple([ 4, 3, 1, 2 ])
    [2, 1, 3, 4]
    """
    if tuple[0] > tuple[-1]:
        tuple = list(reversed(tuple))
    return tuple

class Singleton:
    class Singleton(object):
     """ A Pythonic Singleton """
     def __new__(cls, *args, **kwargs):
         if '_inst' not in vars(cls):
             cls._inst = object.__new__(cls, *args, **kwargs)
         return cls._inst

class Storage(Singleton):
    """
    Class with storage capacity of tables and molecules.

    >>> S = Storage(None)
    >>> q = S.do("SELECT dist3d(0,0,0,0,0,0)")
    >>> q.next()[0] == 0.0
    True
    """
    def __init__(self, path=resource_filename(__name__, 'data'),
                 archive=":memory:", profile=False):
        self._debug = { 'path': path, 'archive': archive }
        self._profile = profile
        if 'path' not in vars(self):
            self.__setup__(path, archive)

    def __setup__(self, path, archive):
        self.path = path
        self.conn = sqlite3.connect(archive)

        # Create functions
        F = sqlfunctions.find_functions()
        for f in F:
            len, sep, name = f.partition('_')
            self.conn.create_function(name, int(len), F[f])

        self.cursor = self.conn.cursor()

        # Create tables if not exists
        if path != None:
            self.structfilename = os.path.join(self.path, 'csv', 'tables.csv')
            structfile = open(self.structfilename)
            d = csv.DictReader(structfile)
            for row in d:
                self.do("CREATE TABLE IF NOT EXISTS %s ( %s )" %(row['name'],
                                                                 row['columns']))
                if row['defaultdata'] != '':
                    self.loadcsv(row['name'], os.path.join(self.path, 'csv',
                                                           row['defaultdata']))

            # Create indexes
            indexfile = open(os.path.join(self.path, 'csv', 'indexes.csv'))
            d = csv.DictReader(indexfile)
            for row in d:
                self.do("CREATE INDEX IF NOT EXISTS %s ON %s ( %s )" % 
                        (row['name'], row['table'], row['columns']))

    def loadcsv(self, table, file):
        """
        Carga un archivo csv en una tabla.
        """
        self._debug['file'] = file
        data = csv.DictReader(open(file))
        self.loaditr(table, data)
        del self._debug['file']
        return

    def loaditr(self, table, iter):
        """
        Carga de un iterador los valores en una tabla.
        """
        fields = None
        #self.begin()
        for i in iter:
            if fields == None: fields = ','.join(['?'] * len(i.keys()))
            self.loaddict(table, i, fields)
        #self.commit()
        return

    def loaddict(self, table, dict, fields=None):
        """
        Carga una linea descripta en un directorio en una tabla.
        """
        if fields == None: fields = ','.join(['?'] * len(dict.keys()))
        try:
            self.do("INSERT OR IGNORE INTO %s (%s) VALUES (%s)" % (table, ','.join(dict.keys()),
                                                    fields), *dict.values())
        except Exception, inst:
            raise RuntimeError(inst.message + '\n' + str(dict) + \
                               'file' in self._debug \
                                and '\n' + self._debug['file'] \
                                or '')
        return

    def loadprm(self, file, prmid=str()):
        """
        Load prm file into the storage.

        >>> import StringIO
        >>> prmline =  "atom      1     1    CT      \\"Amber CT\\"                  7     14.003     4\\n"
        >>> prmline += "atom      2     2    C       \\"Amber C\\"                   6     12.000     3\\n"
        >>> prmline += "atom      3     3    CA      \\"Amber CA\\"                  6     12.000     3\\n"
        >>> prmline += "atom      4     4    CM      \\"Amber CM\\"                  1      1.008     3\\n"
        >>> prmline += "atom      5     5    CC      \\"Amber CC\\"                  8     15.995     3\\n"
        >>> prmline += "biotype      1    CT      \\"Amber CT\\"                  1\\n"
        >>> prmline += "biotype      2    C       \\"Amber C\\"                   2\\n"
        >>> prmline += "biotype      3    CA      \\"Amber CA\\"                  3\\n"
        >>> prmfile = StringIO.StringIO(prmline)
        >>> S = Storage()
        >>> S.loadprm(prmfile)
        >>> d = S.do("SELECT atomid,name,number FROM ffatom")
        >>> d.next() == (1, u'CT', 7) and d.next() == (2, u'C', 6) 
        True
        >>> d = S.do("SELECT biotypeid,name,atomid FROM ffbiotype")
        >>> d.next() == (1, 'CT', 1) and d.next() == (2, 'C', 2) 
        True
        >>> S = Storage()
        >>> S.loadprm(os.path.join(S.path, 'test', 'amber99.prm'))
        >>> d = S.do("SELECT atomid,name,number FROM ffatom")
        >>> d.next() == (1, u'N', 7) and d.next() == (2, u'CT', 6) 
        True
        """
        import io.prm_ff as prmff

        if type(file) == str:
            file = open(file)
        prm = prmff.read(file)
        tbl = prmff.readtables(self.structfilename)

        for table, attribs in [ (t,a) for t,a in tbl.items() if t in prm ]:
            for key, data in prm[table]:
                for datum in data:
                    item = key + datum
                    self.loaddict("ff"+table, dict(zip(attribs, item)))

        # Complete tables
        self.do("""
             INSERT INTO ffresidues
             SELECT description, sort_group(group_concat(name, ','))
             FROM ffbiotype
             GROUP BY description
             """)

        return

    def loadpdb(self, file, pdbid=None):
        """
        Load pdb file into the storage

        WARN: No check if the pdbid yet exists.

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
        >>> S = Storage()
        >>> S.loadpdb(atomfile)
        ''
        >>> d = S.do("SELECT x,y,z FROM molecule")
        >>> d.next() == (11.859999999999999, 13.207000000000001, 12.724)
        True
        >>> d = S.do("SELECT serial,serialcon FROM connect")
        >>> d.next() == (1, 2) and d.next() == (1, 6)
        True

        >>> S = Storage()
        >>> S.loadpdb(os.path.join(S.path, 'test', '1PPE_l_u.pdb'))
        '1PPE_l_u.pdb'
        >>> d = S.do("SELECT x,y,z FROM molecule")
        >>> d.next() == (0.60899999999999999, 18.920000000000002, 11.647)
        True

        >>> S.loadpdb('1J2J')
        '1J2J'
        >>> d = S.do("SELECT serial,serialcon FROM connect WHERE pdbid='1J2J'")
        >>> d.next()
        (109, 1664)
        >>> d.next()
        (242, 1664)
        >>> d.next()
        (1664, 109)
        >>> d.next()
        (1664, 242)
        >>> d.next()
        (1664, 1668)

        """
        import io.pdb_ff as pdbff
        from os.path import basename
        from sys import stderr as err
        from StringIO import StringIO

        #print >> err, "[%f] loadpdb start" % time()

        if type(file) == str:
            if len(file) == 4:
                from urllib import urlretrieve
                from gzip import open as gzopen
                from zlib import decompress
                pdbid = file
                file, headers = urlretrieve("http://www.pdb.org/pdb/files/%s.pdb.gz" % pdbid)
                file = gzopen(file)
            else:
                if pdbid == None:
                    pdbid = basename(file)
                file = open(file)
        else:
            pdbid = str()

        l = file.readlines()

        self.loaditr("molecule", pdbff.reader(iter(l), pdbid=pdbid).filter(['ATOM', 'HETATM']))
        self.loaditr("connect", pdbff.connectorIterator(pdbff.reader(iter(l), pdbid=pdbid)))

        # Calculate other information
        self.do("""
                INSERT INTO residues
                SELECT pdbid, resSeq, sort_group(group_concat(name, ','))
                FROM molecule
                WHERE pdbid=="%s"
                GROUP BY resSeq
                """ % pdbid)

        #print >> err, "[%f] loadpdb end" % time()

        return pdbid

    def loadpqr(self, file, pdbid=''):
        """
        Load prq file into the storage

        >>> import StringIO
        >>> atomline  = "ATOM   3166  HB3 GLU   208     -25.296  -8.084  54.423 -0.0173 1.4870\\n"
        >>> atomline += "ATOM   3167  HB2 GLU   208     -24.501  -7.989  53.000 -0.0173 1.4870\\n"
        >>> atomline += "ATOM   3168  N   CYS   209     -27.344  -5.648  53.887 -0.4157 1.8240\\n"
        >>> atomline += "ATOM   3169  CA  CYS   209     -27.274  -4.213  54.184  0.0429 1.9080\\n"
        >>> atomfile = StringIO.StringIO(atomline)
        >>> S = Storage()
        >>> S.loadpqr(atomfile)
        ''
        >>> d = S.do("SELECT x,y,z,charge,occupancy FROM molecule")
        >>> d.next() == (-25.295999999999999, -8.0839999999999996, 54.423000000000002, -0.0173, 1.4870000000000001)
        True

        """
        import io.pqr_ff as pqrff
        from os.path import basename
        from sys import stderr as err
        from StringIO import StringIO

        if type(file) == str:
            file = open(file)

        l = file.readlines()

        self.loaditr("molecule", pqrff.reader(iter(l), pdbid=pdbid).filter(['ATOM']))
        return pdbid

    def begin(self):
        self.cursor.execute('BEGIN;')

    def commit(self):
        self.cursor.execute('COMMIT;')

    def do(self, sql, *args):
        """
        Ejecuta un comando sql.
        """
        if self._profile:
            t = time()
            print >> stderr, "starting SQL execution time:", t

        try:
            if type(sql) != list:
                sql = [ sql ]
            for command in sql:
                self.cursor.execute(command, args)
        except sqlite3.OperationalError, s:
            raise RuntimeError(s.message + '\n' + command + '\n' + str(args))

        if self._profile:
            print >> stderr, "SQL time:", time()-t

        return iter(self.cursor)

    def get_coords(self, pdbid=None):
        """
        Obtiene las coordinadas de una proteina en el storage.
        """
        from numpy import array
        if pdbid != None:
            d = self.do("SELECT x,y,z,1.0 FROM molecule WHERE pdbid==?", pdbid)
        else:
            d = self.do("SELECT x,y,z,1.0 FROM molecule")
        return array(list(d))

    def get_charges(self, pdbid=None):
        """
        Obtiene las cargas de proteina.
        """
        # TODO: Tiene que devolver las cargas con el mismo orden que las coordenadas
        from numpy import array
        if pdbid != None:
            d = self.do("SELECT charge FROM molecule, charges WHERE pdbid==?", pdbid)
        else:
            d = self.do("SELECT charge FROM molecule, charges WHERE ")
        return array(list(d))

    def get_atoms(self, atts=['serial', 'x','y','z', 'ffcharge.charge', 'radix', 'epsilon', 'resName']):
        '''
        Get full atom information.

        >>> S = Storage()

        > >> S.loadpdb(os.path.join(S.path, 'test', 'obj02.pdb'))
             'obj02.pdb'
        >>> S.loadpdb(os.path.join(S.path, 'test', '1PPE_l_u.pdb'))
        '1PPE_l_u.pdb'
        >>> S.loadprm(os.path.join(S.path, 'test', 'amber99.prm'))
        >>> len( S.get_atoms() ) == len( S.get_coords() )
        True
        '''

        query = ["""
        CREATE TEMPORARY TABLE IF NOT EXISTS Tresidues
        (resSerial INTEGER, triname TEXT, atoms TEXT)
        """, """
        DELETE FROM Tresidues
        """, """
        CREATE TEMPORARY TABLE IF NOT EXISTS Tdescripts
        (description TEXT, triname TEXT, atoms TEXT)
        """, """
        DELETE FROM Tdescripts
        """, """
        CREATE TEMPORARY TABLE IF NOT EXISTS Taa
        (TaaSerial INTEGER, TaaName TEXT)
        """, """
        DELETE FROM Taa
        """, """
        INSERT INTO Tresidues
        SELECT resSeq, resName, group_concat(
                 ifnull((SELECT replace FROM tinkerbiocorrection WHERE atom=molecule.name AND residue=molecule.resName),
                 molecule.name) , ',')
        FROM molecule
        GROUP BY resSeq
        """, """
        INSERT INTO Tdescripts
        SELECT description, triname, group_concat(ffbiotype.name, ',')
        FROM ffbiotype, aminoacids
        WHERE description = aminoacids.longname
        GROUP BY description
        """, """
        INSERT INTO Taa
        SELECT TR.resSerial, TD.description
        FROM Tresidues as TR, Tdescripts as TD
        WHERE eqitems(TR.atoms, TD.atoms) =
                 (SELECT eqitems(Tresidues.atoms, Tdescripts.atoms) as eqv
                    FROM Tresidues, Tdescripts
                    WHERE TR.resSerial = Tresidues.resSerial
                    AND Tresidues.triname = Tdescripts.triname
                    ORDER BY eqv DESC
                    LIMIT 1)
        AND TD.triname = TR.triname
        GROUP BY TR.resSerial
        """, """
        SELECT serial, %s
        FROM molecule, Taa, ffbiotype, ffatom, ffcharge, ffvdw
        WHERE ifnull((SELECT replace FROM tinkerbiocorrection
            WHERE atom=molecule.name AND residue=molecule.resName),
                molecule.name) = ffbiotype.name
              AND molecule.resSeq = Taa.TaaSerial
              AND ffbiotype.description = Taa.TaaName
              AND ffbiotype.atomid = ffatom.atomid
              AND ffbiotype.atomid = ffcharge.atomid
              AND ffatom.classid = ffvdw.classid
        """ % (','.join(atts))]

        d = list(self.do(query))
        c = list(self.do("SELECT serial, name, resname FROM molecule"))

        if len(d) != len(c):
            di = [ i[0] for i in d ]
            ci = [ i[0] for i in c ]
            mask = [ i for i in ci if not i in di ]
            if len(mask) == 0:
                import pdb; pdb.set_trace()
                print ci, di
                raise RuntimeError("Check for duplicated entries of atoms in the forcefield.")
            else:
                d = [ (serial, name, resname) for serial, name, resname in c if serial in mask ]
                raise AssignError("Trouble to assign atoms.", d)

        return list(d)

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

