from mtk.storage import Storage
import os.path

S = Storage()
S.loadpdb(os.path.join('mtk', 'data', 'test', 'obj02.pdb'))
S.loadprm(os.path.join('mtk', 'data', 'test', 'amber99.prm'))

query_biocorrection = lambda molecule : """
     SELECT replace
     FROM tinkerbiocorrection
     WHERE atom=%s.name
       AND residue=%s.resName
""" % (molecule, molecule)

query_moleculename = lambda molecule, ffbiotype : """
ifnull((%s), %s.name) = %s.name
""" % (query_biocorrection(molecule), molecule, ffbiotype)

query_countequals = lambda molecule, ffbiotype, resSeq : """
SELECT count(*)
FROM molecule as %s, ffbiotype as %s
WHERE M.resSeq=%s
  AND %s
GROUP BY B.description
ORDER BY count(*) DESC
LIMIT 1
""" % (molecule, ffbiotype, resSeq, query_moleculename(molecule, ffbiotype))



query = """
SELECT B.description as description, count(*)
FROM molecule as M, ffbiotype as B
WHERE M.resSeq=8
  AND M.name = B.name 
GROUP BY B.description
HAVING count(*) = ( SELECT count(*)
            FROM molecule as M, ffbiotype as B
            WHERE M.resSeq=8
              AND M.name = B.name
            GROUP BY B.description
            ORDER BY count(*) DESC
            LIMIT 1
        )
"""

query = """
	SELECT *
	FROM molecule, aminoacids, (SELECT B.description as description, M.resSeq as resId
			FROM molecule as M, ffbiotype as B
			WHERE M.resSeq=8
			  AND M.name = B.name 
			GROUP BY B.description
			HAVING count(*) = ( SELECT count(*)
						FROM molecule as M, ffbiotype as B
						WHERE M.resSeq=8
						  AND M.name = B.name
						GROUP BY B.description
						ORDER BY count(*) DESC
						LIMIT 1	)
			)
	WHERE description = longname
	  AND molecule.resSeq=8
	  AND molecule.resName=aminoacids.triname
"""

query = """
	SELECT *
	FROM molecule, aminoacids, (SELECT B.description as description, M.resSeq as resId
			FROM molecule as M, ffbiotype as B
			WHERE M.resSeq=8
			  AND M.name = B.name 
			GROUP BY B.description
			HAVING count(*) = ( SELECT count(*)
						FROM molecule as M, ffbiotype as B
						WHERE M.resSeq=8
						  AND M.name = B.name
						GROUP BY B.description
						ORDER BY count(*) DESC
						LIMIT 1	)

	WHERE description = longname
	  AND molecule.resSeq=8
	  AND molecule.resName=aminoacids.triname
"""

query = """
SELECT serial, resSeq
FROM molecule
  AND Mb.resSeq = %s
"""


query = """
SELECT Bb.biotypeid, Mb.serial, Mb.resSeq, AA.triname
FROM molecule as Mb, ffbiotype as Bb, aminoacids as AA
WHERE Mb.name = Bb.name
  AND Bb.description = AA.longname
  AND Mb.resname = AA.triname
GROUP BY Bb.description, Mb.resSeq
HAVING count(*) = ( SELECT count(*)
            FROM molecule as Ma, ffbiotype as Ba
            WHERE Ma.name = Ba.name
			  AND Ma.resSeq=Mb.resSeq
            GROUP BY Ba.description
            ORDER BY count(*) DESC
            LIMIT 1	)
ORDER BY Mb.serial
"""

q = S.do(query)

for i in q:
    print i
