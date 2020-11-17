
# parsetab.py
# This file is automatically generated. Do not edit.
_tabversion = '3.2'

_lr_method = 'LALR'

_lr_signature = '\xe2\xd4J\x82\x96-\x8e\xdf\x9a\xa2\x04Xa\xadj\xa2'
    
_lr_action_items = {'REAL':([37,38,53,57,58,60,61,62,65,66,68,72,73,74,79,80,81,82,83,86,87,88,89,90,91,92,94,95,99,100,101,102,103,104,105,106,107,108,110,114,117,119,120,121,122,123,125,126,127,130,131,132,133,134,135,136,138,139,140,141,143,144,145,146,148,150,151,152,153,154,155,157,158,159,161,162,163,166,167,168,173,191,192,194,195,197,203,204,207,208,209,212,213,215,219,222,],[97,98,113,116,117,119,-58,120,121,122,124,125,126,127,130,133,134,135,136,138,139,140,141,142,143,144,130,130,148,149,150,151,152,153,154,155,156,157,158,130,161,163,164,165,166,-59,167,168,169,172,-64,130,174,175,176,177,179,180,181,182,183,184,130,130,186,187,188,189,190,191,192,193,194,130,196,-60,197,199,200,201,-65,207,208,209,-61,210,130,-67,216,217,218,-64,219,221,222,224,]),'$end':([54,64,71,76,84,],[0,-4,-2,-3,-1,]),'OPBEND':([0,],[1,]),'ANGLEF':([0,],[2,]),'ANGLE':([0,],[3,]),'ANGANG':([0,],[14,]),'RADIUSTYPE':([0,],[5,]),'EOL':([0,6,12,17,24,63,69,70,75,93,97,98,109,111,112,113,116,124,131,132,142,145,146,149,156,159,163,164,165,169,173,174,175,176,177,179,180,181,182,183,184,186,187,188,189,190,193,195,196,198,199,200,201,204,205,210,212,214,216,217,218,220,223,224,],[-5,64,71,76,84,-9,-7,-11,-46,-52,-12,-14,-6,-10,-8,-13,-36,-39,-64,-53,-26,-56,-55,-37,-41,-54,-28,-17,-47,-21,-65,-22,-44,-43,-45,-24,-23,-25,-34,-20,-19,-38,-50,-51,-33,-42,-49,-61,-32,-18,-48,-40,-35,-67,-15,-27,-63,-57,-29,-30,-31,-66,-16,-68,]),'IMPROPER':([0,],[39,]),'POLARIZE':([0,],[7,]),'PIATOM':([0,],[8,]),'PITORS':([0,],[9,]),'EPSILONRULE':([0,],[11,]),'COMMENT':([0,],[12,]),'VDW':([0,],[4,]),'HBOND':([0,],[15,]),'MULTIPOLE':([0,],[16,]),'CHARGE':([0,],[47,]),'IMPTORS':([0,],[19,]),'DIPOLE4':([0,],[21,]),'DIPOLE5':([0,],[22,]),'DIPOLE3':([0,],[23,]),'INTEGER':([1,2,3,4,7,8,9,13,14,15,16,18,19,20,21,22,23,25,26,27,28,29,30,31,32,34,35,36,39,40,41,42,43,44,45,46,47,48,50,55,56,59,61,67,77,78,96,115,118,128,129,160,164,170,171,172,178,202,206,211,221,222,],[56,59,59,61,61,61,67,56,61,67,56,77,56,67,67,67,67,61,67,67,67,59,59,67,61,56,56,61,56,56,67,67,59,67,59,59,61,67,59,56,115,118,-58,123,128,129,147,160,162,170,171,195,198,202,203,204,205,211,215,-62,223,204,]),'BOND':([0,],[20,]),'BOND4':([0,],[26,]),'BOND5':([0,],[27,]),'BOND3':([0,],[28,]),'UREYBRAD':([0,],[29,]),'ELECTNEG':([0,],[30,]),'VDWPR':([0,],[31,]),'VDW14':([0,],[32,]),'METAL':([0,],[33,]),'TORSION4':([0,],[34,]),'TORSION5':([0,],[35,]),'ATOM':([0,],[36,]),'VDW_14_SCALE':([0,],[37,]),'VDWTYPE':([0,],[10,]),'ID':([5,10,11,49,51,52,61,85,147,],[63,69,70,109,111,112,-58,137,185,]),'DIELECTRIC':([0,],[38,]),'DESC':([33,137,185,],[93,178,206,]),'CHG_14_SCALE':([0,],[53,]),'PIBOND5':([0,],[41,]),'PIBOND4':([0,],[42,]),'STRBND':([0,],[43,]),'DIPOLE':([0,],[44,]),'ANGLE5':([0,],[45,]),'ANGLE4':([0,],[46,]),'TORTORS':([0,],[18,]),'PIBOND':([0,],[48,]),'FORCEFIELD':([0,],[49,]),'ANGLE3':([0,],[50,]),'RADIUSSIZE':([0,],[51,]),'BIOTYPE':([0,],[25,]),'RADIUSRULE':([0,],[52,]),'OPDIST':([0,],[40,]),'STRTORS':([0,],[13,]),'TORSION':([0,],[55,]),}

_lr_action = { }
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = { }
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'torpar':([79,94,95,114,132,145,146,159,203,213,],[131,131,131,131,173,173,173,173,212,173,]),'torparams':([79,94,95,114,203,],[132,145,146,159,213,]),'tortorparams':([203,],[214,]),'simplevalue':([0,],[17,]),'multiplevalues':([0,],[6,]),'tortorpar':([213,],[220,]),'idx5':([18,],[78,]),'idx3':([2,3,29,30,43,45,46,50,],[58,60,89,90,103,105,106,110,]),'idx2':([9,15,20,21,22,23,26,27,28,31,41,42,44,48,],[68,74,80,81,82,83,86,87,88,91,101,102,104,108,]),'idx1':([4,7,8,14,25,32,36,47,],[62,65,66,73,85,92,96,107,]),'expression':([0,],[54,]),'empty':([0,],[24,]),'idx4':([1,13,16,19,34,35,39,40,55,],[57,72,75,79,94,95,99,100,114,]),}

_lr_goto = { }
for _k, _v in _lr_goto_items.items():
   for _x,_y in zip(_v[0],_v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = { }
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> expression","S'",1,None,None,None),
  ('expression -> empty EOL','expression',2,'p_expression_none','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',116),
  ('expression -> COMMENT EOL','expression',2,'p_expression_none','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',117),
  ('expression -> simplevalue EOL','expression',2,'p_expression','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',123),
  ('expression -> multiplevalues EOL','expression',2,'p_expression','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',124),
  ('empty -> <empty>','empty',0,'p_empty','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',130),
  ('simplevalue -> FORCEFIELD ID','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',135),
  ('simplevalue -> VDWTYPE ID','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',136),
  ('simplevalue -> RADIUSRULE ID','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',137),
  ('simplevalue -> RADIUSTYPE ID','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',138),
  ('simplevalue -> RADIUSSIZE ID','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',139),
  ('simplevalue -> EPSILONRULE ID','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',140),
  ('simplevalue -> VDW_14_SCALE REAL','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',141),
  ('simplevalue -> CHG_14_SCALE REAL','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',142),
  ('simplevalue -> DIELECTRIC REAL','simplevalue',2,'p_simplevalue','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',143),
  ('multiplevalues -> BIOTYPE idx1 ID DESC INTEGER','multiplevalues',5,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',149),
  ('multiplevalues -> ATOM idx1 INTEGER ID DESC INTEGER REAL INTEGER','multiplevalues',8,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',150),
  ('multiplevalues -> VDW idx1 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',151),
  ('multiplevalues -> VDW idx1 REAL REAL INTEGER','multiplevalues',5,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',152),
  ('multiplevalues -> VDW14 idx1 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',153),
  ('multiplevalues -> VDWPR idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',154),
  ('multiplevalues -> HBOND idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',155),
  ('multiplevalues -> BOND idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',156),
  ('multiplevalues -> BOND5 idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',157),
  ('multiplevalues -> BOND4 idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',158),
  ('multiplevalues -> BOND3 idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',159),
  ('multiplevalues -> ELECTNEG idx3 REAL','multiplevalues',3,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',160),
  ('multiplevalues -> ANGLE idx3 REAL REAL REAL REAL','multiplevalues',6,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',161),
  ('multiplevalues -> ANGLE idx3 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',162),
  ('multiplevalues -> ANGLE5 idx3 REAL REAL REAL REAL','multiplevalues',6,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',163),
  ('multiplevalues -> ANGLE4 idx3 REAL REAL REAL REAL','multiplevalues',6,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',164),
  ('multiplevalues -> ANGLE3 idx3 REAL REAL REAL REAL','multiplevalues',6,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',165),
  ('multiplevalues -> ANGLEF idx3 REAL REAL REAL','multiplevalues',5,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',166),
  ('multiplevalues -> STRBND idx3 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',167),
  ('multiplevalues -> UREYBRAD idx3 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',168),
  ('multiplevalues -> ANGANG idx1 REAL REAL REAL','multiplevalues',5,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',169),
  ('multiplevalues -> OPBEND idx4 REAL','multiplevalues',3,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',170),
  ('multiplevalues -> OPDIST idx4 REAL','multiplevalues',3,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',171),
  ('multiplevalues -> IMPROPER idx4 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',172),
  ('multiplevalues -> PITORS idx2 REAL','multiplevalues',3,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',173),
  ('multiplevalues -> STRTORS idx4 REAL REAL REAL','multiplevalues',5,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',174),
  ('multiplevalues -> CHARGE idx1 REAL','multiplevalues',3,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',175),
  ('multiplevalues -> DIPOLE idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',176),
  ('multiplevalues -> DIPOLE5 idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',177),
  ('multiplevalues -> DIPOLE4 idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',178),
  ('multiplevalues -> DIPOLE3 idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',179),
  ('multiplevalues -> MULTIPOLE idx4','multiplevalues',2,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',180),
  ('multiplevalues -> POLARIZE idx1 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',181),
  ('multiplevalues -> PIATOM idx1 REAL REAL REAL','multiplevalues',5,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',182),
  ('multiplevalues -> PIBOND idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',183),
  ('multiplevalues -> PIBOND5 idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',184),
  ('multiplevalues -> PIBOND4 idx2 REAL REAL','multiplevalues',4,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',185),
  ('multiplevalues -> METAL DESC','multiplevalues',2,'p_multiplevalues','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',186),
  ('multiplevalues -> IMPTORS idx4 torparams','multiplevalues',3,'p_multiplevalues_torparams','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',192),
  ('multiplevalues -> TORSION idx4 torparams','multiplevalues',3,'p_multiplevalues_torparams','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',193),
  ('multiplevalues -> TORSION5 idx4 torparams','multiplevalues',3,'p_multiplevalues_torparams','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',194),
  ('multiplevalues -> TORSION4 idx4 torparams','multiplevalues',3,'p_multiplevalues_torparams','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',195),
  ('multiplevalues -> TORTORS idx5 INTEGER INTEGER INTEGER tortorparams','multiplevalues',6,'p_multiplevalues_params','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',201),
  ('idx1 -> INTEGER','idx1',1,'p_idx1','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',207),
  ('idx2 -> INTEGER INTEGER','idx2',2,'p_idx1','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',208),
  ('idx3 -> INTEGER INTEGER INTEGER','idx3',3,'p_idx1','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',209),
  ('idx4 -> INTEGER INTEGER INTEGER INTEGER','idx4',4,'p_idx1','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',210),
  ('idx5 -> INTEGER INTEGER INTEGER INTEGER INTEGER','idx5',5,'p_idx1','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',211),
  ('tortorparams -> torpar','tortorparams',1,'p_tortorparams_item','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',217),
  ('torparams -> torpar','torparams',1,'p_tortorparams_item','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',218),
  ('torparams -> torparams torpar','torparams',2,'p_torparams_list','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',224),
  ('tortorparams -> torparams tortorpar','tortorparams',2,'p_torparams_list','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',225),
  ('torpar -> REAL REAL INTEGER','torpar',3,'p_torpar','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',231),
  ('tortorpar -> REAL REAL REAL','tortorpar',3,'p_torpar','/home/crocha/Proyectos/crocha_doctorado/apps/python-mtk/trunk/mtk/io/prm_ff.py',232),
]
