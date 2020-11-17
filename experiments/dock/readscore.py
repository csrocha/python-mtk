
import os.path
import os
import sys
from mtk.dock import load as loaddock
from numpy import array, max, argmax

SS = []
for root, dirs, files in os.walk('.'):
    for fn in files:
        if '.csv' in fn:
            R, T, S = loaddock(fn)
            SS += S
            break
SS = array(SS)

for i in SS:
    print i.real

