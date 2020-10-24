# -*- coding: utf-8 -*-

"""IO functions for trajectories."""

import re
import numpy as np
from .molecule import *

__all__ = ['Frame', 'read_cp2k']

header = re.compile(r' *i = *([0-9]+), time = *([-0-9.][^,]*), E = *([-0-9.][^ ]*)\n\Z')

class Frame:
    def __init__(self, index, energy, sys=None, x=None):
        self.index = index
        self.sys = sys
        self.x   = x
        self.energy = energy

    def validate(self):
        assert len(self.x) == self.sys.atoms()
        assert len(self.x.shape) == 2
        assert self.x.shape[1] == 3

def parse_sys(topol, names, crds):
    s = Sys()
    y = []
    for i, (p, m) in enumerate(topol): # permute atoms into molecules
        assert len(p) == m.atoms(), "Invalid permutation for molecule %d"%(i+1)
        x = np.empty((len(p), 3), np.float)
        for i,j in enumerate(p):
            x[i] = crds[j]
            assert m.names[i] == names[j]
        s.append(m)
        y.append(x)

    return s, np.vstack(y)

def read_cp2k(cp2k, topol):
    # Reading data from an xyz trajectory output by cp2k
    #  topol is a list of tuples, [ (index_array, mol) ].
    #    Note: be careful to write a 1-atom index array with a comma: (n,)
    #
    #  each index array has atom numbers needed to make a molecule
    #  for example, a water and an O_2 formatted as:
    #    HW, HW, OW, O1, O2
    #    
    #  would look like:
    #    topol = [ ((2,1,0), Mol(["O", "H", "H"], 0, 1)),
    #              ((3, 4), Mol(["O", "O"], 0, 3)) ]

    dataFrameList = []
    with open(cp2k) as file_in:
        for line in file_in:
            # start new frame
            m = header.match(line)
            if m:
                # process last frame
                if len(dataFrameList) > 0:
                    sys, x = parse_sys(topol, names, crds)
                    dataFrameList[-1].sys = sys
                    dataFrameList[-1].x   = x
                    dataFrameList[-1].validate()
                i = int(m[1])
                ener = float(m[3])
                dataFrameList.append(Frame(i, ener))
                names = []
                crds = []
            else:
                tok = line.strip().split()
                if len(tok) < 4: continue
                names.append(tok[0])
                crds.append(list(map(float, tok[1:4])))

    if len(dataFrameList) > 0:
        sys, x = parse_sys(topol, names, crds)
        dataFrameList[-1].sys = sys
        dataFrameList[-1].x   = x
        dataFrameList[-1].validate()

    return dataFrameList
