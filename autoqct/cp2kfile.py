# -*- coding: utf-8 -*-

"""IO functions for trajectories."""

import re
import numpy as np
from .molecule import *

__all__ = ['Frame', 'read_cp2k']

header = re.compile(r' *i = *([0-9]+), time = *([-0-9.][^,]*), E = *([-0-9.][^ ]*)\n\Z')

class Frame:
    def __init__(self, index, sys, energy):
        self.index = index
        self.sys = sys
        self.energy = energy

def parse_sys(topol, names, crds):
    s = Sys()
    for (p, m) in topol: # permute atoms into molecules
        m = m.copy()
        m.x = np.zeros((len(p), 3), np.float)
        for i,j in enumerate(p):
            m.x[i] = crds[j]
        s.append(m)

    return s

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
    #    topol = [ ((2,1,0), Mol(["O", "H", "H"], np.zeros(3,3), 0, 1)),
    #              ((3, 4), Mol(["O", "O"], np.zeros((2,3)), 0, 3)) ]

    file_in = open(cp2k)

    dataFrameList = []
    for line in file_in:
        # start new frame
        m = header.match(line)
        if m:
            # process last frame
            if len(dataFrameList) > 0:
                dataFrameList[-1].sys = parse_sys(topol, names, crds)
            i = int(m[1])
            ener = float(m[3])
            dataFrameList.append(Frame(i, None, ener))
            names = []
            crds = []
        else:
            tok = line.strip().split()
            if len(tok) < 4: continue
            names.append(tok[0])
            crds.append(list(map(float, tok[1:4])))

    return dataFrameList
