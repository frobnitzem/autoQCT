#!/usr/bin/env python3

import numpy as np
from autoqct import *

water = Mol(["O", "H", "H"], np.zeros((3,3)), 0)
chloride = Mol(["Cl"], np.zeros((1,3)), -1)

class ChlorideWater(QCT):
    molecules = [
        [(6,), chloride],
        [(0,2,3), water],
        [(1,4,5), water]
    ]
    theory='wB97X-D'
    basis='aug-cc-pvdz'
    def clustered(self, x):
        return all( dist(x[0], x[j]) < 2.3
                    for j in [2,3, 5,6]
                  )

autoQCT(ChlorideWater)

