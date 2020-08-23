#!/usr/bin/env python3

import numpy as np
from autoqct import *

water = Mol(["O", "H", "H"], 0)
chloride = Mol(["Cl"], -1)

class ChlorideWater(QCT):
    molecules = [
        [(6,), chloride],
        [(0,2,5), water],
        [(1,3,4), water]
    ]
    theory='wB97X-D'
    basis='aug-cc-pvdz'
    radii = []
    H = [[2,3], [5,6]]
    def clustered(self, x):
        # all waters have >= 1 H close to Cl-
        return all( any(dist(x[0], x[j]) < 2.3 for j in Hs)
                    for Hs in self.H
                  )

autoQCT(ChlorideWater)

