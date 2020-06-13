import os
from .cp2kfile import *
from .qm_input import *

__all__ = ['QCT', 'autoQCT']

class QCT:
    def __init__(self):
        if not hasattr(self, 'molecules'):
            raise NameError("Incomplete user code: 'molecules' is not defined.")

        try:
            self.trj = read_cp2k('cp2k_cluster.xyz', self.molecules)
        except FileNotFoundError:
            raise ValueError('cp2k_cluster.xyz not found.')

    def validate_input(self):
        print("Read trajectory with %d frames."%(len(self.trj)))
        print("First frame has %d molecules:"%(len(self.trj[0].sys)))
        print()
        # debug only
        name = self.__class__.__name__

        print(mk_com(self.trj[0].sys, '%s_N1'%name, 4, '#P B3LYP/cc-pVDZ OPT FREQ NOSYMMETRY'))
        return True

def autoQCT(Q, target='valid'):
    q = Q()
    q.validate_input()
