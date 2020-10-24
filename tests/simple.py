import unittest
from autoqct import *
import numpy as np

class TestSlurm(unittest.TestCase):
    def setUp(self):
        self.wat = Mol(["O", "H", "H"], np.identity(3), 0, 1)
        self.cl  = Mol(["Cl"], np.zeros((1,3)), -1)

    def test_molecule(self):
        self.assertEqual(0, self.wat.tot_charge(), "Incorrect Mol.tot_charge()")
        self.assertEqual(1, self.wat.tot_spin(), "Incorrect Mol.tot_spin()")

    def test_com(self):
        self.assertEqual(mk_com(self.wat, "wat_en", 16, "#N B3LYP/aug-cc-pvdz SP NOSYMMETRY"), '%chk=H.chk\n%NProcShared=16\n#N B3LYP/aug-cc-pvdz SP NOSYMMETRY\n\nH\n\n0 1\n   O      1.000000     0.000000     0.000000\n   H      0.000000     1.000000     0.000000\n   H      0.000000     0.000000     1.000000\n')

    def test_slurm(self):
        slurm = mk_slurm(["1.com", "2.com", "3.com"],partition='defq',qos='normal',time='1:00:00',nodes='1',ntasks_per_node='1',cpus_per_task='20',gres='mic:0',mem='64000',export='ALL')
        self.assertEqual(slurm, '#!/bin/bash\n#SBATCH -J BatchGaussian\n#SBATCH --partition=defq\n#SBATCH --qos=normal\n#SBATCH --time=1:00:00\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=20\n#SBATCH --gres=mic:0\n#SBATCH --mem=64000\n#SBATCH --export=ALL\n\nmodule load gaussian09\n\nmpirun --bynode -np 1 g09 1.com : g09 2.com : g09 3.com\n')

    # failing
    def test_io(self):
        topol = [
            [(6,), self.cl],
            [(0,2,3), self.wat],
            [(1,4,5), self.wat]
        ]
        x = read_cp2k("examples/chloride/cp2k_cluster.xyz", topol)
        self.assertEqual(701, len(x))

if __name__ == '__main__':
    unittest.main()
