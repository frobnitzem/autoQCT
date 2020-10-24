import unittest
from autoqct import *
import numpy as np

class TestSlurm(unittest.TestCase):
    def setUp(self):
        self.wat = Mol(["O", "H", "H"], 0, 1)
        self.cl  = Mol(["Cl"], -1)

    def test_molecule(self):
        self.assertEqual(0, self.wat.charge(), "Incorrect Mol.charge()")
        self.assertEqual(1, self.wat.spin(), "Incorrect Mol.spin()")

    def test_io(self):
        topol = [
            [(6,), self.cl],
            [(0,2,3), self.wat],
            [(1,4,5), self.wat]
        ]
        frames = read_cp2k("examples/chloride/cp2k_cluster.xyz", topol)
        self.assertEqual(701, len(frames))

if __name__ == '__main__':
    unittest.main()
