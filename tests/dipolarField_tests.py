#!/usr/bin/env python
"""
Test code for DipolarField
"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import unittest
import numpy as np
from ase import io
from soprano.selection import AtomSelection
from soprano.calculate.dipolar import DipolarField
import ase

class TestDipolarField(unittest.TestCase):
    si_atoms = ase.io.read("test_data/Si.json", index=":", format="json")[0]
    muon_pos = [0.1, 0.1, 0.1]
    field = DipolarField(atoms=si_atoms, mu_pos=muon_pos, cutoff=11)

    def setUp(self) -> None:
        self.field = TestDipolarField.field.clone()

    def tearDown(self) -> None:
        del self.field

    def test_init(self):
        """
        TODO: Add tests for all attr that are calculated and not set directly from args
        :return:
        """
        # Check properties are set correctly
        self.assertTrue(np.equal(np.array(TestDipolarField.muon_pos), self.field.mu_pos).all)
        self.assertTrue(hasattr(self.field, "cell"))
        self.assertTrue(hasattr(self.field, "scell"))
        self.assertTrue(hasattr(self.field, "nmr_gammas"))
        self.assertTrue(hasattr(self.field, "nmr_spin"))

    def test_add_supercell_atoms(self):
        """
        """
        # Add one unit cell in each direction (Make 3x3x3 lattice of unit cells)
        self.field.add_supercell_atoms(n=1)
        # Check that duplicate atoms are ignored - cant have two atoms at same site
        self.assertEqual(len(self.field.atom_pos), 27)
        # Increase size of supercell
        N = 3
        self.field.add_supercell_atoms(n=N)
        self.assertEqual(len(self.field.atom_pos), len(range(-N, N+1))**3)

    def test_get_dipolar_tensors(self):
        num_cells = 1
        num_atoms = (2*num_cells + 1)**3
        self.field.add_supercell_atoms(n=num_cells)
        pairs = [(i, j)
                 for i in range(num_atoms)
                 for j in range(num_atoms)
                 if i < j]
        dip_ten = self.field.get_dipolar_tensors()
        pairs = set(pairs)
        keys = set(dip_ten.keys())
        # Ensure that the number of dipolar couplings matches the number of possible atom interactions
        # Gets very expensive for large num_cells
        self.assertTrue(keys == pairs)

    def test_plot(self):
        pass
        # import matplotlib.pyplot as plt
        # self.field.add_supercell_atoms(n=1)
        # dip_ten = self.field.get_dipolar_tensors()
        #
        # from ase.visualize import view
        # hyd1 = ase.Atoms(["H"], scaled_positions=[[0, 0, 0]], cell=[1, 1, 1])
        # hyd1 *= 2
        # del hyd1[[1, 3, 4, 5, 6, 7]]
        # hyd1.center()
        # print(hyd1)
        # view(hyd1)
        # # full_tensor = [np.sum(x) ]
        # # for i in dip_ten.values():
        # #     print(i)
        # #     break
        # # print(np.sum(dip_ten.values(), axis=0))
        # # print(dip_ten[(0, 1)])
        # # print(dip_ten[(0, 2)])




if __name__ == '__main__':
    unittest.main()
