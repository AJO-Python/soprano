#!/usr/bin/env python
"""
Test code for the AtomsProperty class and its children
"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import glob
from ase import Atoms
from ase.io import read, write
sys.path.insert(0, os.path.abspath(
                   os.path.join(os.path.dirname(__file__), "../")))  # noqa
from soprano.collection import AtomsCollection
from soprano.properties import AtomsProperty
import unittest
import numpy as np

_TESTDATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "test_data")


class TestPropertyLoad(unittest.TestCase):

    def test_dummyprop(self):

        # Define a dummy derived property and test it

        class DummyProperty(AtomsProperty):

            default_name = 'dummy'
            default_params = {'mul': 2.0}

            @staticmethod
            def extract(s, mul):
                return s.positions.shape[0]*mul

        # Now two atoms objects to test it on
        a1 = Atoms('C')
        a2 = Atoms('CC')

        dummyDoubled = DummyProperty(name='doubledummy', mul=4)

        self.assertEqual(DummyProperty.get(a1), 2.0)
        self.assertEqual(dummyDoubled(a2), 8.0)

        # Also check that improper parameters are rejected
        self.assertRaises(ValueError, DummyProperty, wrong='this is wrong')

        # Test behaviour on a collection instead
        c1 = AtomsCollection([a1, a2])
        c2 = AtomsCollection([a2, a1])

        DummyProperty.get(c1, store_array=True)
        dummyDoubled(c2, store_array=True)
        self.assertTrue(np.all(c1.get_array('dummy') == [2, 4]))
        self.assertTrue(np.all(c2.get_array('doubledummy') == [8, 4]))    

    def test_basicprop(self):

        from soprano.utils import cart2abc
        from soprano.properties.basic import (LatticeCart, LatticeABC,
                                              CalcEnergy)

        cell = np.array([[1, 0, 0], [0, 2, 0], [0, 0.5, 3.0]])
        a = Atoms(cell=cell)
        linLatt = LatticeCart(shape=(9,))
        degLatt = LatticeABC(deg=True)

        ans = cell.reshape((9,))
        self.assertTrue(np.all(linLatt(a) == ans))
        ans = cart2abc(cell)
        ans[1, :] *= 180.0/np.pi
        self.assertTrue(np.all(degLatt(a) == ans))

    def test_propertymap(self):

        from soprano.properties.basic import NumAtoms

        num_atoms = np.random.random_integers(1, 10, 10)
        coll = AtomsCollection([Atoms('H'*n) for n in num_atoms])

        num_atoms_prop = coll.all.map(NumAtoms.get)

        self.assertTrue((num_atoms == num_atoms_prop).all)

    def test_linkageprops(self):

        from soprano.properties.linkage import (LinkageList,
                                                Molecules, MoleculeNumber,
                                                MoleculeMass,
                                                MoleculeCOMLinkage,
                                                MoleculeRelativeRotation)
        from soprano.properties.transform import Rotate

        a = read(os.path.join(_TESTDATA_DIR, 'mol_crystal.cif'))

        mols = Molecules.get(a)

        # for i, m in enumerate(mols):
        #     ai = m.subset(a)
        #     ai.arrays['positions'] += np.tensordot(ai.get_cell(),
        #                                            ai.arrays['cell_indices'],
        #                                            axes=(0,1)).T
        #     #write('mol{0}.cell'.format(i+1), ai)
        #     #print(ai.get_chemical_symbols())
        #     print(m.get_array('cell_indices'))

        self.assertTrue(MoleculeNumber.get(a) == 4)
        self.assertTrue(np.isclose(MoleculeMass.get(a), 142.06788).all())
        self.assertTrue(len(MoleculeCOMLinkage.get(a)) == 6)

        #print(MoleculeRelativeRotation.get(a))

    def test_transformprops(self):

        from ase.quaternions import Quaternion
        from soprano.selection import AtomSelection        
        from soprano.properties.transform import (Translate, Rotate, Mirror)

        a = Atoms('CH', positions=[[0,0,0],[0.5,0,0]])

        sel = AtomSelection.from_element(a, 'C')
        transl = Translate(selection=sel, vector=[0.5,0,0])        
        rot = Rotate(selection=sel, center=[0.25,0.0,0.25],
                     quaternion=Quaternion([np.cos(np.pi/4.0),
                                            0,
                                            np.sin(np.pi/4.0),
                                            0]))
        mirr = Mirror(selection=sel, plane=[1, 0, 0, -0.25])

        aT = transl(a)
        aR = rot(a)
        aM = mirr(a)

        self.assertAlmostEqual(np.linalg.norm(aT.get_positions()[0]),
                               np.linalg.norm(aT.get_positions()[1]))
        self.assertAlmostEqual(np.linalg.norm(aR.get_positions()[0]),
                               np.linalg.norm(aR.get_positions()[1]))
        self.assertAlmostEqual(np.linalg.norm(aM.get_positions()[0]),
                               np.linalg.norm(aM.get_positions()[1]))

if __name__ == '__main__':
    unittest.main()
