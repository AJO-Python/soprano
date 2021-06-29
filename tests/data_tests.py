#!/usr/bin/env python
"""
Test code for data retrieval
"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import unittest
import numpy as np
from ase.data import vdw_radii as vdw_radii_ase
from soprano.data import (vdw_radius,
                          nmr_gamma,
                          nmr_spin,
                          nmr_quadrupole,
                          )
from soprano.data.nmr import (_get_nmr_data,
                              _get_isotope_data,
                              _el_iso,
                             )


class TestData(unittest.TestCase):

    def test_get_nmr_data(self):
        output = _get_nmr_data()
        self.assertEqual(type(output), dict)
        # Sanity check that the json loaded has correct Muon data
        self.assertIn("Mu", output)
        self.assertEqual(output["Mu"]["1"]["gamma"], -851615463.149621)
        from soprano.data.nmr import _nmr_data
        # Check that the import fails if the file cannot be loaded
        # _nmr_data = None
        # self.assertRaises(RuntimeError, lambda: _get_nmr_data)

    def test_vdw(self):

        self.assertEqual(vdw_radius('C', 'csd'), 1.77)
        self.assertEqual(vdw_radius('C', 'jmol'), 1.95)
        self.assertEqual(vdw_radius('C', 'ase'), vdw_radii_ase[6])

    def test_gamma(self):
        """nmr_gamma, nmr_spin, nmr_quadrupole can all be refactored into a single function with arg "key"=("gamma"|"I"|"Q")
        """
        # Single values
        self.assertEqual(nmr_gamma('H'), 267522128.0)
        self.assertEqual(nmr_gamma('H', 2), 41066279.1)

        # Array inputs
        gamma_list = nmr_gamma(['H', 'I'], [1, 127])
        self.assertEqual(gamma_list[0], 267522128.0)
        self.assertEqual(gamma_list[1], 53895730.0)

        # Error handling for mismatched element and isotope dims
        self.assertRaises(ValueError, lambda: nmr_gamma(['H', 'S'], [2]))

    def test_spin(self):
        self.assertEqual(nmr_spin('H'), 0.5)
        self.assertEqual(nmr_spin('H', 2), 1)

    def test_quadrupole(self):
        self.assertEqual(nmr_quadrupole('H'), 0)
        self.assertEqual(nmr_quadrupole('H', 2), 2.86)

    def test_get_isotope_data(self):
        """TODO: Refactor to remove the isotope: dict arg. Isotope_list does the same
        thing but allows for multiple isotopes of the same element
        I think it is redundant but its used in a lot of places so will leave for now
        """
        H_iso1_gamma = 267522128.0
        H_iso2_gamma = 41066279.1
        # Basic functionality
        data = _get_isotope_data(elems='H',
                                 key='gamma',
                                 isotopes=None,
                                 isotope_list=None,
                                 use_q_isotopes=True)
        self.assertIsNotNone(data)
        # Element does not exist
        self.assertRaises(RuntimeError, _get_isotope_data,
                          elems='NotAnElement',
                          key='gamma',
                          isotopes=None,
                          isotope_list=None,
                          use_q_isotopes=False)
        # One element of list does not exist
        self.assertRaises(RuntimeError, _get_isotope_data,
                          elems=['NotAnElement', 'H'],
                          key='gamma',
                          isotopes=None,
                          isotope_list=None,
                          use_q_isotopes=False)
        # Incorrect key requested
        self.assertRaises(RuntimeError, _get_isotope_data,
                          elems=['H', 'H'],
                          key='fake_key',
                          isotopes=None,
                          isotope_list=None,
                          use_q_isotopes=True)
        # Isotopes as a dictionary
        data = _get_isotope_data(elems='H',
                                 key='gamma',
                                 isotopes={'H': 2},
                                 isotope_list=None,
                                 use_q_isotopes=False)
        self.assertEqual(data[0], H_iso2_gamma)
        # Isotope as an ordered list
        data = _get_isotope_data(elems=['H', 'H'],
                                 key='gamma',
                                 isotopes=None,
                                 isotope_list=[2, 2],
                                 use_q_isotopes=False)
        self.assertEqual(data[0], H_iso2_gamma)
        # Multiple isotopes of the same element (isotope_list)
        data = _get_isotope_data(elems=['H', 'H'],
                                 key='gamma',
                                 isotopes=None,
                                 isotope_list=[1, 2],
                                 use_q_isotopes=False)
        self.assertEqual(data[0], H_iso1_gamma)
        self.assertEqual(data[1], H_iso2_gamma)
        # Use Q isotopes
        data = _get_isotope_data(elems=['H', 'H'],
                                 key='gamma',
                                 isotopes=None,
                                 isotope_list=None,
                                 use_q_isotopes=True)
        self.assertEqual(data[0], H_iso2_gamma)
        # Prioritise Q_iso over isotope_list (Q_iso=2 for H)
        data = _get_isotope_data(elems=['H', 'H'],
                                 key='gamma',
                                 isotopes=None,
                                 isotope_list=[1, 1],
                                 use_q_isotopes=True)
        self.assertEqual(data[0], H_iso2_gamma)

    def test_el_iso(self):
        # Default isotope used
        element, isotope = _el_iso('H')
        self.assertEqual(element, 'H')
        self.assertEqual(isotope, '1')
        # Request specific isotope
        element, isotope = _el_iso('3H')
        self.assertEqual(element, 'H')
        self.assertEqual(isotope, '3')
        # Giving more than one element
        self.assertRaises(ValueError, _el_iso, '1H1H')
        # Element not in database
        self.assertRaises(ValueError, _el_iso, 'NotAnElement')
        # Isotope doesnt exist
        self.assertRaises(ValueError, _el_iso, '5H')

if __name__ == '__main__':
    unittest.main()
