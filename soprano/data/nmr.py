# Soprano - a library to crack crystals! by Simone Sturniolo
# Copyright (C) 2016 - Science and Technology Facility Council

# Soprano is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Soprano is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
NMR Data

Data on NMR relevant properties of elements and isotopes - spin, gyromagnetic
ratio and quadrupole moment.
"""

import re
import json
import pkgutil
import numpy as np
import scipy.constants as cnst

# EFG conversion constant.
# Units chosen so that EFG_TO_CHI*Quadrupolar moment*Vzz = Hz
EFG_TO_CHI = cnst.physical_constants['atomic unit of electric field '
                                     'gradient'][0]*cnst.e*1e-31/cnst.h

try:
    _nmr_data = pkgutil.get_data('soprano',
                                 'data/nmrdata.json').decode('utf-8')
    _nmr_data = json.loads(_nmr_data)
except IOError:
    _nmr_data = None


def _get_nmr_data() -> dict:
    """
    Helper method to ensure the proper NMR data is loaded from
    soprano/data/nmrdata.json
    :return: Data on NMR relevant properties of elements and isotopes
        - spin, gyromagnetic ratio and quadrupole moment
    """
    if _nmr_data is not None:
        return _nmr_data
    else:
        raise RuntimeError('NMR data not available. Something may be '
                           'wrong with this installation of Soprano')


def _get_isotope_data(elems, key, isotopes=None, isotope_list=None,
                      use_q_isotopes=False) -> list:
    """
    Provides isotope_data about isotopes from soprano/isotope_data/nmrdata.json
    TODO: Add support for multiple isotopes of same element
        Is this covered by using isotope_list instead?

    :param elems: Elements to get isotope_data for
    :type elems: Union(list[str] | str)
    :param key: Value to retrieve - ["gamma", "I", "Q"]
    :type key: str
    :param isotopes: Elements with their corresponding isotopes Z numbers
    :type isoptopes: dict
    :param isotope_list: Ordered list of isotopes corresponding to elems, optional
    :type isotope_list: list
    :param use_q_isotopes: Whether to use the quadrupole
    :type use_q_isotopes: bool
    :return: Isotope isotope_data for the requested elements and parameters
    :rtype: list
    """

    if isinstance(elems, str):
        elems = [elems]  # It's a single element

    isotope_data = np.zeros(len(elems))
    nmr_data = _get_nmr_data()
    for i, element in enumerate(elems):
        if element not in nmr_data:
            # Non-existing element
            raise RuntimeError('No NMR data on element {0}'.format(element))

        element_data = nmr_data[element]
        # Uses most common isotope if not otherwise specified
        iso_Z = element_data['iso']
        # Refactored to use elifs - in soprano/properties/nmr/efg docstrings explain
        # precedence. Having "if" instead of "elif" did not respect this priority
        # Uses the quadrupole isotope as the Z number
        if use_q_isotopes and element_data['Q_iso'] is not None:
            iso_Z = element_data['Q_iso']
        # Uses the isotopes parameter if its set
        elif (isotopes is not None) and (element in isotopes):
            iso_Z = isotopes[element]
        # Uses the isotope list parameter if its set
        elif isotope_list is not None and isotope_list[i] is not None:
            iso_Z = isotope_list[i]

        try:
            # Retrieve requested isotope parameter to return
            isotope_data[i] = element_data[str(iso_Z)][key]
        except KeyError:
            raise RuntimeError('Data {0} does not exist for isotope {1} of '
                               'element {2}'.format(key, iso_Z, element))
    return isotope_data


def _el_iso(sym):
    """ Utility function: split isotope and element in conventional
    representation.
    """

    nmr_data = _get_nmr_data()
    match = re.findall('([0-9]*)([A-Za-z]+)', sym)
    if len(match) != 1:
        raise ValueError(f'Invalid isotope symbol - Use a single element, not a formula')
    elif match[0][1] not in nmr_data:
        raise ValueError('Invalid element symbol')

    el = match[0][1]
    # What about the isotope?
    iso = str(nmr_data[el]['iso']) if match[0][0] == '' else match[0][0]

    if iso not in nmr_data[el]:
        raise ValueError('No data on isotope {0} for element {1}'.format(iso,
                                                                         el))

    return el, iso


def nmr_gamma(elems, iso=None):
    """Gyromagnetic ratio for an element
    N.B This cannot handle a system with the different isotopes of the same element
        dict keys being the element means the isotopes overload each other and only
        one dict entry is made
        - Is this ever going to be an edge case or can we ignore it?

    Return the gyromagnetic ratio for the given element and isotope, in
    rad/(s*T)

    | Args:
    |   elems (str|list(str)):  element symbol(s)
    |   iso (int|list(int)):  isotopes. Default is the most abundant one.

    | Returns:
    |   gamma (float):  gyromagnetic ratio in rad/(s*T)    
    """
    if isinstance(iso, int) or isinstance(iso, str):
        iso = [iso]  # It's a single element
    if isinstance(elems, str):
        elems = [elems]  # It's a single element
    if iso:
        if len(iso) != len(elems):
            raise ValueError("Elements and isotopes dimension mismatch")
    return _get_isotope_data(elems, 'gamma', isotope_list=iso)


def nmr_spin(elems, iso=None):
    """Nuclear spin for an element

    Return the nuclear spin for the given element and isotope, in
    Bohr magnetons

    | Args:
    |   elems (str):  element symbol
    |   iso (int):  isotope. Default is the most abundant one.

    | Returns:
    |   I (float):  nuclear spin in Bohr magnetons
    """
    if isinstance(iso, int) or isinstance(iso, str):
        iso = [iso]  # It's a single element
    if isinstance(elems, str):
        elems = [elems]  # It's a single element
    if iso:
        if len(iso) != len(elems):
            raise ValueError("Elements and isotopes dimension mismatch")
    return _get_isotope_data(elems, 'I', isotope_list=iso)


def nmr_quadrupole(elems, iso=None):
    """Quadrupole moment for an element

    Return the quadrupole moment for the given element and isotope, in
    barns

    | Args:
    |   elems (str):   element symbol
    |   iso (int):  isotope. Default is the most abundant one.

    | Returns:
    |   Q (float):  quadrupole moment in barns
    """
    if isinstance(iso, int) or isinstance(iso, str):
        iso = [iso]  # It's a single element
    if isinstance(elems, str):
        elems = [elems]  # It's a single element
    if iso:
        if len(iso) != len(elems):
            raise ValueError("Elements and isotopes dimension mismatch")
    return _get_isotope_data(elems, 'Q', isotope_list=iso)
