import scipy as sci
import numpy as np
import matplotlib.pyplot as plt

from ase import Atoms
from ase.io import read, write
class DipolarField():
    """A class for sampling Dipolar Field distributions from atomic structures"""

    def __init__(self, info):
        """Init method for class"""
        self.info = info

    @staticmethod
    def _from_atoms(atoms):
        """init method for use with ASE.AtomsCollection
        """
        return atoms
        
        
    def frac_to_abs(self, coords):
        """Converts fractional coordinates to absolute"""
        return coords
    
    def abs_to_frac(self, coords):
        """Converts fractional coordinates to absolute"""        
        return coords
    
    def set_gyro_ratios(self):
        """Sets the gyromagnetic ratios of the spins depending on the Element/isotope list given
        """ 
        for symbol in self.symbols:
            self.gyro_ratios.append(get_gyro_from_symbol)


if __name__ == "__main__":
    si_atoms = read("soprano/calculate/dipolar/Si.json", index=":", format="json")
    print(si_atoms)
    si_obj = DipolarField._from_atoms(si_atoms)
