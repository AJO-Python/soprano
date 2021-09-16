import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cnst

from ase import Atom
from ase import Atoms
from ase.io import read, write
from ase.visualize import view

from soprano.utils import minimum_supcell, supcell_gridgen


class DipolarLattice:
    def __init__(
        self,
        positions: np.ndarray,
        cell: Atoms.cell,
        gammas: list,
        cutoff: float = 10,
        ordered: bool = False,
    ):
        """
        :param positions: 3d locations of dipoles
        :type positions: np.ndarray
        :param cell: Describes the unit cell of the lattice
        :type cell: ase.atoms.cell
        :param gammas: Gyromagnetic ratio. Either a single value (homogenous dipoles) or of the same length as positions
        :type gammas: list(str)
        :param ordered: False if random dipole orientations
        :type ordered: bool
        """
        self.positions = np.array(positions)
        self.cell = cell
        self.gammas = gammas
        self.cutoff = cutoff
        self.scell = minimum_supcell(cutoff, self.cell)
        self.scell_grid_frac, self.scell_grid_cart = supcell_gridgen(
            self.cell, self.scell
        )
        print(self.scell_grid_cart.shape)
        # print(self.cutoff)

    @staticmethod
    def from_atoms(atoms: Atoms):
        """Generates a DipolarLattice from an ASE.Atoms collection

        :param atoms: Atoms struture to generate spins and lattice from
        :type atoms: Atoms
        :return: DipolarLattice
        :rtype: DipolarLattice
        """
        # TODO: write extract_info_from_atoms()
        positions, cell, gammas = extract_info_from_atoms(atoms)
        return DipolarLattice(positions, cell, gammas)

    def extract_info_from_atoms(self, atoms: Atoms):
        positions = atoms.get_positions()
        for atom in atoms.get_initial_magnetic_moments:
            pass

    def get_field_distribution(self):
        """get field distribtion"""
        pass

    def get_coords_inside_sphere(self):
        print("positions: ", self.positions)
        r = self.positions[:, None, :]
        # print(self.positions)
        # print(r)


if __name__ == "__main__":
    from ase.build import bulk
    from ase.visualize import view

    cell = bulk("Au").cell
    # print(cell)
    dipole_pos = (
        [0, 0, 0],
        [0, 2.04, 0],
        [2.04, 2.04, 0],
        [2.04, 0, 0],
    )
    # lattice_cell = [1, 1, 1]

    test_atoms = Atoms(
        "HHHH",
        positions=dipole_pos,
        cell=cell,
    )
    cell = test_atoms.get_cell()
    # view(test_atoms)
    lattice = DipolarLattice(positions=dipole_pos, cell=cell, gammas=[1])

    lattice.get_coords_inside_sphere()
    print(lattice)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    for i in lattice.scell_grid_cart[::1]:
        ax.scatter(i[0], i[1], i[2], marker=".", alpha=0.2)
    for molecule in test_atoms.positions:
        ax.scatter(*molecule, marker="o", c="r")
    plt.show()
