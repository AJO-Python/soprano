import numpy as np
import scipy.constants as cnst

from ase import Atom
from ase import Atoms
from ase.io import read, write
from ase.visualize import view

class DipolarField():
    """
    Provides information about the magnetic state of an ASE Atoms collection
    Focus on dipolar
    """
    def __init__(self, atom):
        """
        :param Atom atom: ASE atom object
        """
        self.atom = atom

    
    def get_field_at_point(self, point):
        """Gets resultant field from the Atom magmom at a given location
        B(r) = (mu_0/4pi) * [ ((3r(m dot r))/r^5) - (m/r^3) ]
                ^                   ^                   ^
                mag_perm            term1               term2
        :param array point: 
        """
        moment = self.atom.magmom
        mag_perm = cnst.mu_0 / (4 * cnst.pi)
        relative_loc = np.subtract(point, self.atom.position)
        rel_magnitude = np.linalg.norm(relative_loc)
        m_dot_r = np.dot(moment, relative_loc)

        term1 = np.array(
            (3 * relative_loc * m_dot_r)
            / (rel_magnitude ** 5)
                        )
        term2 = np.array(moment / (rel_magnitude ** 3))
        
        B_field = np.subtract(term1, term2)
        B_field *= mag_perm
        return B_field


class DipolarFields():
    """
    Provides information about the magnetic state of an ASE Atoms collection
    Takes each atom as a point dipole with moment == atom.magmom
    Does NOT account for boundary conditions or unit/super cells
    """

    def __init__(self, atoms):
        """
        :param Atom atom: ASE atom object
        """
        self.atoms = atoms
    
    def atom_position_1d(self, dim="x"):
        """
        Gets list of atom positions in single dimension

        :param dim: coordinate plane to get
        :type dim: string, default="x", options=["x", "y", "z"]
        :return: positions of all atoms in 1d
        :rtype: list
        """
        dims={"x": 0, "y": 1, "z": 2}
        try:
            return [pos[dims[dim.lower()]] for pos in self.atoms.positions]
        except KeyError:
            raise('Please specify a dimension ["x", "y", "z"]')


    def get_field_at_point(self, point):
        """Gets resultant field from the Atom magmom at a given location
        B(r) = (mu_0/4pi) * [ ((3r(m dot r))/r^5) - (m/r^3) ]
                ^                   ^                   ^
                mag_perm            term1               term2
        :param array point: 
        """
        B_field = np.empty(shape=3)
        for position, moment in zip(self.atoms.positions, self.atoms.get_initial_magnetic_moments()):
            mag_perm = cnst.mu_0 / (4 * cnst.pi)
            relative_loc = point - position
            rel_magnitude = np.linalg.norm(relative_loc)
            m_dot_r = np.dot(moment, relative_loc)

            term1 = np.array(
                (3 * relative_loc * m_dot_r)
                / (rel_magnitude ** 5)
                            )
            term2 = np.array(moment / (rel_magnitude ** 3))
            
            B_field += term1 - term2
        B_field *= mag_perm
        return B_field

if __name__=="__main__":
    import matplotlib.pyplot as plt
    import matplotlib

    d = 0.9
    atom = Atom(symbol="H", position=(0, 0, 0), magmom=[0, 1, 0])
    #atoms = Atoms(symbols="H", positions=[(0, 0, 0), (0.3, 0, 0)], cell=[d, 0, 0], pbc=[1, 0, 0])

    atoms_pos = [()]
    atoms = Atoms(symbols=("H", "H", "H", "H"),
            positions=[(-0.5, 0, 0), (0.5, 0, 0), (-0.5, 1, 0), (0.5, 1, 0)],
            magmoms=[(0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0)]
            )
    hydrogen = DipolarField(atom)
    multi_atoms = DipolarFields(atoms)
    from ase.visualize import view
    view(multi_atoms.atoms)
    nx, ny = 64, 64
    x = np.linspace(-5, 5, nx)
    y = np.linspace(-5, 5, ny)
    fields = {"x": np.empty(shape=(nx, ny)), "y": np.empty(shape=(nx, ny)), "z": np.empty(shape=(nx, ny))}

    # Populate field arrays
    for i, x_val in enumerate(x):
        for j, y_val  in enumerate(y):
            fields["x"][j][i], fields["y"][j][i], fields["z"][j][i] = multi_atoms.get_field_at_point([x_val, y_val, 0.1])
    
    # Plot streamplot
    color = 2 * np.log(np.hypot(fields["x"], fields["y"]))
    plt.streamplot(x, y, fields["x"], fields["y"], color=color, linewidth=1, cmap=plt.cm.inferno, density=2, arrowstyle='->', arrowsize=1.5)
    xs=multi_atoms.atom_position_1d("x")
    ys=multi_atoms.atom_position_1d("y")
    plt.scatter(x=xs, y=ys)
    plt.show()
    
    # Plot field distribution
    fig, axes = plt.subplots(3, 1)
    for ax, field in zip(axes, fields.items()):
        flat = field[1].flatten()
        ax.hist(flat, bins=50, range=np.percentile(flat, [5, 95]))
        ax.set_xlabel(f"Field component {field[0]}")
        ax.set_ylabel(f"Frequency")
        #ax.set_xlim(-5e-10, 5e-10)
    plt.tight_layout()
    plt.show()
