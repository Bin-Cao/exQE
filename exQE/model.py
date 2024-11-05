from ase.db import connect
import spglib
from ase import Atoms
import numpy as np
from ase.io import read,write
from ase.spacegroup import get_spacegroup, Spacegroup
from pymatgen.core.structure import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os
import numpy as np
import matplotlib.pyplot as plt
import math

class Cell(object):
    def __init__(self, ):
        """
        Initialize the Cell class.
        """
        # self.atom = atom  # Placeholder for atom, if needed

    def cif2db(self, loc, cry_name):
        """
        Parse a CIF file and save the structure in a database.

        Parameters:
        loc (str): The location of the CIF file to be read. e.g., './BaTe.cif'.
        name (str): The name of the crystal. e.g., 'BaTe'.

        Returns:
        ASE Atoms object: The atoms from the database .
        """
        try:
            name = cry_name + '.db'
            # Connect to the database file or create it if it doesn't exist
            db = connect(name)
            # Read the CIF file to get the atomic structure
            atoms = read(loc)
            # Write the atomic structure to the database with additional metadata
            db.write(atoms, format='cif', remark_cb='caobin')
            print('The CIF structure has been parsed successfully and saved in the database as conventional unit cell.')

            # Reconnect to the database to retrieve the stored structure
            database = connect("./" + name)
            # Return the atoms with the specified ID
            atom = database.get_atoms(id=1)
            os.remove("./" + name)
            return atom
        
        except Exception as e:
            print(f"Error processing {loc}: {e}")
        

    def conv2prim(self, conv_atom):
        """
        Convert a conventional cell to a primitive cell.

        Parameters:
            conv_atom (Atoms): The conventional atom defined in the atomic simulation unit (asu).

        Returns:
            tuple: Lattice constants, primitive lattice cell matrix in Cartesian coordinates,Atoms attribute
        """
        lattice = conv_atom.get_cell()
        positions = conv_atom.get_scaled_positions()
        numbers = conv_atom.get_atomic_numbers()
        cell = (lattice, positions, numbers)
        primitive_uc = spglib.standardize_cell(cell, to_primitive=True, no_idealize=True)
        prim_lattice, prim_positions, prim_numbers = primitive_uc
        prim_atoms = Atoms(cell=prim_lattice, scaled_positions=prim_positions, numbers=prim_numbers, pbc=True)
        lc = prim_atoms.cell.cellpar()
        lmtx = prim_atoms.get_cell()[:]
        return lc, lmtx,prim_atoms

    def prim2conv(self, prim_atom):
        """
        Convert a primitive cell to a conventional cell.

        Parameters:
            prim_atom (Atoms): The primitive atom defined in the atomic simulation unit (asu).

        Returns:
            tuple: Lattice constants, conventional lattice cell matrix in Cartesian coordinates, Atoms attribute
        """
        lattice = prim_atom.get_cell()
        positions = prim_atom.get_scaled_positions()
        numbers = prim_atom.get_atomic_numbers()
        cell = (lattice, positions, numbers)
        conventional_cell = spglib.standardize_cell(cell, to_primitive=False, no_idealize=True)
        conv_lattice, conv_positions, conv_numbers = conventional_cell
        conventional_atoms = Atoms(cell=conv_lattice, scaled_positions=conv_positions, numbers=conv_numbers, pbc=True)
        lc = conventional_atoms.cell.cellpar()
        lmtx = conventional_atoms.get_cell()[:]
        return lc, lmtx, conventional_atoms

    def slab_cell(self, surface_index, atom, slab_thickness=30.0, min_vacuum_size=20.0, atom_type='primitive'):
        """
        Create a slab cell.

        Parameters:
            surface_index (tuple): Tuple of lattice indices, e.g., (1, 1, 0).
            atom (Atoms): Primitive or conventional atom defined in atomic simulation environment (ase).
            slab_thickness (float): Minimum size of layers containing atoms, in angstroms. Default is 30.0.
            min_vacuum_size (float): Minimum vacuum layer size, in angstroms. Default is 20.0.
            atom_type (str): Type of the atom, 'primitive' or 'conventional'. Default is 'primitive'.

        Returns:
            Slab: slab cell.
        """
        if atom_type == 'primitive':
            _, _, atoms = self.prim2conv(atom)
            structure = _atom2str(atoms)
        elif atom_type == 'conventional':
            structure = _atom2str(atom)
        else:
            raise ValueError("Error: unknown atom type!")
       
        # Perform space group analysis on the crystal structure
        #analyzer = SpacegroupAnalyzer(structure)
        #optimized_structure = analyzer.get_primitive_standard_structure()
        

        # Generate the slab using the optimized structure
        slab_generator = SlabGenerator(
            initial_structure=structure,
            miller_index=surface_index,
            min_slab_size=slab_thickness,
            min_vacuum_size=min_vacuum_size,
            center_slab=True,
            in_unit_planes=True,
            primitive=False
        )
       

        # Get the slab from the generator
        slab = slab_generator.get_slab()
        # slab.remove_sites([len(slab) - 1])

        print('\n\nThe QUANTUM-ESPRESSO slab template: \n')

        lattice = slab.lattice
        A, B, C = lattice.matrix

        print("CELL_PARAMETERS (angstrom) % /conventional matrix in cartesian coordinates")
        print(f"  {A[0]:.10f}  {A[1]:.10f}  {A[2]:.10f}")
        print(f"  {B[0]:.10f}  {B[1]:.10f}  {B[2]:.10f}")
        print(f"  {C[0]:.10f}  {C[1]:.10f}  {C[2]:.10f}")

        # 找到最接近 (0.5, 0.5, 0.5) 的原子
        target_coords = np.array([0.5, 0.5, 0.5])
        min_distance = float('inf')
        closest_site_index = -1

        for i, site in enumerate(slab.sites):
            frac_coords = np.array(site.frac_coords)
            distance = np.linalg.norm(frac_coords - target_coords)
            if distance < min_distance:
                min_distance = distance
                closest_site_index = i

        print("ATOMIC_POSITIONS (crystal) % /fractional positions")
        for i, site in enumerate(slab.sites):
            species = site.species_string
            frac_coords = site.frac_coords
            if i == closest_site_index:
                print(f"  {species}  {frac_coords[0]:.10f}  {frac_coords[1]:.10f}  {frac_coords[2]:.10f}  0 0 0")
            else:
                print(f"  {species}  {frac_coords[0]:.10f}  {frac_coords[1]:.10f}  {frac_coords[2]:.10f}")

        print("\n\nATOMIC_POSITIONS (angstrom) % /atom positions in cartesian coordinates")
        for i, site in enumerate(slab.sites):
            species = site.species_string
            cart_coords = site.coords
            if i == closest_site_index:
                print(f"  {species}  {cart_coords[0]:.10f}  {cart_coords[1]:.10f}  {cart_coords[2]:.10f}  0 0 0")
            else:
                print(f"  {species}  {cart_coords[0]:.10f}  {cart_coords[1]:.10f}  {cart_coords[2]:.10f}")

        return slab,

    def plot_band_structure(self, file_path, Ef, ylim, vline, k_points=('L', r'${\Gamma}$', 'X', 'U|K', r'${\Gamma}$')):
        """
        Plot the band structure from the given file.

        Parameters:
        file_path (str): Path to the file containing band structure data.
        Ef (float): Fermi energy.
        ylim (list): Minimum and Maximum value for the y-axis.
        vline (list): Positions of vertical lines, according to the grids on k space.
        k_points (tuple): Labels for high-symmetry points in k space.

        Example:
        import exQE.model as qe
        file_path = "cuband.dat"
        Ef = 12.5663
        ylim = [-10, 30]
        vline = [0, 50, 100, 120, 170]
        k_points = ('L', r'${\Gamma}$', 'X', 'U|K', r'${\Gamma}$')

        qe.Cell().plot_band_structure(file_path, Ef, ylim, vline, k_points)

        Returns:
        None
        """
        # Open the file
        with open(file_path) as feig:
            # Read the first line and extract number of bands and k-points
            first_line = feig.readline()
            nbnd = int(first_line.split(',')[0].split('=')[1])
            nks = int(first_line.split(',')[1].split('=')[1].split('/')[0])

            # Initialize a 2D array to store band data
            eig = np.zeros((nks, nbnd), dtype=float)

            # Read band data
            for i in range(nks):
                feig.readline()  # Skip the k-point line
                count = 0
                for j in range((nbnd - 1) // 10 + 1):
                    data_line = feig.readline()
                    for value in data_line.split():
                        eig[i][count] = float(value)
                        count += 1

        # Create a figure and set its size
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)

        # Set the x and y axis limits
        ax.set_xlim([0, nks - 2])  # Avoid additional lines at the end
        ax.set_ylim(ylim)

        # Set the y-axis label
        ax.set_ylabel(r'E (eV)', fontsize=16, fontweight='bold')

        # Plot the band structure
        for i in range(nbnd):
            ax.plot(eig[:, i] - Ef, color='b', linewidth=1.5)

        # Draw vertical lines at specified positions
        for x in vline:
            ax.axvline(x, ymin=ylim[0], ymax=ylim[1], linewidth=1.2, color='black')

        # Set x-axis ticks and labels
        ax.set_xticks(vline)
        ax.set_xticklabels(k_points, fontsize=14, fontweight='bold')

        # Set other properties for scientific journal style
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        # Save the figure
        plt.savefig('band_structure.png', dpi=300, bbox_inches='tight')

        # Display the plot
        plt.show()

    def plot_density_of_states(self, file_path, Ef, xlim, ylim):
        """
        Plot the density of states from the given file.

        Parameters:
        file_path (str): Path to the file containing density of states data.
        Ef (float): Fermi energy.
        xlim (list): Minimum and Maximum value for the x-axis.
        ylim (list): Minimum and Maximum value for the y-axis.

        Example:
        import exQE.model as qe
        file_path = "cu.dat"
        Ef = 12.5663
        xlim = [0,6]
        ylim = [-10, 30]

        qe.Cell().plot_density_of_states(file_path, Ef, xlim, ylim)

        Returns:
        None
        """
        # Create a figure and set its size
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)

        # Set the x and y axis limits
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # Set the x-axis label
        ax.set_xlabel(r'DOS (States/atom/eV)', fontsize=16, fontweight='bold')

        # Set the y-axis label
        ax.set_ylabel(r'Energy (eV)', fontsize=16, fontweight='bold')

        # Load the density of states data, ignoring lines starting with '#'
        dos_data = np.loadtxt(file_path, comments='#')

        # Plot the density of states
        ax.plot(dos_data[:, 1], dos_data[:, 0] - Ef, color='b', linewidth=1.5)

        # Set other properties for scientific journal style
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        # Save the figure
        plt.savefig('DOS.png', dpi=300, bbox_inches='tight')

        # Display the plot
        plt.show()

################################################################
########################## Functions ###########################
################################################################


def _atom2str(atoms):
    """
    Convert ASE Atoms object to Pymatgen Structure object.

    Args:
        atoms (Atoms): ASE Atoms object containing information about the structure.

    Returns:
        Structure: Pymatgen Structure object created from the ASE Atoms object.
    """
    # Extract lattice information
    cell = atoms.get_cell()
    
    # Extract element symbols and positions
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_scaled_positions()
    
    # Create Pymatgen Lattice object using the extracted cell information
    lattice = Lattice(cell)
    
    # Create Pymatgen Structure object using the extracted information
    structure = Structure(lattice, symbols, positions)
    
    return structure

########################## uncalled functions ###########################
def miller_to_cartesian(lattice_cell, miller_indices):
    """
    Convert Miller indices to Cartesian coordinates.

    Parameters:
        lattice_cell (array): 3x3 array of lattice vectors.
        miller_indices (list): List of Miller indices, e.g., [h, k, l].

    Returns:
        array: Cartesian coordinates of the lattice vector.
    """
    h, k, l = miller_indices
    cart_vector = np.dot(miller_indices, lattice_cell)
    return cart_vector / np.linalg.norm(cart_vector)


def rotate_hkl_to_z(lattice_cell, miller_indices):
    """
    Compute a rotation matrix to rotate a crystal direction (hkl) to the z-axis.

    Parameters:
        lattice_cell (array): 3x3 array of lattice vectors.
        miller_indices (list): List of Miller indices, e.g., [h, k, l].

    Returns:
        array: 3x3 rotation matrix.
    """
    
    cart_vector = miller_to_cartesian(lattice_cell, miller_indices)

    # z-axis unit vector
    z_axis = np.array(lattice_cell[2], dtype=float)

    # Calculate the rotation axis (cross product of hkl vector and z-axis)
    axis = np.cross(cart_vector, z_axis)
    

    # If hkl is already in the z-axis direction, return the identity matrix
    if np.linalg.norm(axis) <= 1e-5:
        return np.eye(3)
    
    axis = axis / np.linalg.norm(axis)
    # Calculate the rotation angle (angle between hkl vector and z-axis)
    angle = np.arccos(np.dot(cart_vector, z_axis))

    # Calculate the rotation matrix using Rodrigues' rotation formula
    Rk = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])

    rotation_matrix = np.eye(3) + np.sin(angle) * Rk + (1 - np.cos(angle)) * np.dot(Rk, Rk)
    return rotation_matrix

def grep_latt_atoms(atomic_numbers, kinds):
    """
    Return a list of atomic indices in the periodic table for given kinds.

    Parameters:
        atomic_numbers (list): List of atomic lables of all atoms in the unit cell.
        kinds (list): List of kinds (indices) to retrieve atomic numbers for.

    Returns:
        list: List of atomic lables corresponding to the specified kinds.
    """
    atoms_list = []
    for kind in kinds:
        atoms_list.append(atomic_numbers[kind])
    return atoms_list



def calc_cell_volume(a, b, c, alpha, beta, gamma):
    # convert to rad
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)
    
    volume = (a * b * c) * math.sqrt(1 - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma))
    
    return volume
