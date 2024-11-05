from .model import Cell
import os

class Generator(object):
    def __init__(self, pseudo_dir, qe_dir):
        """
        Initialize the Generator object with paths to pseudo potential and quantum-espresso directories.

        :param pseudo_dir: str, the directory of pseudo potential files folder
        :param qe_dir: str, the directory of quantum-espresso folder

        example:

        import exQE.inputFile as file
        pseudo_dir = '../../SSSP'
        qe_dir = '../../qe-7.2'
        generator = file.Generator(pseudo_dir, qe_dir)

        generator.relax('Cu.cif',type='vc-relax')
        generator.relax('Cu.cif', type='relax')
        generator.scf('Cu.cif',type='scf')
        generator.scf('Cu.cif', type='nscf')
        generator.bands('Cu.cif',)
        generator.dos()
        """
        self.pseudo_dir = pseudo_dir
        self.qe_dir = qe_dir
        print('Welcome to apply exQE for DFT calculations!')
        print('quantum-espresso: https://www.quantum-espresso.org/Doc/INPUT_PW.html')

    def relax(self, cif, type='vc-relax', outdir='./exQE'):
        """
        Perform a vc-relax calculation using the given CIF file

        :param cif: str, the location of the CIF file
        :param type: str, the calculation type in QE, default is 'vc-relax'
        :param outdir: str, the directory of output files, default is './exQE'
        """
        if type == 'vc-relax':
            print('Performing vc-relax calculation...')
        elif type == 'relax':
            print('Performing relax calculation...')
        else:
            print('Invalid calculation type, please choose from vc-relax or relax.')
            return
        # Convert CIF file to atomic database and primitive cell
        cov_atom = Cell().cif2db(cif, 'None')
        _, _, pr_atom = Cell().conv2prim(cov_atom)
        print('The primitive unit cell is created')

        # Extract necessary information from the primitive cell
        atoms = pr_atom.get_chemical_symbols()
        fra_coords = pr_atom.get_scaled_positions()
        cell_matrix = pr_atom.get_cell()

        # Prepare the input file content for Quantum Espresso
        input_file_content = f"""
&CONTROL
    calculation='{type}', 
    disk_io='low', 
    prefix='caobin',
    pseudo_dir='{self.pseudo_dir}', 
    outdir='{outdir}', 
    verbosity='high',
    tprnfor=.true., 
    tstress=.true., 
    forc_conv_thr=1.0d-5
/
&SYSTEM
    ibrav=0,
    nat={len(atoms)}, 
    ntyp={len(set(atoms))},
    occupations='smearing', 
    smearing='gauss', 
    degauss=1.0d-2,
    ecutwfc=50, 
    ecutrho=500,
/
&ELECTRONS
    electron_maxstep=100,
    conv_thr=1.0d-9,
    mixing_mode='plain',
    mixing_beta=0.8d0,
    diagonalization='david',
/
&IONS
    ion_dynamics='bfgs'
/
&CELL
    press_conv_thr=0.1
/
ATOMIC_SPECIES
"""

        # Add atomic species information
        for atom in set(atoms):
            atomic_mass = finder_ele_mass(atom)
            pseudo_file = finder_pseudo(self.pseudo_dir, atom)
            input_file_content += f"{atom} {atomic_mass} {pseudo_file}\n"

        # Add cell parameters
        input_file_content += "CELL_PARAMETERS (angstrom)\n"
        for row in cell_matrix:
            input_file_content += ' '.join([f"{x:.9f}" for x in row]) + '\n'

        # Add atomic positions
        input_file_content += "ATOMIC_POSITIONS (crystal)\n"
        for atom, coord in zip(atoms, fra_coords):
            input_file_content += f"{atom} " + ' '.join([f"{x:.9f}" for x in coord]) + '\n'

        # Add k-points
        input_file_content += "K_POINTS {automatic}\n"
        input_file_content += "13 13 13 0 0 0\n"

        if type == 'vc-relax':
            # Save the input file
            with open(f"./vc_relax.in", 'w') as file:
                file.write(input_file_content)

            print(f"Input file 'vc_relax.in' is saved in the current folder")
            print(f"Use command to run the file : mpirun -n 16 {self.qe_dir}/bin/pw.x < vc_relax.in > vc_relax.out")
            print("Check the total energy in relax.out by : grep '!' vc_relax.out")
        else:
            # Save the input file
            with open(f"./relax.in", 'w') as file:
                file.write(input_file_content)

            print(f"Input file 'relax.in' is saved in the current folder")
            print(f"Use command to run the file : mpirun -n 16 {self.qe_dir}/bin/pw.x < relax.in > relax.out")
            print("Check the total energy in relax.out by : grep '!' relax.out")



    def scf(self, cif, type='scf', outdir='./exQE'):
        """
        Perform a self-consistent calculation using the given CIF file

        :param cif: str, the location of the CIF file
        :param type: str, the calculation type in QE, default is 'scf'
        :param outdir: str, the directory of output files, default is './exQE'
        """
        if type == 'scf':
            print('Performing scf calculation...')
        elif type == 'nscf':
            print('Performing no-scf calculation...')
        else:
            print('Invalid calculation type, please choose from scf or nscf.')
            return
        # Convert CIF file to atomic database and primitive cell
        cov_atom = Cell().cif2db(cif, 'None')
        _, _, pr_atom = Cell().conv2prim(cov_atom)
        print('The primitive unit cell is created')

        # Extract necessary information from the primitive cell
        atoms = pr_atom.get_chemical_symbols()
        fra_coords = pr_atom.get_scaled_positions()
        cell_matrix = pr_atom.get_cell()

        # Prepare the input file content for Quantum Espresso
        input_file_content = f"""
&CONTROL
    calculation='{type}', 
    disk_io='low', 
    prefix='caobin',
    pseudo_dir='{self.pseudo_dir}', 
    outdir='{outdir}', 
    verbosity='high',
    tprnfor=.true., 
    tstress=.true., 
    forc_conv_thr=1.0d-5
/
&SYSTEM
    ibrav=0,
    nat={len(atoms)}, 
    ntyp={len(set(atoms))},
    occupations='smearing', 
    smearing='gauss', 
    degauss=1.0d-9,
    ecutwfc=50, 
    ecutrho=500,
/
&ELECTRONS
    electron_maxstep=100,
    conv_thr=1.0d-9,
    mixing_mode='plain',
    mixing_beta=0.8d0,
    diagonalization='david',
/
&IONS
    ion_dynamics='bfgs'
/
&CELL
    press_conv_thr=0.1
/
ATOMIC_SPECIES
"""

        # Add atomic species information
        for atom in set(atoms):
            atomic_mass = finder_ele_mass(atom)
            pseudo_file = finder_pseudo(self.pseudo_dir, atom)
            input_file_content += f"{atom} {atomic_mass} {pseudo_file}\n"

        # Add cell parameters
        input_file_content += "CELL_PARAMETERS (angstrom)\n"
        for row in cell_matrix:
            input_file_content += ' '.join([f"{x:.9f}" for x in row]) + '\n'

        # Add atomic positions
        input_file_content += "ATOMIC_POSITIONS (crystal)\n"
        for atom, coord in zip(atoms, fra_coords):
            input_file_content += f"{atom} " + ' '.join([f"{x:.9f}" for x in coord]) + '\n'

        # Add k-points
        if type == 'scf':
            input_file_content += "K_POINTS {automatic}\n"
            input_file_content += "13 13 13 0 0 0\n"
            # Save the input file
            with open(f"./scf.in", 'w') as file:
                file.write(input_file_content)

            print(f"Input file 'scf.in' is saved in the current folder")
            print(f"Use command to run the file : mpirun -n 16 {self.qe_dir}/bin/pw.x < scf.in > scf.out")
        else:
            input_file_content += "K_POINTS {automatic}\n"
            input_file_content += "53 53 53 0 0 0\n"
            # Save the input file
            with open(f"./nscf.in", 'w') as file:
                file.write(input_file_content)

            print(f"Input file 'nscf.in' is saved in the current folder")
            print(f"Use command to run the file : mpirun -n 16 {self.qe_dir}/bin/pw.x < nscf.in > nscf.out")
       


    def bands(self, cif, outdir='./exQE'):
        """
        Perform a self-consistent calculation using the given CIF file

        :param cif: str, the location of the CIF file
        :param outdir: str, the directory of output files, default is './exQE'
        """
       
        print('Performing energyband calculation...')

        # Convert CIF file to atomic database and primitive cell
        cov_atom = Cell().cif2db(cif, 'None')
        _, _, pr_atom = Cell().conv2prim(cov_atom)
        print('The primitive unit cell is created')

        # Extract necessary information from the primitive cell
        atoms = pr_atom.get_chemical_symbols()
        fra_coords = pr_atom.get_scaled_positions()
        cell_matrix = pr_atom.get_cell()

        # Prepare the input file content for Quantum Espresso
        input_file_content = f"""
&CONTROL
    calculation='bands', 
    disk_io='low', 
    prefix='caobin',
    pseudo_dir='{self.pseudo_dir}', 
    outdir='{outdir}', 
    verbosity='high',
    tprnfor=.true., 
    tstress=.true., 
    forc_conv_thr=1.0d-5
/
&SYSTEM
    ibrav=0,
    nat={len(atoms)}, 
    ntyp={len(set(atoms))},
    occupations='smearing', 
    smearing='gauss', 
    degauss=1.0d-9,
    ecutwfc=50, 
    ecutrho=500,
/
&ELECTRONS
    electron_maxstep=100,
    conv_thr=1.0d-9,
    mixing_mode='plain',
    mixing_beta=0.8d0,
    diagonalization='david',
/
&IONS
    ion_dynamics='bfgs'
/
&CELL
    press_conv_thr=0.1
/
ATOMIC_SPECIES
"""

        # Add atomic species information
        for atom in set(atoms):
            atomic_mass = finder_ele_mass(atom)
            pseudo_file = finder_pseudo(self.pseudo_dir, atom)
            input_file_content += f"{atom} {atomic_mass} {pseudo_file}\n"

        # Add cell parameters
        input_file_content += "CELL_PARAMETERS (angstrom)\n"
        for row in cell_matrix:
            input_file_content += ' '.join([f"{x:.9f}" for x in row]) + '\n'

        # Add atomic positions
        input_file_content += "ATOMIC_POSITIONS (crystal)\n"
        for atom, coord in zip(atoms, fra_coords):
            input_file_content += f"{atom} " + ' '.join([f"{x:.9f}" for x in coord]) + '\n'

        # Add k-points
        input_file_content += "K_POINTS {crystal_b}\n"
        input_file_content += "6\n"
        input_file_content += "0.5 0.5 0.5 50\n"
        input_file_content += "0.0 0.0 0.0 50\n"
        input_file_content += "0.5 0.0 0.5 20\n"
        input_file_content += "0.625 0.25 0.625 0\n"
        input_file_content += "0.375 0.375 0.75 50\n"
        input_file_content += "0.0 0.0 0.0 1\n"

        # Save the input file
        with open(f"./band.in", 'w') as file:
            file.write(input_file_content)

        print(f"Input file 'band.in' is saved in the current folder")
        print(f"Use command to run the file : mpirun -n 16 {self.qe_dir}/bin/pw.x < band.in > band.out")


    def dos(self,outdir='./exQE'):
        """
        Processing dos files

        :param outdir: str, the directory of output files, default is './exQE'
        """
       
        print('Performing dos calculation...')
        input_file_content = f"""
&DOS
    prefix='caobin',
    outdir='{outdir}', 
    ngauss=1,
    degauss=1.5d-2,
    DeltaE=1.0d-2, 
    fildos='DOS.dos' 
/
"""
        
        with open(f"./dos.in", 'w') as file:
            file.write(input_file_content)

        print(f"Input file 'dos.in' is saved in the current folder")
        print(f"Use command to run the file : mpirun -n 16 {self.qe_dir}/bin/dos.x < dos.in > dos.out")

    def ph(self, cif, outdir='./exQE'):
        """
        Processing phonons templates

        :param outdir: str, the directory of output files, default is './exQE'
        """
        print('Performing phonons calculation...')

        # Convert CIF file to atomic database and primitive cell
        cov_atom = Cell().cif2db(cif, 'None')
        _, _, pr_atom = Cell().conv2prim(cov_atom)
        print('The primitive unit cell is created')

        # Extract necessary information from the primitive cell
        atoms = pr_atom.get_chemical_symbols()
        unique_atoms = set(atoms)

        input_file_content = f"""
% q-grid
&INPUTPH
    tr2_ph=1.0d-12,
    prefix='caobin',
    outdir='{outdir}', 
    fildyn='phonon.dyn',
    ldisp=.true.,
    nq1=6, nq2=6, nq3=6
"""
        # Add atomic species information
        for i, atom in enumerate(unique_atoms, start=1):
            atomic_mass = finder_ele_mass(atom)
            input_file_content += f"    amass({i}) = {atomic_mass},\n"

        input_file_content += "/\n"


        with open(f"./ph.in", 'w') as file:
            file.write(input_file_content)

        print(f"Input file 'ph.in' is saved in the current folder")
        print(f"Use command to run the file : mpirun -n 16 {self.qe_dir}/bin/ph.x < ph.in > ph.out")


    def q2r(self,):
        """
        Mapping photons (q) to real space 
        """
       
        print('Mapping photons (q) to real space...')
        input_file_content = f"""
&INPUT
    fildyn='phonon.dyn',
    zasr='simple',
    flfrc='phonon.fc'
/
"""
        
        with open(f"./qr.in", 'w') as file:
            file.write(input_file_content)

        print(f"Input file 'qr.in' is saved in the current folder")
        print(f"Use command to run the file : {self.qe_dir}/bin/q2r.x < qr.in > qr.out")


    def matdyn(self, cif,):
        """
        Processing phonons files using the given CIF file
        """
        print('Processing phonons files ...')

        # Convert CIF file to atomic database and primitive cell
        cov_atom = Cell().cif2db(cif, 'None')
        _, _, pr_atom = Cell().conv2prim(cov_atom)
        print('The primitive unit cell is created')

        # Extract necessary information from the primitive cell
        atoms = pr_atom.get_chemical_symbols()
        unique_atoms = set(atoms)

        input_file_content = f"""

&INPUTPH
    asr='simple',
    flfrc='phonon.fc',
    flfrq='phonon.freq',
    q_in_band_form=.true.,
"""
        # Add atomic species information
        for i, atom in enumerate(unique_atoms, start=1):
            atomic_mass = finder_ele_mass(atom)
            input_file_content += f"    amass({i}) = {atomic_mass},\n"

        input_file_content += "/"
        input_file_content += """
6
  gG    40
  X     10
  U      0
  K     40
  gG    40
  L     1
"""


        with open(f"./matdyn.in", 'w') as file:
            file.write(input_file_content)

        print(f"Input file 'matdyn.in' is saved in the current folder")
        print(f"Use command to run the file : mpirun -n 16 {self.qe_dir}/bin/matdyn.x < matdyn.in > matdyn.out")



############################functions####################################

def finder_ele_mass(element):
    """
    Find the atomic mass of an element.

    :param element: str, the chemical symbol of the element
    :return: float, the atomic mass of the element
    """
    # Dictionary of atomic masses (example data)
    atomic_masses = {
        'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81,
        'C': 12.01, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
        'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
        'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
        'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
        'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
        'Ga': 69.723, 'Ge': 72.63, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
        'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
        'Nb': 92.906, 'Mo': 95.95, 'Tc': 98, 'Ru': 101.07, 'Rh': 102.91,
        'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71,
        'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91,
        'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24,
        'Pm': 145, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93,
        'Dy': 162.50, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05,
        'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21,
        'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59,
        'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Th': 232.04, 'Pa': 231.04,
        'U': 238.03, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247, 'Bk': 247,
        'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 262,
        'Rf': 267, 'Db': 270, 'Sg': 271, 'Bh': 270, 'Hs': 277, 'Mt': 278,
        'Ds': 281, 'Rg': 282, 'Cn': 285, 'Nh': 286, 'Fl': 289, 'Mc': 290,
        'Lv': 293, 'Ts': 294, 'Og': 294
    }
    return atomic_masses.get(element, 0.0)

def finder_pseudo(pseudo_dir, element):
    """
    Find the pseudo potential file for an element by searching in the pseudo_dir directory.

    :param element: str, the chemical symbol of the element
    :return: str, the pseudo potential file name for the element
    """
    for root, dirs, files in os.walk(pseudo_dir):
        for file in files:
            if element in file:
                return file
    return f"{element}.pseudo.nofile"