import math
import os
import shutil
import subprocess
import warnings
from warnings import warn

import matplotlib.pyplot as plt
import mbuild as mb
import mosdef_gomc.formats.gmso_charmm_writer as mf_charmm
import mosdef_gomc.formats.gmso_gomc_conf_writer as gomc_control
import numpy as np
import pandas as pd
import unyt as u
from mbuild.utils.conversion import OPLS_to_RB
from mosdef_gomc.utils.conversion import OPLS_to_periodic
from scipy.optimize import curve_fit
from unyt.dimensions import angle, energy, length, temperature

import mosdef_dihedral_fit.utils.file_read_and_write as mdf_frw
import mosdef_dihedral_fit.utils.math_operations as mdf_math

warnings.filterwarnings("ignore")


def fit_dihedral_with_gomc(
    fit_dihedral_atom_types,
    mol2_file,
    forcefield_file,
    temperature_unyt_units,
    gomc_binary_path,
    qm_log_file_dict,
    manual_dihedral_atom_numbers_list=None,
    zero_dihedral_atom_types=None,
    qm_engine="gaussian",
    combining_rule=None,
    atom_type_naming_style="general",
    gomc_cpu_cores=1,
    r_squared_min=0.98,
    r_squared_rtol=2.5e-02,
    opls_force_k0_zero=False,
):
    """Fit the desired dihedral to a MM force field, based on QM data.

    This fits the desired dihedral to a Molecular Mechanics (MM) field,
    based on Quantum Mechanics (QM) simulation data. The fit is to the
    OPLS dihedral style, but the periodic/CHARMM dihedral and RB
    torsions are also fit up to the 4th cosine power, inline with the
    OPLS dihedral. The user has options to remove some QM data points if
    they desire.

    NOTE: The Molecular Mechanics (MM) energy calculation use GOMC, which
    does not currently calculate any impropers, even if input in the XML
    or GOMC force field file.  Therefore, systems with any impropers may
    not be fit properly or the user should use caution when using this
    function (GOMC) to do the dihedral fitting.

    NOTE: This dihedral fitting process can accommodate fitting more
    than 1 dihedral fit of the same type simultaneously.

    NOTE: Not all cos power terms are allowed to be utilized in the
    dihedral fits. This function tests if a cos power is able to be
    used for the fit, and only uses the valid power in the fitting
    process.

    NOTE: The 'extracted_gaussian_data' and 'GOMC_simulations'
    folder are deleted at the beginning of this function,
    and recreated while running this function to ensure only the
    lasted data is in these folders.

    NOTE: The OPLS dihedral equation

    .. math::
    opls_dihedral_n_1 &= 1/2 *(
        &= k0
        &= + k1 * (1 + cos(1 * phi))
        &= + k2 * (1 - cos(2 * phi))
        &= + k3 * (1 + cos(3 * phi))
        &= + k4 * (1 - cos(4 * phi))
        &= )

    Parameters
    ----------
    fit_dihedral_atom_types: list of four (4) strings (Example: ['HC', 'CT, 'CT, 'HC'])
        The atom types/classes (strings in the list) of the dihedral which is
        being fitted with non-zero k-values.

        NOTE: The extracted atom types/classes can be determined also by
        looking at the 'GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp' and
        'GOMC_pdb_psf_ff_files_dihedrals_zeroed.inp' files in the
        'GOMC_simulations' folder. These files can also be checked to
        confirm it is zeroing the correct dihedrals.
    mol2_file: str
        The mol2 file which matches the element, atom type, bonded connnections,
        the 'EXACT ATOM ORDER AND CONFIGURATION AS IN THE QM SIMULATION INPUT FILES'.
        This is required to know the MM bonding in the atoms, because QM simulations
        do not explictly specify the system bonds.
    forcefield_file: str
        Apply a foyer or gmso forcefield to the output file by selecting a
        force field XML file with its path or by using the standard force
        field name provided the `foyer` package. This force field file must
        contain all the non-bonded, bonds, angles, dihedrals parameters for
        the system.

        * Example str for FF file: 'path_to file/trappe-ua.xml'
    temperature_unyt_units: unyt.unyt_quantity
        The temperature of the system that was performed for the Quantum Mechanics
        (QM) simulation.
    gomc_binary_path: str
        The path or directory of the GOMC binary file "GOMC_CPU_NVT" (GOMC >= v2.75),
        which is used to perform the Molecular Mechanics (MM) energy calculations.
        This does not include the "GOMC_CPU_NVT" in this variable.

        Example: '/home/brad/Programs/GOMC/GOMC_2_76/bin'
    qm_log_file_dict: dict, {str: [int>=0, ..., int>=0]}
        * qm_engine="gaussian"
            This is a dictionary comprised of a key (string) of the QM log file path and name
            (Gaussian 16 log file only), and a list of integers, which are the QM optimization
            parameters to remove from the written data, in order of reading from each file.
            These can be seen in the order of the dictionary file name (strings).
            These removed parameters allow users to remove any bad or repeated data
            points for the QM log file when needed.

            Example 1: {'path/gaussian_log_file.log': []}

            Uses all the optimized data points from the 'path/gaussian_log_file.log' file.

            Example 2: {'path/gaussian_log_file.log': [0, 23]}
            Uses all data points from the 'path/gaussian_log_file.log' file, except points
            0 and 23.  NOTE: Python counting starts at 0.

        * qm_engine="gaussian_style_final_files"
            This is a dictionary comprised of a key (string) of the  file paths to the
            Gaussian 16 style final formatted files, and a list of integers, which are the
            QM optimization parameters to remove from the written data, in order of reading
            from each folder. These can be seen in the order of the dictionary file name (strings).
            These removed parameters allow users to remove any bad or repeated data points
            for the QM log file when needed.
            NOTE: The energy and dihedral angle file in this directory need to be
            named 'dihedral.txt' for the energy and dihedral angle values (one 1 per directory).

            Example of energy and dihedral angle file ('dihedral.txt'):

            | # Scan of Total Energy
            | # X-Axis:  Scan Coordinate
            | # Y-Axis:  Total Energy (Hartree)
            | #                  X                   Y
            |                  0.0     -267.0062955742
            |                 10.0     -267.0062900424

            NOTE: The coordinate files in this directory need to be
            named 'dihedral_coords_position_XXXX.txt' for the each angles coordinate values.
            There are as XXX file in this directory where XXX is the number of dihedral angles.
            The file numbering starts at 1 so the files are named 'dihedral_coords_position_1.txt'
            to 'dihedral_coords_position_XXXX.txt'

            Example of coordinate file ('dihedral_coords_position_1.txt'):

            | Row	Highlight	Display	Tag	Symbol	X	Y	Z
            | 1       No      Show    1       C       0.077153        -0.010211       0.106889
            | 2       No      Show    2       C       -1.455163       0.076994        0.364648
            | 3       No      Show    3       C       -2.162794       1.205823        -0.378912
            | 4       No      Show    4       O       0.614863        1.022719        -0.303596
            | 5       No      Show    5       O       0.581656        -1.105138       0.370604
            | 6       No      Show    6       H       -1.703737       2.157201        -0.140757
            | 7       No      Show    7       H       -2.079381       1.073202        -1.454515
            | 8       No      Show    8       H       -1.898266       -0.885627       0.121028
            | 9       No      Show    9       H       -1.593015       0.205080        1.439694
            | 10      No      Show    10      H       -3.224767       1.255506        -0.130085

            Example 1: {'path_to_gaussian_style_final_files': []}
            Uses all the optimized data points from the 'path/gaussian_log_file.log' file.

            Example 2: {'path_to_gaussian_style_final_files': [0, 23]}
            Uses all data points from the 'path/gaussian_log_file.log' file, except points
            0 and 23.  NOTE: Python counting starts at 0.

    manual_dihedral_atom_numbers_list: list of 4 integers, default=None
        NOTE: Only needed for qm_engine="gaussian_style_final_files"

        This is a list of the dihedral atom numbers in order that were used for the dihedral
        fit. This information needs to be correct and in order to produce correct results.
        The values must be the same in all the combined files.
    zero_dihedral_atom_types: nest list with lists of four (4) strings, default=None
        The nests list(s) of the other dihedrals, that need to have their k-values zeroed to
        properly fit the the 'fit_dihedral_atom_types' dihedral.

        Example: [['CT', 'CT, 'CT, 'HC'], ['NT', 'CT, 'CT, 'HC']]

        NOTE: The extracted atom types/classes can be determined also by
        looking at the 'GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp' and
        'GOMC_pdb_psf_ff_files_dihedrals_zeroed.inp' files in the
        'GOMC_simulations' folder. These files can also be checked to
        confirm it is zeroing the correct dihedrals.
    qm_engine: str ('gaussian' or 'gaussian_style_final_files'), default='gaussian'
        The Quantum Mechanics (QM) simulation engine utilized to produce the files listed
        in the 'qm_log_file_dict' variable(s).

        'gaussian' is used for the QM Gaussian 16 log files.

        'gaussian_style_final_files' are the file paths to the Gaussian 16 style final
        formatted files, which are the QM optimization parameters, in order of reading
        from each folder.
        NOTE: The energy and dihedral angle file in this directory need to be
        named 'dihedral.txt' for the energy and dihedral angle values (one 1 per directory).

    combining_rule: str ('geometric' or 'lorentz'), default = None
        The Lennard-Jones or VDW sigmas 'combining_rule'. in the foyer or GMSO XML file.

        The 'combining_rule' in the foyer or GMSO XML file for the Lennard-Jones or VDW sigmas
        (Rmin --> sigmas for Exp6). If this is None, it will use whatever is specified in
        the XML file, or the default foyer or GMSO values. BEWARE, if it is not specified
        XML file, it has a default. If this is None, it will use whatever is specified in
        the XML file, or the default foyer or GMSO values.

        'geometric' is the geometric mean used to combine the
        Lennard-Jones or VDW sigmas (Rmin --> sigmas for Exp6), as required by OPLS force field.

        'lorentz' is the arithmetic mean used to combine the
        Lennard-Jones or VDW sigmas (Rmin --> sigmas for Exp6), as required by
        TraPPE, Amber, or CHARMM force fields.

        'None', the default setting used to combine the Lennard-Jones or VDW sigmas
        (Rmin --> sigmas for Exp6), is pulled from the force field XML, but if it is
        not present it defaults to 'geometric' via MoSDeF's default settings.

        NOTE: The mixing rules for the other non-Lennard-Jones or non-VDW sigmas forces are
        listed below for reference, which applicable in GOMC and here:

        "NOTE: epsilon_ij = 'geometric' for All FFs  --> epsilon_ij = (epsilon_ii * epsilon_jj)**0.5."

        "NOTE: n_ij = 'lorentz' for Mie FFs --> n_ij = (n_ii + n_jj)/2."

        "NOTE: alpha_ij = 'geometric' for Exp FFs --> alpha_ij = (alpha_ii * alpha_jj)**0.5."


    atom_type_naming_style: str, optional, default='all_unique', ('general' or 'all_unique')
        * 'general'

        NOTE: In this case, as long as the force field XML file is correct, we can use the
        'general' convention safely, because there is no potential of atom type/class overlap
        with another force field file.

        The 'general' convention only tests if the sigma, epsilons, mass, and Mie-n values are
        identical between the different molecules (residues in this context) and their applied
        force fields and DOES NOT check that any or all of the bonded parameters have the same
        or conflicting values.

        The 'general' convention is where all the atom classes in the Foyer force field
        XML files are converted to the CHARMM-style atom types (FOYER ATOM CLASSES).
        The 'general' convention ONLY auto-checks that the sigma, epsilon, mass, and Mie-n values
        are the same and does not currently ensure all the bonded parameters are the same
        or conflicting between different force field XML files.
        If the sigma, epsilons, mass, and Mie-n values are the same between force fields
        the general method can be applied; if not, it defaults to the 'all_unique' method.

        --- Example of CHARMM style atom types in an all-atom ethane and ethanol system:

        --- Ethane: alkane carbon = CT, alkane hydrogen = HC

        --- Ethanol: alkane carbon = CT, alkane hydrogen = HC , oxygen in alcohol = OH, hydrogen in alcohol = OH

        This is only permitted when the following is true; otherwise it will default to the the 'all_unique':

        --- All the MoSDeF force field XML's atom classes' non-bonded parameters
        (sigma, epsilon, mass, and Mie-n power constant) values are THE SAME.

        --- If the general CHARMM style atom type in any residue/molecule's gomc_fix_bonds_angles,
        gomc_fix_bonds, or gomc_fix_angles NOT IN any other residue/molecule, the 'all_unique' type
        will be used.

        * 'all_unique'

        The 'all_unique' convention is the SAFE way to parameterize the system.
        The MoSDeF force field XML atom names within residue/molecule are all unique,
        where each one adds an alpha numberic value after the MoSDeF force field XML atom classes to
        ensure uniqueness within the within residue/molecule.
        The OPLS atom types/name do not require all the sigma, epsilon, mass values to be the same,
        but have less bonded class parameters.

        Example of CHARMM style atom types in an all-atom ethane and ethanol system:

        --- Ethane: alkane carbon type 0 = CT0, alkane hydrogen type 0 = HC0

        --- Ethanol: alkane carbon type 1 = CT1, alkane carbon type 2 = CT2,
        alkane hydrogen type 1 = HC1 , oxygen in alcohol type 0 = OH0, hydrogen in alcohol type 0 = OH0

        This is selected when auto-selected when:

        --- All the MoSDeF force field XML's atom classes' non-bonded parameters
        (sigma, epsilon, mass, and Mie-n power constant) values are NOT THE SAME.

        --- If the general CHARMM style atom type in any residue/molecule's gomc_fix_bonds_angles,
        gomc_fix_bonds, or gomc_fix_angles are IN any other residue/molecule.
    gomc_cpu_cores: int>0, default=1
        The number of CPU-cores that are used to perform the GOMC simulations, required
        for the Molecular Mechanics (MM) energy calulations.
    r_squared_min: float (0 < r_squared_min < 1), default=0.98
        The minimum R**2 (R-squared) value to test the validity of the fit with the
        new dihedral fitted constants, as fitted in the
        QM - MM energy data vs. the dihedral function fit, mentioned below.

        This mean that any R**2 (R-squared) value in the
        'opls_all_summed_dihedrals_k_constants_figure.pdf' plot less than this value
        will not be check for accuracy because the  R**2 (R-squared) values are compared
        differently as a check.  These are compared between the following calculations:

        QM - MM energy data vs. the dihedral function fit:
        For the MM calculations, the 'fit_dihedral_atom_types' and
        'zero_dihedral_atom_types' are dihedral energies are set to zero, which is
        during the fitting process with 1 or more of the same dihedrals being fit simultaneously.

        QM vs. the MM energy data:
        For the MM calculations, the 'fit_dihedral_atom_types' are set to the values which were
        fit for the specific cosine combinations during the fitting process with 1 or more of the
        same dihedrals being fit simultaneously, and the 'zero_dihedral_atom_types' are
        dihedral energies are set to zero.

        NOTE: This value may need adjusted to get the dihedral fit to solve correctly.
    r_squared_rtol: float (0 < r_squared_min < 1), default=2.5e-02
        Where the QM data is defined as the actual data; this is the difference
        of the dihedral's calculated R-squared values between:
        * The QM-MM fitting process, where the fit MM dihedral k-values are zero (0).
        * The MM calculations where the fit k-value are entered in the MM data and
        compared to the QM data.

        NOTE: This value may need adjusted to get the dihedral fit to solve correctly.
    opls_force_k0_zero: bool, default=False
        The k0 constant is from the equation listed above.
        If True, this is force sets the k0 constant in the opls equation to zero (0),
        which is the original OPLS form. This means that the dihedral energy fit must
        be zero (0) at dihedral angles of (-180 and 180 degrees), which could mean
        the dihedral produces both positive (+) and negative energy values (-).
        NOTE: Using this option may not allow some dihedrals to be fit properly.

        If False, the k0 constant is allowed to be fitted to a non-zero constant.
        This allows the k0 constant to be zero (0) or a scaler value. In this case,
        the dihedral's k0 constant is shifted to make the minimum of the dihedral
        equation equal to zero (0). This is the more standard overall dihedral form
        (i.e. not forcing the equation k0 constant to be zero) as it allows more
        dihedrals to be fit properly, and this energy shifting is in almost all newer
        dihedral forms now.

    Outputs
    -------
    Files containing the following information in the following relative locations:

    GOMC_simulations/GOMC_pdb_psf_ff_files.pdb
        The PDB file used in the GOMC simulations generated via MoSDeF-GOMC.
    GOMC_simulations/GOMC_pdb_psf_ff_files.pdf
        The PSF file used in the GOMC simulations generated via MoSDeF-GOMC.
    GOMC_simulations/GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp
        The original force field file generated via MoSDeF-GOMC with the
        force field parameters exactly as listed in the provided force field
        XML file.
    GOMC_simulations/GOMC_pdb_psf_ff_files_dihedrals_zeroed.inp
        The modified force field file with the fit and other selected dihedrals
        zeroed out.
    GOMC_simulations/GOMC_pdb_psf_ff_files_OPLS_fit_XXX_dihedral.inp
        The modified force field file with the fit is set to the fitted k-values
        ('fit_dihedral_atom_types') with all the other selected dihedrals
        ('zero_dihedral_atom_types') being zeroed out, where the XXX
        in the file name being the different  possibilities of cos power combinations
        of the dihedral fit. These may be 1 or more of these files/cos power
        combinations.
    GOMC_simulations/GOMC_zeroed_dihedral_coords_XXX.conf
        The GOMC control files for all the phi dihedral angles selected
        from the QM simulations where the 'fit_dihedral_atom_types' and the
        'zero_dihedral_atom_types' k-values are all set to zero.  The XXX is the
        integer number of the dihedrals starting at 1, adding 1 for every
        additional phi dihedral angles selected from the QM simulations.

        These are simulations are to get the total energy difference between
        QM and MM simulation with all the 'fit_dihedral_atom_types' and the
        'zero_dihedral_atom_types' k-values are all set to zero.
    GOMC_simulations/GOMC_OPLS_fit_YYY_dihedral_coords_XXX.conf
        The GOMC control files for all the phi dihedral angles selected
        from the QM simulations where the 'fit_dihedral_atom_types' are set to
        the solved k-values for the cos power equation combination, which is listed
        as the variable YYY, and the 'zero_dihedral_atom_types' k-values are all set
        to zero.  The XXX is the integer number of the dihedrals starting at 1, adding
        1 for every additional phi dihedral angles selected from the QM simulations.

        These are simulations are to get the total energy difference between
        QM and MM simulation with the 'fit_dihedral_atom_types' set to
        the solved k-values for each YYY combination of cos powers.  There
        will be a set of phi angles for every YYY combination.

        These GOMC simulations are used to validate the fit when the k-values
        for each individual dihedral is automatically entered in the force
        field file. NOTE: In the original fitting process, there can be more
        than 1 dihedral of the same type fit simultaneously.
    GOMC_simulations/output_GOMC_zeroed_dihedral_coords_XXX.txt
        The log output for the 'GOMC_zeroed_dihedral_coords_XXX.conf' simulation.
        The XXX is the integer number of the dihedrals starting at 1, adding 1 for
        every additional phi dihedral angles selected from the QM simulations.
    GOMC_simulations/output_GOMC_OPLS_fit_YYY_dihedral_coords_XXX.txt
        The log output for the 'GOMC_OPLS_fit_YYY_dihedral_coords_XXX.conf' simulation.
        The XXX is the integer number of the dihedrals starting at 1, adding 1 for
        every additional phi dihedral angles selected from the QM simulations.
        The variable YYY is the cos power equation combinations used for the
        'fit_dihedral_atom_types' to fit the k-values.
    extracted_gaussian_data/dihedral.txt
        The QM data in a Gaussian-style output file which holds the scanned
        dihedral angles, in degrees, and the optimized energy value, in Hartree units,
        for the molecule/fragment.
    extracted_gaussian_data/dihedral_coords_position_XXX.txt
        The optimized QM dihedral coordinates in a Gaussian-style output file.
        The XXX is the integer number of the dihedrals starting at 1, adding 1 for
        every additional phi dihedral angles selected from the QM simulations.
        The coordinates are in Angstroms.
    xyz_restart_xsc_coor_files/dihedral_coords_position_XXX.xyz
        The optimized QM dihedral coordinates in the '.xyz' format. This permits
        VMD to convert the '.xyz' format to a GOMC-usable '.coor' format, which
        contains all the coordinate precision of the QM data and allows GOMC to
        restart the simulation with all this precision. Otherwise, GOMC would
        need to use the 3-decimal precision of the PDB format.
        The XXX is the integer number of the dihedrals starting at 1, adding 1 for
        every additional phi dihedral angles selected from the QM simulations.
        The coordinates are in Angstroms.
    xyz_restart_xsc_coor_files/dihedral_coords_position_XXX.coor
        The GOMC-usable '.coor' format, allowing the optimized QM dihedral coordinates
        to be used in full percision when restarting the GOMC simulation.
        Otherwise, GOMC would need to use the 3-decimal precision of the PDB format.
        The XXX is the integer number of the dihedrals starting at 1, adding 1 for
        every additional phi dihedral angles selected from the QM simulations.
        The coordinates are in Angstroms.
    all_normalized_energies_in_kcal_per_mol.txt
        For each dihedral angle (degrees), the QM, MM, and difference in dihedral
        energies (kcal/mol) when the 'fit_dihedral_atom_types' and the
        'zero_dihedral_atom_types' k-values are all set to zero. This also
        contains the sum of all the C1, C2, C3, and C4 values for every
        'fit_dihedral_atom_types' in the system, where there may be multiple of the
        same dihedral in the same molecule.

        This file contains the corrected dihedral fit points between QM and MM
        after all the non-bonded, 1-4 dihedral scaling, and
        bonded (no MM dihedral bonded interactions) interactions of the
        MM force fields are taken into account.  Specifically, the
        'Gaussian_minus_GOMC_E_kcal_per_mol' column contains these points,
        which can be used to perform alternative dihedral fits outside of
        this software package.  For example, higher cosine power fits or
        the CHARMM style periodic fits with no C0/k0 term (CHARMM has the
        C0/k0 term being a harmonic dihedral, not a periodic dihedral).

        OPLS dihedral equation in with C1, C2, C2, and C4 instead of the cos terms:

        OPLS_energy = k1 * C1 + k2 * C2 + k3 * C3 + k4 * C4

        C1 = (1 + cos(1 * phi))

        C2 = (1 - cos(2 * phi))

        C3 = (1 + cos(3 * phi))

        C4 = (1 - cos(4 * phi))

    all_normalized_energies_OPLS_fit_YYY_in_kcal_per_mol.txt
        For each dihedral angle (degrees), the QM, MM, and difference in dihedral
        energies (kcal/mol) when the 'fit_dihedral_atom_types' are set to
        the solved k-values for the cos power equation combination, which is listed
        as the variable YYY. This also contains fit k1, k2, k3, and k4
        values (kcal/mol), are the same in every dihedral angle row.

        These files contains the corrected dihedral fit points between QM and MM
        after all the non-bonded, 1-4 dihedral scaling, and
        bonded (including the MM fitted dihedral YYY energies interactions)
        interactions of the MM force fields are taken into account.
        The 'Gaussian_minus_GOMC_E_kcal_per_mol' column contains these
        points, which can be used to determine the goodness of the fit; if all
        these values are zero (0) or close to zero (0), it means the fit is good.

        OPLS dihedral equation in with C1, C2, C2, and C4 instead of the cos terms:

        OPLS_energy = k1 * C1 + k2 * C2 + k3 * C3 + k4 * C4

        C1 = (1 + cos(1 * phi))

        C2 = (1 - cos(2 * phi))

        C3 = (1 + cos(3 * phi))

        C4 = (1 - cos(4 * phi))

    gomc_raw_energies_in_Kelvin.txt
        This is the initial raw energy, in Kelvin, extracted from all the
        'output_GOMC_zeroed_dihedral_coord_XXX.txt' files.
    all_normalized_energies_OPLS_fit_YYY_in_kcal_per_mol.txt
        This is the initial raw energy, in Kelvin, extracted from all the
        'output_GOMC_OPLS_fit_YYY_dihedral_coord_XXX.txt' files, where YYY
        is an additional output file for each OPLS cosine combination.
    opls_dihedral_k_constants_fit_energy.txt
        The OPLS dihedral style k-value solutions given the 'fit_dihedral_atom_types'
        input for each valid cosine power combination.  There could be 1 or
        more valid cosine power combinations.
        The k-values are in kcal/mol energy units.
    periodic_dihedral_k_constants_fit_energy.txt
        The periodic/CHARMM dihedral style k-value solutions given the
        'fit_dihedral_atom_types' input for each valid cosine power combination.
        There could be 1 or more valid cosine power combinations.
        The periodic/CHARMM style dihedral k-value were analytically converted
        from the OPLS dihedral form.
        The k-values are in kcal/mol energy units.
    RB_torsion_k_constants_fit_energy.txt
        The RB torsion style k-value solutions given the
        'fit_dihedral_atom_types' input for each valid cosine power combination.
        There could be 1 or more valid cosine power combinations.
        The RB torsion's k-value were analytically converted
        from the OPLS dihedral form.
        The k-values are in kcal/mol energy units.
    opls_all_summed_dihedrals_k_constants_figure.pdf
        For each valid cosine power combinations, the total OPLS dihedral energies
        between QM and MM simulations are plotted, where the MM simulation's fitted
        k-values ('fit_dihedral_atom_types') and all the other selected dihedrals
        ('zero_dihedral_atom_types') are zeroed out.
        This is the total energy difference in QM vs MM energy
        for all the 'fit_dihedral_atom_types' summed together, as there can
        be 1 or multiple of this dihedral type in a molecule, which were
        all fit simultaneously.
    opls_all_single_fit_dihedral_k_constants_figure.pdf
        The final OPLS dihedral fit for the individual dihedrals
        (individual dihedrals in the 'fit_dihedral_atom_types') for each
        each valid cosine power combination. These plotted dihedral energies would be
        the energy from a single (1) dihedral in the system, or the dihedral energy
        (via k-values) that would be entered into the force field parameters.
        This plot allows the users to compare the different valid cosine power
        combinations and their R-squared values.
        If multiple R-squared values are both nearly perfect and the plots
        look different, this may be a sign that something is not correct with 1
        or more of the nearly perfect R-squared fitted values, or fitting procedure
        itself.
    """
    # check if 'mol2_file' file is correct format
    if not isinstance(mol2_file, str):
        raise TypeError(
            "ERROR: Please enter mol2 file ('mol2_file') as a string."
        )

    extension_ff_name = os.path.splitext(mol2_file)[-1]
    if extension_ff_name != ".mol2":
        raise ValueError(
            "ERROR: Please enter enter mol2 file ('mol2_file') name with the .mol2 extension."
        )

    if not os.path.exists(mol2_file):
        raise ValueError(
            f"ERROR: The {mol2_file} file ('mol2_file') does not exists."
        )

    # check if 'forcefield_file' file is correct format
    if not isinstance(forcefield_file, str):
        raise TypeError(
            "ERROR: Please enter xml file ('forcefield_file') as a string."
        )

    extension_ff_name = os.path.splitext(forcefield_file)[-1]
    if extension_ff_name != ".xml":
        raise ValueError(
            "ERROR: Please enter enter xml file ('forcefield_file') name with the .xml extension."
        )

    if not os.path.exists(forcefield_file):
        raise ValueError(
            f"ERROR: The {forcefield_file} file ('forcefield_file') does not exists."
        )

    if (
        qm_engine == "gaussian"
        and manual_dihedral_atom_numbers_list is not None
    ):
        warn(
            "WARNING: When reading the qm_engine = 'gaussian' files, the "
            "'manual_dihedral_atom_numbers_list' is set to None, and will not be used, "
            "because the the gaussian log files already contain this information."
        )
        manual_dihedral_atom_numbers_list = None

    # test the temperature_unyt_units input
    print_error_value = f"ERROR: The 'temperature_unyt_units' is not temperature of type {type(u.unyt_quantity)}."
    if isinstance(temperature_unyt_units, u.unyt_quantity):
        if temperature == temperature_unyt_units.units.dimensions:
            temperature_unyt_units = temperature_unyt_units.to("K")

        else:
            raise ValueError(print_error_value)

    else:
        raise TypeError(print_error_value)

    # test the qm_log_file_dict input
    print_error_value = (
        "ERROR: The 'qm_log_file_dict' is not a dict "
        "with a string keys and list of int>=0 as the values. Example: "
        "{'path/HC_CT_CT_HC_part_1.log'): [], 'path/HC_CT_CT_HC_part_2.log'): [0, 5]}"
    )
    if isinstance(qm_log_file_dict, dict):
        for key_j, value_j in qm_log_file_dict.items():
            if isinstance(key_j, str):
                if not os.path.exists(key_j):
                    raise ValueError(
                        f"ERROR: The {key_j} file ('qm_log_file_dict') does not exists."
                    )

            else:
                raise TypeError(print_error_value)

            if isinstance(value_j, list):
                for int_j in value_j:
                    if not isinstance(int_j, int) or int_j < 0:
                        raise TypeError(print_error_value)

            else:
                raise TypeError(print_error_value)

    else:
        raise TypeError(print_error_value)

    # check if 'gomc_binary_path' leads to the file is correct format GOMC_CPU_NVT
    if not isinstance(gomc_binary_path, str):
        raise TypeError(
            "ERROR: Please enter the 'gomc_binary_path' file as a string."
        )

    if not os.path.exists(f"{gomc_binary_path}/{'GOMC_CPU_NVT'}"):
        raise ValueError(
            f"ERROR: The 'gomc_binary_path' file does not exist or contain the GOMC 'GOMC_CPU_NVT' file."
        )

    # test the 'zero_dihedral_atom_types' input
    print_error_value = (
        "ERROR: The 'zero_dihedral_atom_types' is not None or a list containing "
        "lists with 4 strings each. Example: "
        "[['CT', 'CT, 'CT, 'HC'], ['NT', 'CT, 'CT, 'HC']]."
    )
    if isinstance(zero_dihedral_atom_types, (list, type(None))):
        if isinstance(zero_dihedral_atom_types, list):
            for list_j in zero_dihedral_atom_types:
                if isinstance(list_j, list) and len(list_j) == 4:
                    for str_j in list_j:
                        if not isinstance(str_j, str):
                            raise TypeError(print_error_value)

                else:
                    raise TypeError(print_error_value)

    else:
        raise TypeError(print_error_value)

    # Check if 'combining_rule' is correct
    # Set the VDWGeometricSigma variable from the combining_rule
    if combining_rule == None:
        VDWGeometricSigma = None
    elif combining_rule == "geometric":
        VDWGeometricSigma = True
    elif combining_rule == "lorentz":
        VDWGeometricSigma = False
    else:
        raise ValueError(
            "ERROR: Please enter the 'combining_rule' file as a string of 'geometric' or 'lorentz' or None."
        )

    # test the 'atom_type_naming_style' input
    if isinstance(atom_type_naming_style, str):
        if not atom_type_naming_style in ["general", "all_unique"]:
            raise ValueError(
                f"ERROR: The 'atom_type_naming_style' = {atom_type_naming_style}, which is not "
                f"any of the available options. "
                f"The options are 'general' or 'all_unique'."
            )

    else:
        raise TypeError(
            f"ERROR: The 'atom_type_naming_style' is a {type(atom_type_naming_style)}, but it needs to be a str."
        )

    # test the 'r_squared_min' input
    if isinstance(r_squared_min, float):
        if not (r_squared_min > 0 and r_squared_min < 1):
            raise ValueError(
                f"ERROR: The 'r_squared_min'= {r_squared_min}, "
                f"but it must be a 0<float<1."
            )

    else:
        raise TypeError(
            f"ERROR: The 'r_squared_min' is a {type(r_squared_min)}, "
            f"but it must be a 0<float<1."
        )

    # test the 'r_squared_rtol' input
    if isinstance(r_squared_rtol, float):
        if not (r_squared_rtol > 0 and r_squared_rtol < 1):
            raise ValueError(
                f"ERROR: The 'r_squared_rtol' = {r_squared_rtol}, "
                f"but it must be a 0<float<1."
            )

    else:
        raise TypeError(
            f"ERROR: The 'r_squared_rtol' is a {type( r_squared_rtol)}, "
            f"but it must be a 0<float<1."
        )

    # test the 'gomc_cpu_cores' input
    if isinstance(gomc_cpu_cores, int):
        if gomc_cpu_cores <= 0:
            raise ValueError(
                f"ERROR: The 'gomc_cpu_cores' = {gomc_cpu_cores}, and it must be an int > 0."
            )

    else:
        raise TypeError(
            f"ERROR: The 'gomc_cpu_cores' is a {type(gomc_cpu_cores)}, but it needs to be a int."
        )

    # check values write the qm data files data out
    if isinstance(qm_engine, str):
        if qm_engine == "gaussian":
            mdf_frw.write_qm_data_files(
                qm_log_file_dict,
                manual_dihedral_atom_numbers_list=manual_dihedral_atom_numbers_list,
                qm_engine=qm_engine,
            )
        elif qm_engine == "gaussian_style_final_files":
            mdf_frw.write_qm_data_files(
                qm_log_file_dict,
                manual_dihedral_atom_numbers_list=manual_dihedral_atom_numbers_list,
                qm_engine=qm_engine,
            )

        else:
            raise ValueError(
                f"ERROR: The 'qm_engine' = {qm_engine}, which is not "
                f"any of the available options. "
                f"The options are 'gaussian' or 'gaussian_style_final_files'."
            )

    else:
        raise TypeError(
            f"ERROR: The 'qm_engine' is a {type(qm_engine)}, but it needs to be a str."
        )

    # check if opls_force_k0_zero' is a bool
    if not isinstance(opls_force_k0_zero, bool):
        raise TypeError(
            "ERROR: Please enter the 'opls_force_k0_zero' as a bool."
        )

    # **************************************************************
    # **************************************************************
    # Use the existing Gaussian (QM) file for their energy data.
    # Build the system using a mol2 file via mBuild and MoSDeF.
    # Use GOMC to run the MoSDeF build with the .coor and .xsc files
    # to yield the exact configuration of the QM equilibrium fragments.
    # Note: the MoSDeF XML file must set the target dihedral to zero (0) energy.
    # Subract the QM - GOMC eneries, normalize the dihedral to min of zero (0) energy,
    # fit the various dihedral combinations to standard OPLS form (no f0 constant),
    # then convert each fit to the periodic/CHARMM dihedral and RB torsions forms.
    # (START)
    # **************************************************************
    # **************************************************************

    # **************************************************************
    # make the PDB, PSF and FF files for all the dihedral angles (START)
    # **************************************************************

    # delete existing make a gomc simulation folder ('GOMC_simulations') and move there
    gomc_runs_folder_name = "GOMC_simulations"

    if os.path.isdir(gomc_runs_folder_name):
        shutil.rmtree(gomc_runs_folder_name)
    os.mkdir(gomc_runs_folder_name)

    # The gomc raw energy filename in Kelvin Energies
    gomc_raw_energy_filename = "gomc_raw_energies_in_Kelvin.txt"

    # The combined GOMC and Gaussian dihedral angles and energies in kcal/mol
    gomc_gaussian_kcal_per_mol_energy_filename = (
        "all_normalized_energies_in_kcal_per_mol.txt"
    )

    output_gomc_pdb_psf_ff_file_name_str = f"GOMC_pdb_psf_ff_files"
    seed_no = 12345

    # load bis(ethylhexylamido) with Gd
    fragment = mb.load(mol2_file, smiles=False)
    fragment.name = "TMP"

    residues_list = [fragment.name]
    fix_residues = None
    fix_residue_in_box = None

    # Build the methanol liquid box_1

    # liquid_box_0_length_nm must be a value <= 999.8 nm and an interger value in angstroms <= 9998 Ang
    liquid_box_0_length_nm = 999.8

    # Started building the fragment for the GOMC simulation with dihedral k values = 0
    box_0_liq = mb.fill_box(
        compound=[fragment],
        n_compounds=[1],
        box=[
            liquid_box_0_length_nm,
            liquid_box_0_length_nm,
            liquid_box_0_length_nm,
        ],
        seed=seed_no,
    )

    # Started building the Charmm object for the GOMC simulation with dihedral k values = 0
    # build the charmm object
    charmm = mf_charmm.Charmm(
        box_0_liq,
        f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}",
        structure_box_1=None,
        filename_box_1=None,
        ff_filename=f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_per_xml",
        forcefield_selection=forcefield_file,
        residues=residues_list,
        bead_to_atom_name_dict=None,
        fix_residue=fix_residues,
        fix_residue_in_box=fix_residue_in_box,
        gomc_fix_bonds_angles=None,
        atom_type_naming_style=atom_type_naming_style,
    )

    # Started writing the PDB, PSF, and FF files for the GOMC simulation with dihedral k values = 0
    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()

    # **************************************************************
    # make the PDB, PSF and FF files for all the dihedral angles (END)
    # **************************************************************

    # **************************************************************
    # make the GOMC control file with the selected dihedral angles set to zero (START)
    # **************************************************************
    mdf_frw.change_gomc_ff_file_dihedral_values(
        f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_per_xml.inp",
        f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_zeroed.inp",
        fit_dihedral_atom_types,
        fit_dihedral_opls_k_0_1_2_3_4_values=[0, 0, 0, 0, 0],
        zero_dihedral_atom_types=zero_dihedral_atom_types,
    )
    # **************************************************************
    # make the GOMC control file with the selected dihedral angles set to zero (END)
    # **************************************************************

    # **************************************************************
    # Create the xyz files in a folder so VMD can read them with
    # the correct precision, which will be converted to a .coor
    # restart file so GOMC can do the MM energy calculations with
    # the correct precision.
    # (START)
    # **************************************************************

    # delete the 'xyz_and_coor_files' folder, if it exists, and create a new 'xyz_and_coor_files' folder
    xyz_xsc_coor_files_directory = "xyz_restart_xsc_coor_files"
    if os.path.isdir(xyz_xsc_coor_files_directory):
        shutil.rmtree(xyz_xsc_coor_files_directory)
    os.mkdir(xyz_xsc_coor_files_directory)

    # write all the xyz coordinate from the gaussian optimized coordinate file in the 'xyz_files' folder
    [
        atom_pdb_names_list,
        elementpdb_names_list,
    ] = mdf_frw.get_atom_names_and_elements_from_pdb(
        f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}.pdb"
    )

    qm_energy_file_dir_and_name = "extracted_gaussian_data/dihedral.txt"
    qm_parital_coordinate_file_starting_dir_and_name = (
        "extracted_gaussian_data/dihedral_coords_position_"
    )
    qm_coordinate_file_extension = "txt"

    # check the gaussian file is correct
    mdf_frw.check_gaussian_angle_energy_file_correct(
        qm_energy_file_dir_and_name
    )

    # Read the gaussian data and extract angles and number of scans (number of angles and degress analyzed)
    gaussian_raw_degrees_list = (
        pd.DataFrame(
            pd.read_csv(qm_energy_file_dir_and_name, sep="\s+", header=3)
        )
        .iloc[:, 0]
        .tolist()
    )
    total_qm_scans = len(gaussian_raw_degrees_list)

    mdf_frw.write_xyz_file_from_gaussian_coordinates(
        elementpdb_names_list,
        qm_parital_coordinate_file_starting_dir_and_name,
        qm_coordinate_file_extension,
        xyz_xsc_coor_files_directory,
        total_qm_scans,
    )

    # Using vmd write the GOMC restart .coor files required for the
    mdf_frw.write_restart_coor_from_xyz_file(
        xyz_xsc_coor_files_directory, total_qm_scans
    )

    # **************************************************************
    # Create the xyz files in a folder so VMD can read them with
    # the correct precision, which will be converted to a .coor
    # restart file so GOMC can do the MM energy calculations with
    # the correct precision.
    # (END)
    # **************************************************************

    # **************************************************************
    # make the GOMC control file for each dihedral angle (START)
    # **************************************************************

    # write the GOMC control files
    for scan_iter in range(1, len(gaussian_raw_degrees_list) + 1):
        read_gomc_restart_file_coor_dir_and_name = f"../{xyz_xsc_coor_files_directory}/dihedral_coords_position_{scan_iter}.coor"
        read_gomc_restart_file_xsc_dir_and_name = (
            f"../{xyz_xsc_coor_files_directory}/starting_point.xsc"
        )

        control_file_name_str = f"GOMC_zeroed_dihedral_coords_{scan_iter}.conf"
        output_name_control_file_name_str = (
            f"output_GOMC_zeroed_dihedral_coords_{scan_iter}.txt"
        )

        # cutoff and tail correction
        Rcut = int(liquid_box_0_length_nm * 10 / 2) * u.angstrom
        RcutLow = 0 * u.angstrom
        LRC = False
        Exclude = "1-3"

        # MC move ratios
        DisFreq = 0.5
        VolFreq = 0
        RotFreq = 0.5
        RegrowthFreq = 0
        IntraSwapFreq = 0
        CrankShaftFreq = 0
        SwapFreq = 0

        # NVT GOMC control file writing for the GOMC simulation with dihedral k values = 0
        # calc MC steps for gomc equilb
        MC_steps = 2
        EqSteps = 1
        AdjSteps = 1
        output_true_list_input = [True, 1]
        output_false_list_input = [False, 1]

        gomc_control.write_gomc_control_file(
            charmm,
            f"{gomc_runs_folder_name}/{control_file_name_str}",
            "NVT",
            MC_steps,
            temperature_unyt_units,
            ff_psf_pdb_file_directory=None,
            check_input_files_exist=False,
            Parameters=f"{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_zeroed.inp",
            Restart=True,
            Checkpoint=False,
            ExpertMode=False,
            Coordinates_box_0=f"{output_gomc_pdb_psf_ff_file_name_str}.pdb",
            Structure_box_0=f"{output_gomc_pdb_psf_ff_file_name_str}.psf",
            binCoordinates_box_0=read_gomc_restart_file_coor_dir_and_name,
            extendedSystem_box_0=read_gomc_restart_file_xsc_dir_and_name,
            binVelocities_box_0=None,
            Coordinates_box_1=None,
            Structure_box_1=None,
            binCoordinates_box_1=None,
            extendedSystem_box_1=None,
            binVelocities_box_1=None,
            input_variables_dict={
                "PRNG": seed_no,
                "Pressure": None,
                "Ewald": True,
                "ElectroStatic": True,
                "VDWGeometricSigma": VDWGeometricSigma,
                "Rcut": Rcut,
                "RcutLow": RcutLow,
                "LRC": LRC,
                "Exclude": Exclude,
                "DisFreq": DisFreq,
                "VolFreq": VolFreq,
                "RotFreq": RotFreq,
                "RegrowthFreq": RegrowthFreq,
                "IntraSwapFreq": IntraSwapFreq,
                "CrankShaftFreq": CrankShaftFreq,
                "SwapFreq": SwapFreq,
                "OutputName": output_name_control_file_name_str,
                "EqSteps": EqSteps,
                "AdjSteps": AdjSteps,
                "PressureCalc": output_true_list_input,
                "RestartFreq": output_false_list_input,
                "CheckpointFreq": output_false_list_input,
                "ConsoleFreq": output_true_list_input,
                "BlockAverageFreq": output_false_list_input,
                "HistogramFreq": output_false_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_false_list_input,
                "CBMC_First": 12,
                "CBMC_Nth": 10,
                "CBMC_Ang": 50,
                "CBMC_Dih": 50,
            },
        )

        # **************************************************************
        # make the GOMC control file for each dihedral angle (END)
        # **************************************************************

        # *********************************
        # Write the restart .xsc file for GOMC (START)
        # *********************************
        gomc_restart_xsc_txt_file = open(
            f"{xyz_xsc_coor_files_directory}/starting_point.xsc", "w"
        )
        gomc_restart_xsc_txt_file.write(
            f"# GOMC extended system configuration output file\n"
        )
        gomc_restart_xsc_txt_file.write(
            f"#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w\n"
        )

        gomc_restart_xsc_txt_file.write(
            f"1 {Rcut.to_value('angstrom') * 2} 0 0 "
            f"0 {Rcut.to_value('angstrom') * 2} 0 "
            f"0 0 {Rcut.to_value('angstrom') * 2} "
            f"{Rcut.to_value('angstrom')} "
            f"{Rcut.to_value('angstrom')} "
            f"{Rcut.to_value('angstrom')} "
            f"0 0 0 "
            f"0 0 0"
            f"\n"
        )
        gomc_restart_xsc_txt_file.close()
        # *********************************
        # Write the restart .xsc file for GOMC (END)
        # *********************************

        # **************************************************************
        # Run GOMC and get the system/fragment energy (START)
        # **************************************************************

        # only run NVT as we only want/need the initial energy and want the box a constant size
        run_gomc_command = (
            f"cd {gomc_runs_folder_name} && {gomc_binary_path}/GOMC_CPU_NVT +p{gomc_cpu_cores} "
            f"{control_file_name_str} > {output_name_control_file_name_str}"
        )

        exec_gomc_run_command = subprocess.Popen(
            run_gomc_command, shell=True, stderr=subprocess.STDOUT
        )

        os.wait4(
            exec_gomc_run_command.pid, os.WSTOPPED
        )  # pauses python until exec_gomc_run_command sim done

        # **************************************************************
        # Run GOMC and get the system/fragment energy (END)
        # **************************************************************

        # open the gomc raw energy file so it can be writen over in the loop
        read_gomc_log_file = open(
            f"{gomc_runs_folder_name}/{output_name_control_file_name_str}", "r"
        ).readlines()
        get_e_titles = True
        for log_file_iter, log_file_line_iter in enumerate(read_gomc_log_file):
            log_file_splitline_iter = log_file_line_iter.split()

            # scan_iter starts at 1
            dihedral_angle_degrees = gaussian_raw_degrees_list[scan_iter - 1]

            # only open the gomc raw energy file and write header for 1st iteration (1)
            if len(log_file_splitline_iter) >= 2:
                if (
                    scan_iter == 1
                    and log_file_splitline_iter[0] == "ETITLE:"
                    and get_e_titles is True
                ):
                    gomc_combined_raw_energy = open(
                        gomc_raw_energy_filename, "w"
                    )
                    # remove the wrongly entered 5 spaces in before "ETITLE:
                    # (This will be fixed in GOMC so it is not required)
                    extra_spaces_for_header_space_gomc_bug = 5
                    gomc_combined_raw_energy.write(
                        f'{"Dihedral_Position": <19} '
                        f'{"Dihedral_Degrees": <19} '
                        f"{log_file_line_iter[extra_spaces_for_header_space_gomc_bug:]}"
                    )
                    if get_e_titles is True:
                        get_e_titles = False

                # get only the initial configuration energy lin (i.e., step 0 line)
                if (
                    log_file_splitline_iter[0] == "ENER_0:"
                    and log_file_splitline_iter[1] == "0"
                ):
                    gomc_combined_raw_energy.write(
                        f"{scan_iter: <19}"
                        f"{dihedral_angle_degrees: <19} "
                        f"{log_file_line_iter[5:]}"
                    )

    # This fails when the GOMC simulations did not run
    try:
        # close the gomc_combined_raw_energy file
        gomc_combined_raw_energy.close()

    except:
        raise ValueError(
            "ERROR: The GOMC simulations did not run. There is likely an error in creating the "
            "required GOMC files or user inputs to the desired files."
        )

    # **************************************************************
    # **************************************************************
    # Extract Energies from GOMC remove duplicate 0 point and
    # Change GOMC energies from K to kcal/mol
    # Change Gaussian energies from Hartree to kcal/mol
    # Write the restart .xsc file for GOMC.
    # Get GOMC and Gaussian data.
    # (START)
    # **************************************************************
    # **************************************************************

    conversion_hartree_to_kcal_per_mol = 627.509474063
    conversion_K_to_kcal_per_mol = (1 * u.Kelvin).to_value(
        "kcal/mol", equivalence="thermal"
    )

    # *********************************
    # get GOMC data (START)
    # *********************************
    # extract the raw data
    GOMC_data_df = pd.DataFrame(
        pd.read_csv(gomc_raw_energy_filename, sep="\s+")
    )
    GOMC_data_dihedral_degrees_list = GOMC_data_df.loc[
        :, "Dihedral_Degrees"
    ].tolist()
    GOMC_data_total_energy_K_list = GOMC_data_df.loc[:, "TOTAL"].tolist()

    # convert from Kelvin to kcal/mol normalize so the min value is 0
    GOMC_data_total_energy_kcal_per_mol_list = [
        i * conversion_K_to_kcal_per_mol for i in GOMC_data_total_energy_K_list
    ]
    GOMC_data_total_energy_kcal_per_mol_normalize_list = [
        i - min(GOMC_data_total_energy_kcal_per_mol_list)
        for i in GOMC_data_total_energy_kcal_per_mol_list
    ]

    # *********************************
    # get GOMC data (END)
    # *********************************

    # *********************************
    # get Gaussian data (START)
    # *********************************

    # extract the raw data
    gaussian_data_df = pd.DataFrame(
        pd.read_csv(qm_energy_file_dir_and_name, sep="\s+", header=3)
    )
    gaussian_data_dihedral_degrees_list = gaussian_data_df.iloc[:, 0].tolist()
    gaussian_data_total_energy_Hartree_list = gaussian_data_df.iloc[
        :, 1
    ].tolist()

    # convert from Hartree to kcal/mol energy units
    gaussian_data_total_energy_kcal_per_mol_list = [
        i * conversion_hartree_to_kcal_per_mol
        for i in gaussian_data_total_energy_Hartree_list
    ]

    # normalize so the min value is 0
    gaussian_data_total_energy_kcal_per_mol_normalize_list = [
        i - min(gaussian_data_total_energy_kcal_per_mol_list)
        for i in gaussian_data_total_energy_kcal_per_mol_list
    ]

    # get the Gaussian minus GOMC total energy and then it normalized
    Gaussian_minus_GOMC_data_dihedral_degrees_list = (
        GOMC_data_dihedral_degrees_list
    )
    Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list = [
        gaussian_data_total_energy_kcal_per_mol_normalize_list[i]
        - GOMC_data_total_energy_kcal_per_mol_normalize_list[i]
        for i in range(
            0, len(gaussian_data_total_energy_kcal_per_mol_normalize_list)
        )
    ]

    Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list = [
        Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list[i]
        - min(Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list)
        for i in range(
            0, len(Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list)
        )
    ]

    # *********************************
    # get Gaussian data (END)
    # *********************************

    # *********************************
    # get all other dihedral angles (phi) that match the atom type of the scanned dihedral,
    # and all sum 1/2*(k1 * scalar_i) = sum 1/2*(k1 * (1 +/- cos(n * phi)) values
    # (START)
    # *********************************
    # extract the data from the QM log file
    if qm_engine == "gaussian":
        [
            matching_dihedral_types_by_atom_numbers_list,
            matching_dihedral_types_by_atom_type_list,
            all_matching_dihedral_coordinates_angstroms_added_to_k_values_list,
            all_matching_dihedral_phi_degrees_added_to_k_values_list,
            all_sum_opls_const_1_plus_or_minus_cos_n_list,
        ] = mdf_frw.get_matching_dihedral_info_and_opls_fitting_data(
            fit_dihedral_atom_types,
            f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}.psf",
            qm_log_file_dict,
            mol2_file,
            qm_engine=qm_engine,
        )

    elif qm_engine == "gaussian_style_final_files":
        [
            matching_dihedral_types_by_atom_numbers_list,
            matching_dihedral_types_by_atom_type_list,
            all_matching_dihedral_coordinates_angstroms_added_to_k_values_list,
            all_matching_dihedral_phi_degrees_added_to_k_values_list,
            all_sum_opls_const_1_plus_or_minus_cos_n_list,
        ] = mdf_frw.get_matching_dihedral_info_and_opls_fitting_data(
            fit_dihedral_atom_types,
            f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}.psf",
            qm_log_file_dict,
            mol2_file,
            qm_engine=qm_engine,
            manual_dihedral_atom_numbers_list=manual_dihedral_atom_numbers_list,
        )
    else:
        raise ValueError(
            f"ERROR: The entered qm_engine = {qm_engine} and the only valid options are {['gaussian']}"
        )

    # get the individual sum_opls_const_1_plus_or_minus_cos_n_list ones for fitting and plotting
    const_1_minus_Cos_0_phi_data_lists = []
    const_1_plus_Cos_1_phi_data_lists = []
    const_1_minus_Cos_2_phi_data_lists = []
    const_1_plus_Cos_3_phi_data_lists = []
    const_1_minus_Cos_4_phi_data_lists = []
    for (
        const_1_plus_or_minus_cos_i
    ) in all_sum_opls_const_1_plus_or_minus_cos_n_list:
        const_1_minus_Cos_0_phi_data_lists.append(
            const_1_plus_or_minus_cos_i[0]
        )
        const_1_plus_Cos_1_phi_data_lists.append(const_1_plus_or_minus_cos_i[1])
        const_1_minus_Cos_2_phi_data_lists.append(
            const_1_plus_or_minus_cos_i[2]
        )
        const_1_plus_Cos_3_phi_data_lists.append(const_1_plus_or_minus_cos_i[3])
        const_1_minus_Cos_4_phi_data_lists.append(
            const_1_plus_or_minus_cos_i[4]
        )

    # *********************************
    # get all other dihedral angles (phi) that match the atom type of the scanned dihedral,
    # and all sum 1/2*(k1 * scalar_i) = sum 1/2*(k1 * (1 +/- cos(n * phi)) values
    # (END)
    # *********************************

    # Check if all the columns are the same length for GOMC and Gaussian data
    if (
        not len(GOMC_data_dihedral_degrees_list)
        == len(GOMC_data_total_energy_kcal_per_mol_list)
        == len(GOMC_data_total_energy_kcal_per_mol_normalize_list)
        == len(gaussian_data_dihedral_degrees_list)
        == len(gaussian_data_total_energy_kcal_per_mol_list)
        == len(gaussian_data_total_energy_kcal_per_mol_normalize_list)
        == len(Gaussian_minus_GOMC_data_dihedral_degrees_list)
        == len(Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list)
        == len(
            Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
        )
        == len(all_sum_opls_const_1_plus_or_minus_cos_n_list)
    ):
        raise ValueError(
            "ERROR: The GOMC and gaussian outputs do not match in length. "
            "This could mean something is changed and wrong in the code, "
            "or GOMC is outputting multiple Initial eneries in the log file "
            ", in this case use a new version of GOMC."
        )

    # Check if all the angles match between sorted GOMC and Gaussian data
    for _ in range(0, len(GOMC_data_dihedral_degrees_list)):
        if not len(GOMC_data_dihedral_degrees_list) == len(
            gaussian_data_dihedral_degrees_list
        ):
            raise ValueError(
                "ERROR: The GOMC and gaussian output angles are not in the same angles in order."
            )

    # Check if all the angles match between sorted GOMC and Gaussian data
    for k_angle in range(0, len(GOMC_data_dihedral_degrees_list)):
        if not len(GOMC_data_dihedral_degrees_list) == len(
            gaussian_data_dihedral_degrees_list
        ):
            raise ValueError(
                "ERROR: The GOMC and gaussian output angles are not in the same angles in order."
            )
        if k_angle == 0:
            # write out the GOMC and Gaussian data in a file
            gomc_gaussian_kcal_per_mol_energy_data_txt_file = open(
                gomc_gaussian_kcal_per_mol_energy_filename, "w"
            )
            gomc_gaussian_kcal_per_mol_energy_data_txt_file.write(
                f"{'Dihedral_Degrees': <30} "
                f"{'GOMC_E_kcal_per_mol': <30} "
                f"{'Gaussian_E_kcal_per_mol': <30} "
                f"{'Gaussian_minus_GOMC_E_kcal_per_mol': <40} "
                f"{'const_1_minus_Cos_0_phi': <30} "
                f"{'const_1_plus_Cos_1_phi': <30} "
                f"{'const_1_minus_Cos_2_phi': <30} "
                f"{'const_1_plus_Cos_3_phi': <30} "
                f"{'const_1_minus_Cos_4_phi': <30} "
                f" \n"
            )

        gomc_gaussian_kcal_per_mol_energy_data_txt_file.write(
            f"{Gaussian_minus_GOMC_data_dihedral_degrees_list[k_angle]: <30} "
            f"{GOMC_data_total_energy_kcal_per_mol_normalize_list[k_angle]: <30} "
            f"{gaussian_data_total_energy_kcal_per_mol_normalize_list[k_angle]: <30} "
            f"{Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list[k_angle]: <40} "
            f"{const_1_minus_Cos_0_phi_data_lists[k_angle]: <30} "
            f"{const_1_plus_Cos_1_phi_data_lists[k_angle]: <30} "
            f"{const_1_minus_Cos_2_phi_data_lists[k_angle]: <30} "
            f"{const_1_plus_Cos_3_phi_data_lists[k_angle]: <30} "
            f"{const_1_minus_Cos_4_phi_data_lists[k_angle]: <30} "
            f" \n"
        )
    gomc_gaussian_kcal_per_mol_energy_data_txt_file.close()

    # *********************************
    # get GOMC and Gaussian data (END)
    # *********************************

    # **************************************************************
    # **************************************************************
    # Extract Energies from GOMC remove duplicate 0 point and
    # Change GOMC energies from K to kcal/mol
    # Change Gaussian energies from Hartree to kcal/mol
    # Write the restart .xsc file for GOMC.
    # Get GOMC and Gaussian data.
    # (END)
    # **************************************************************
    # **************************************************************

    # *********************************
    # fit the Gaussian - GOMC dihedral (START)
    # *********************************

    fig1, (ax1) = plt.subplots(1)
    axis_Label_font_size = 12
    legend_font_size = 12
    ax1.set_xlabel("phi (degrees)", size=axis_Label_font_size)
    ax1.set_ylabel("Dihedral Energy (kcal/mol)", size=axis_Label_font_size)

    # sort the by Dihedral_Degrees at the same time for
    # Dihedral_Degrees, GOMC_E_kcal_per_mol, Gaussian_E_kcal_per_mol, and Gaussian_minus_GOMC_E_kcal_per_mol
    (
        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list,
        sorted_GOMC_data_total_energy_kcal_per_mol_normalize_list,
        sorted_gaussian_data_total_energy_kcal_per_mol_normalize_list,
        sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list,
        sorted_all_sum_opls_const_1_plus_or_minus_cos_n_list,
        sorted_const_1_minus_Cos_0_phi_data_lists,
        sorted_const_1_plus_Cos_1_phi_data_lists,
        sorted_const_1_minus_Cos_2_phi_data_lists,
        sorted_const_1_plus_Cos_3_phi_data_lists,
        sorted_const_1_minus_Cos_4_phi_data_lists,
    ) = zip(
        *sorted(
            zip(
                Gaussian_minus_GOMC_data_dihedral_degrees_list,
                GOMC_data_total_energy_kcal_per_mol_normalize_list,
                gaussian_data_total_energy_kcal_per_mol_normalize_list,
                Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list,
                all_sum_opls_const_1_plus_or_minus_cos_n_list,
                const_1_minus_Cos_0_phi_data_lists,
                const_1_plus_Cos_1_phi_data_lists,
                const_1_minus_Cos_2_phi_data_lists,
                const_1_plus_Cos_3_phi_data_lists,
                const_1_minus_Cos_4_phi_data_lists,
            )
        )
    )

    plot_max = int(
        max(
            sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
        )
        + 1.51
    )
    ax1.set_ylim(-int(plot_max * 1.5), plot_max)
    plt.title(
        "OPLS: All viable summed dihedral fits $\Longrightarrow$ $\Sigma$ matching dihedrals. \n"
        "(Non-zero k's = 1 and 3 -> label k_non_0='1_3')"
    )

    # loop thru dihderal_k_zeros_list_k0_k1_k2_k3_k4 list and fit all that are listed in here
    # add add the label designnator k's used (i.e., Non-zero k's = 1 and 3 -> label '1_3')
    # and write the constants out
    end_part_dihedral_k_constants_fit_energy_figure_filename = (
        "all_summed_dihedrals_k_constants_figure.pdf"
    )

    all_individual_fit_dihedral_k_constants_figure_filename = (
        "all_single_fit_dihedral_k_constants_figure.pdf"
    )

    end_part_dihedral_k_constants_fit_energy_filename = (
        "k_constants_fit_energy.txt"
    )
    opls_dihedral_k_constants_fit_energy_kcal_per_mol_txt_file = open(
        f"{'opls_dihedral'}_{end_part_dihedral_k_constants_fit_energy_filename}",
        "w",
    )
    opls_dihedral_k_constants_fit_energy_kcal_per_mol_txt_file.write(
        f"{'non_zero_k_constants': <25} "
        f"{'k0_kcal_per_mol': <25} "
        f"{'k1_kcal_per_mol': <25} "
        f"{'k2_kcal_per_mol': <25} "
        f"{'k3_kcal_per_mol': <25} "
        f"{'k4_kcal_per_mol': <25} "
        f"{'r_squared': <25} "
    )

    # **********************************
    # Select only the OPLS fit cos powers that produce a unique solution
    # without any symmetry issues
    # (START)
    # **********************************

    # set to kn or n variable fit for that you want as a string in a list with the variable name 'fit_k_list'
    # (i.e, n = 1 and 3 --> "1_3", i.e, n = 1 2, and 3 --> "1_2_3"))
    #
    fit_k_list = [
        "1",
        "2",
        "3",
        "4",
        "1_3",
        "2_4",
        "3_4",
        "1_2",
        "1_2_3",
        "1_2_3_4",
    ]

    # Determine which power values you can use.
    # If any of the const_1_plus_Cos_1_phi_data, const_1_minus_Cos_2_phi_data,
    # const_1_plus_Cos_3_phi_data, or const_1_minus_Cos_4_phi_data values
    # are constant for the whole rotation, they can not be used in the dihedral fit.
    # Check if these sums are constant or variable here, and make a list of which
    # cos powers can be used.
    allowed_fitting_powers_1_2_3_and_4 = []
    check_sig_figs = 2
    const_cos_count_powers_1_2_3_and_4_int_list = [0, 0, 0, 0]

    # set the number of differing cos constant terms needed to determine if it
    # can be used in the fitting. If there is only 1 that is different in a large
    # group, it may not be wrong, if it is not minimized in QM....
    number_of_differnt_const_cos_needed_int = int(
        25 / 100 * len(sorted_const_1_plus_Cos_1_phi_data_lists) + 1
    )
    for ck_const_i in range(0, len(sorted_const_1_plus_Cos_1_phi_data_lists)):
        # Check if can use cos power 1
        if (
            bool(
                np.isclose(
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_plus_Cos_1_phi_data_lists[ck_const_i],
                        sig_figs=check_sig_figs,
                    ),
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_plus_Cos_1_phi_data_lists[0],
                        sig_figs=check_sig_figs,
                    ),
                )
            )
            is False
        ):
            const_cos_count_powers_1_2_3_and_4_int_list[0] = (
                1 + const_cos_count_powers_1_2_3_and_4_int_list[0]
            )

        # Check if can use cos power 2
        if (
            bool(
                np.isclose(
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_minus_Cos_2_phi_data_lists[ck_const_i],
                        sig_figs=check_sig_figs,
                    ),
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_minus_Cos_2_phi_data_lists[0],
                        sig_figs=check_sig_figs,
                    ),
                )
            )
            is False
        ):
            const_cos_count_powers_1_2_3_and_4_int_list[1] = (
                1 + const_cos_count_powers_1_2_3_and_4_int_list[1]
            )

        # Check if can use cos power 3
        if (
            bool(
                np.isclose(
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_plus_Cos_3_phi_data_lists[ck_const_i],
                        sig_figs=check_sig_figs,
                    ),
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_plus_Cos_3_phi_data_lists[0],
                        sig_figs=check_sig_figs,
                    ),
                )
            )
            is False
        ):
            const_cos_count_powers_1_2_3_and_4_int_list[2] = (
                1 + const_cos_count_powers_1_2_3_and_4_int_list[2]
            )

        # Check if can use cos power 4
        if (
            bool(
                np.isclose(
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_minus_Cos_4_phi_data_lists[ck_const_i],
                        sig_figs=check_sig_figs,
                    ),
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_minus_Cos_4_phi_data_lists[0],
                        sig_figs=check_sig_figs,
                    ),
                )
            )
            is False
        ):
            const_cos_count_powers_1_2_3_and_4_int_list[3] = (
                1 + const_cos_count_powers_1_2_3_and_4_int_list[3]
            )

        # add the cos power value if enough differences in the constants are seen
        if (
            const_cos_count_powers_1_2_3_and_4_int_list[0]
            >= number_of_differnt_const_cos_needed_int
        ):
            allowed_fitting_powers_1_2_3_and_4.append("1")

        if (
            const_cos_count_powers_1_2_3_and_4_int_list[1]
            >= number_of_differnt_const_cos_needed_int
        ):
            allowed_fitting_powers_1_2_3_and_4.append("2")

        if (
            const_cos_count_powers_1_2_3_and_4_int_list[2]
            >= number_of_differnt_const_cos_needed_int
        ):
            allowed_fitting_powers_1_2_3_and_4.append("3")

        if (
            const_cos_count_powers_1_2_3_and_4_int_list[3]
            >= number_of_differnt_const_cos_needed_int
        ):
            allowed_fitting_powers_1_2_3_and_4.append("4")

    allowed_fitting_powers_1_2_3_and_4.sort()

    # This selects only fit acceptable power based on dihedral symmetry of the molecule/fragment.
    # Not all the powers can be used when fitting a dihedral, so this outputs only the
    # powers that will produce the correct solution (i.e., removes all powers that don't
    # yeild the correct solution).
    fit_k_list_allowed = []
    for mod_i, mod_val in enumerate(fit_k_list):
        power_allowed_iter = True
        for str_i in mod_val:
            if (
                str_i not in ["-", "_"]
                and str_i not in allowed_fitting_powers_1_2_3_and_4
            ):
                power_allowed_iter = False
        if power_allowed_iter is True:
            fit_k_list_allowed.append(mod_val)

    # **********************************
    # Select only the OPLS fit cos powers that produce a unique solution
    # without any symmetry issues
    # (END)
    # **********************************

    # modify 'sorted_const_1_minus_Cos_0_phi_data_lists' based on "opls_force_k0_zero"
    # which allows k0=0 or a k0=constant.
    # If 'sorted_const_1_minus_Cos_0_phi_data_lists' = all 0s, then k0=0,
    # which is changed to 'k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists'.
    # If 'sorted_const_1_minus_Cos_0_phi_data_lists' = all 1s, then k0=constant,
    # which is changed to 'k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists'.
    # It is critical that 'const_1_minus_Cos_0_phi_data' be all (0, 0, .., 0) for k0=0
    # and all (1, 1, .., 1) 'const_1_minus_Cos_0_phi_data' for k0=constant.
    k0_forced_to_0_list = []
    k0_is_constant_list = []
    for determine_k0_m in range(
        0, len(sorted_const_1_minus_Cos_0_phi_data_lists)
    ):
        k0_forced_to_0_list.append(0)
        k0_is_constant_list.append(1)

    k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists = tuple(
        k0_forced_to_0_list
    )
    k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists = tuple(
        k0_is_constant_list
    )

    # Run the fitting for only the allowed power types
    for k_type_i in fit_k_list_allowed:
        # make the list of k_type_i for fitting in the data
        k_type_list_i = []
        for _ in range(
            0, len(sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list)
        ):
            k_type_list_i.append(k_type_i)

        if k_type_i == "1":
            # ***** (k=0  -> forced) ****
            # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )

            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[2] = 0
            parameters_k0_forced_to_0[3] = 0
            parameters_k0_forced_to_0[4] = 0

            # fit the OPLS dihedral  (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )

            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.  (k=constant)
            parameters_k0_is_constant[2] = 0
            parameters_k0_is_constant[3] = 0
            parameters_k0_is_constant[4] = 0

            # fit the OPLS dihedral  (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "2":
            # ***** (k=0  -> forced) ****
            # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.  (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[1] = 0
            parameters_k0_forced_to_0[3] = 0
            parameters_k0_forced_to_0[4] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )

            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)
            parameters_k0_is_constant[1] = 0
            parameters_k0_is_constant[3] = 0
            parameters_k0_is_constant[4] = 0

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "3":
            # ***** (k=0  -> forced) ****
            # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[1] = 0
            parameters_k0_forced_to_0[2] = 0
            parameters_k0_forced_to_0[4] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)
            parameters_k0_is_constant[1] = 0
            parameters_k0_is_constant[2] = 0
            parameters_k0_is_constant[4] = 0

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "4":
            # ***** (k=0  -> forced) ****
            # # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[1] = 0
            parameters_k0_forced_to_0[2] = 0
            parameters_k0_forced_to_0[3] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)
            parameters_k0_is_constant[1] = 0
            parameters_k0_is_constant[2] = 0
            parameters_k0_is_constant[3] = 0

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "1_3":
            # ***** (k=0  -> forced) ****
            # # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[2] = 0
            parameters_k0_forced_to_0[4] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)
            parameters_k0_is_constant[2] = 0
            parameters_k0_is_constant[4] = 0

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "2_4":
            # ***** (k=0  -> forced) ****
            # # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[1] = 0
            parameters_k0_forced_to_0[3] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)
            parameters_k0_is_constant[1] = 0
            parameters_k0_is_constant[3] = 0

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "1_2":
            # ***** (k=0  -> forced) ****
            # # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[3] = 0
            parameters_k0_forced_to_0[4] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)
            parameters_k0_is_constant[3] = 0
            parameters_k0_is_constant[4] = 0

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "3_4":
            # ***** (k=0  -> forced) ****
            # # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[1] = 0
            parameters_k0_forced_to_0[2] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)
            parameters_k0_is_constant[1] = 0
            parameters_k0_is_constant[2] = 0

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "1_2_3":
            # ***** (k=0  -> forced) ****
            # # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)

            parameters_k0_forced_to_0[0] = 0
            parameters_k0_forced_to_0[4] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)
            parameters_k0_is_constant[4] = 0

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        elif k_type_i == "1_2_3_4":
            # ***** (k=0  -> forced) ****
            # # get parameeters and covariance the OPLS dihedral (k=0  -> forced)
            parameters_k0_forced_to_0, covariance_k0_forced_to_0 = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=0  -> forced)
            parameters_k0_forced_to_0[0] = 0

            # fit the OPLS dihedral (k=0  -> forced)
            fit_opls_dihedral_k0_forced_to_0 = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_forced_to_0_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            )

            # ***** (k=constant) ****
            # # get parameeters and covariance the OPLS dihedral (k=constant)
            parameters_k0_is_constant, covariance_k0_is_constant = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list
                ),
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1. (k=constant)

            # fit the OPLS dihedral (k=constant)
            fit_opls_dihedral_k0_is_constant = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list
                    ),
                    np.asarray(
                        k0_is_constant_sorted_const_1_minus_Cos_0_phi_data_lists
                    ),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists),
                ),
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            )

        else:
            raise ValueError(
                f"ERROR: The {k_type_i} selected in the 'fit_k_list' variable is not a valid selection"
            )

        # shift dihedral so always at energy = 0 minimum (changing 'fit_opls_dihedral' + C and 'parameters[0]' + C,
        # and 'sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list' + C)
        theta_spacing_degree = 1

        # Can not do the min shift for k0_forced_to_0 or you add k0
        fitted_opls_dihedral_small_theta_spacing_k0_forced_to_0 = []
        for theta_i in range(-180, 180, 1):
            fitted_opls_dihedral_small_theta_spacing_k0_forced_to_0.append(
                mdf_math.opls_dihedral_n_1_2_3_4(
                    theta_i,
                    parameters_k0_forced_to_0[0],
                    parameters_k0_forced_to_0[1],
                    parameters_k0_forced_to_0[2],
                    parameters_k0_forced_to_0[3],
                    parameters_k0_forced_to_0[4],
                )
            )

        # min shift to 0 for k0_is_constant
        fitted_opls_dihedral_small_theta_spacing_k0_is_constant = []
        for theta_i in range(-180, 180, 1):
            fitted_opls_dihedral_small_theta_spacing_k0_is_constant.append(
                mdf_math.opls_dihedral_n_1_2_3_4(
                    theta_i,
                    parameters_k0_is_constant[0],
                    parameters_k0_is_constant[1],
                    parameters_k0_is_constant[2],
                    parameters_k0_is_constant[3],
                    parameters_k0_is_constant[4],
                )
            )

        min_0_fitted_opls_dihedral_small_theta_spacing_k0_is_constant = np.min(
            fitted_opls_dihedral_small_theta_spacing_k0_is_constant
        )

        # subtract minimum to  k0_is_constant[0] (shift dihedral so always at energy = 0 minimum )
        # remember min times 2 as OPLS_energy = 1/2 (k0+...)
        if opls_force_k0_zero is False:
            parameters_k0_is_constant[0] -= (
                min_0_fitted_opls_dihedral_small_theta_spacing_k0_is_constant
                * 2
            )

        # subtract minimum to  parameters[0] (shift dihedral so always at energy = 0 minimum )
        for dih_energy_m in range(0, len(fit_opls_dihedral_k0_is_constant)):
            fit_opls_dihedral_k0_is_constant[
                dih_energy_m
            ] -= min_0_fitted_opls_dihedral_small_theta_spacing_k0_is_constant

        # calulate TSS, RSS and R**2 (k0_is_constant)
        r_squared_k0_is_constant = mdf_math.get_r_squared(
            sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list,
            fit_opls_dihedral_k0_is_constant,
        )

        # calulate TSS, RSS and R**2 (k0_forced_to_0)
        r_squared_k0_forced_to_0 = mdf_math.get_r_squared(
            sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list,
            fit_opls_dihedral_k0_forced_to_0,
        )

        # Determine what fit to take 'k0_forced_to_0' or 'k0_is_constant' and
        # select the best fit from them
        if (
            opls_force_k0_zero is False
            and r_squared_k0_is_constant >= r_squared_k0_forced_to_0
        ):
            r_squared = r_squared_k0_is_constant
            fit_opls_dihedral = fit_opls_dihedral_k0_is_constant
            parameters = [
                parameters_k0_is_constant[0],
                parameters_k0_is_constant[1],
                parameters_k0_is_constant[2],
                parameters_k0_is_constant[3],
                parameters_k0_is_constant[4],
            ]

        else:
            r_squared = r_squared_k0_forced_to_0
            fit_opls_dihedral = fit_opls_dihedral_k0_forced_to_0
            parameters = [
                parameters_k0_forced_to_0[0],
                parameters_k0_forced_to_0[1],
                parameters_k0_forced_to_0[2],
                parameters_k0_forced_to_0[3],
                parameters_k0_forced_to_0[4],
            ]

        plt.plot(
            np.asarray(sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
            np.asarray(fit_opls_dihedral),
            "-",
            label=f" {k_type_i} | $R^{2}$={np.round(r_squared, decimals=4)}",
        )

        # wrie out the k constants and R^2
        opls_dihedral_k_constants_fit_energy_kcal_per_mol_txt_file.write(
            f"\n{k_type_i: <25} "
            f"{mdf_math.round_to_sig_figs(parameters[0], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(parameters[1], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(parameters[2], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(parameters[3], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(parameters[4], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(r_squared, sig_figs=12): <25} "
        )

    # plot the data point that it is being fit too
    plt.plot(
        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list,
        sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list,
        "o",
        label="Data",
    )

    # close the file
    opls_dihedral_k_constants_fit_energy_kcal_per_mol_txt_file.close()

    major_xticks = np.arange(-180, 180 + 0.0001, 60)
    minor_xticks = np.arange(-180, 180 + 0.0001, 10)
    ax1.set_xticks(major_xticks)
    ax1.set_xticks(minor_xticks, minor=True)
    ax1.set_xlim(-180, 180 + 0.001)

    plt.legend(
        ncol=2,
        loc="lower center",
        fontsize=legend_font_size,
        prop={"size": legend_font_size},
    )

    # plt.show()
    fig1.savefig(
        f"opls_{end_part_dihedral_k_constants_fit_energy_figure_filename}",
        dpi=300,
    )

    # *********************************
    # fit the Gaussian - GOMC dihedral (END)
    # *********************************

    # *********************************
    # Plot all OPLS dihedral to fitted forms together (START)
    # *********************************
    opls_fit_data_df = pd.DataFrame(
        pd.read_csv(
            f"{'opls_dihedral'}_{end_part_dihedral_k_constants_fit_energy_filename}",
            sep="\s+",
            header=0,
        )
    )
    opls_fit_data_non_zero_k_constants_list = opls_fit_data_df.loc[
        :, "non_zero_k_constants"
    ].tolist()
    opls_fit_data_k0_kcal_per_mol_list = opls_fit_data_df.loc[
        :, "k0_kcal_per_mol"
    ].tolist()
    opls_fit_data_k1_kcal_per_mol_list = opls_fit_data_df.loc[
        :, "k1_kcal_per_mol"
    ].tolist()
    opls_fit_data_k2_kcal_per_mol_list = opls_fit_data_df.loc[
        :, "k2_kcal_per_mol"
    ].tolist()
    opls_fit_data_k3_kcal_per_mol_list = opls_fit_data_df.loc[
        :, "k3_kcal_per_mol"
    ].tolist()
    opls_fit_data_k4_kcal_per_mol_list = opls_fit_data_df.loc[
        :, "k4_kcal_per_mol"
    ].tolist()
    opls_fit_data_r_squared_list = opls_fit_data_df.loc[:, "r_squared"].tolist()

    # *********************************
    # Plot all OPLS dihedral to fitted forms together  (END)
    # *********************************

    # *********************************
    # Convert the OPLS dihedral to other forms (START)
    # *********************************

    # create the Periodic / CHARMM dihedrals file
    periodic_dihedral_k_constants_fit_energy_kcal_per_mol_txt_file = open(
        f"{'periodic_dihedral'}_{end_part_dihedral_k_constants_fit_energy_filename}",
        "w",
    )
    periodic_dihedral_k_constants_fit_energy_kcal_per_mol_txt_file.write(
        f"{'non_zero_k_constants': <25} "
        f"{'k0_kcal_per_mol': <25} "
        f"{'k1_kcal_per_mol': <25} "
        f"{'k2_kcal_per_mol': <25} "
        f"{'k3_kcal_per_mol': <25} "
        f"{'k4_kcal_per_mol': <25} "
        f"{'k5_kcal_per_mol': <25} "
        f"{'n0_kcal_per_mol': <25} "
        f"{'n1_kcal_per_mol': <25} "
        f"{'n2_kcal_per_mol': <25} "
        f"{'n3_kcal_per_mol': <25} "
        f"{'n4_kcal_per_mol': <25} "
        f"{'n5_kcal_per_mol': <25} "
        f"{'d0_kcal_per_mol': <25} "
        f"{'d1_kcal_per_mol': <25} "
        f"{'d2_kcal_per_mol': <25} "
        f"{'d3_kcal_per_mol': <25} "
        f"{'d4_kcal_per_mol': <25} "
        f"{'d5_kcal_per_mol': <25} "
        f"{'r_squared': <25} "
    )

    # create the RB torsions file
    RB_torsion_k_constants_fit_energy_kcal_per_mol_txt_file = open(
        f"{'RB_torsion'}_{end_part_dihedral_k_constants_fit_energy_filename}",
        "w",
    )
    RB_torsion_k_constants_fit_energy_kcal_per_mol_txt_file.write(
        f"{'non_zero_k_constants': <25} "
        f"{'k0_kcal_per_mol': <25} "
        f"{'k1_kcal_per_mol': <25} "
        f"{'k2_kcal_per_mol': <25} "
        f"{'k3_kcal_per_mol': <25} "
        f"{'k4_kcal_per_mol': <25} "
        f"{'k5_kcal_per_mol': <25} "
        f"{'r_squared': <25} "
    )

    # loop thru the different 'non_zero_k_constants' for the OPLS dihedral
    for opls_fit_i in range(0, len(opls_fit_data_non_zero_k_constants_list)):
        # list the phi values to check
        phi_check_degrees = 1
        phi_check_number_of_degree_values = int(360 / phi_check_degrees)
        phi_values_for_check_degrees_list = [
            -180 + i * phi_check_degrees
            for i in range(0, phi_check_number_of_degree_values)
        ]

        cos_power_list = [
            opls_fit_data_non_zero_k_constants_list[opls_fit_i]
            for i in range(0, len(phi_values_for_check_degrees_list))
        ]

        # check if the periodic and opls dihedral energies match all the values in the list
        for period_opls_ck_i in range(
            0, len(phi_values_for_check_degrees_list)
        ):
            # calculate the dihedral energy for a check to the converted dihedrals
            opls_dihedral_energy_check = mdf_math.opls_dihedral(
                (
                    cos_power_list[opls_fit_i],
                    phi_values_for_check_degrees_list[opls_fit_i],
                    None,
                    None,
                    None,
                    None,
                    None,
                ),
                opls_fit_data_k0_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k1_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k2_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k3_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k4_kcal_per_mol_list[opls_fit_i],
            )

            # check if the periodic and opls dihedral energies match all the values in the list
            # for period_opls_ck_i in range(0, len(phi_values_for_check_degrees_list)):

            # *********************************
            # Periodic / CHARMM dihedrals form calculations and checks (transformed from OPLS) (START)
            # *********************************
            # Convert the OPLS fit to Periodic/CHARMM style dihedral (function = OPLS_to_periodic(f0, f1, f2, f3, f4))
            periodic_dihedral_k_n_d_values = OPLS_to_periodic(
                opls_fit_data_k0_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k1_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k2_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k3_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k4_kcal_per_mol_list[opls_fit_i],
            )

            periodic_dihedral_energy_check = (
                mdf_math.periodic_dihedral_n_1_2_3_4_5(
                    phi_values_for_check_degrees_list[opls_fit_i],
                    periodic_dihedral_k_n_d_values[0][0],
                    periodic_dihedral_k_n_d_values[1][0],
                    periodic_dihedral_k_n_d_values[2][0],
                    periodic_dihedral_k_n_d_values[3][0],
                    periodic_dihedral_k_n_d_values[4][0],
                    periodic_dihedral_k_n_d_values[5][0],
                    periodic_dihedral_k_n_d_values[0][1],
                    periodic_dihedral_k_n_d_values[1][1],
                    periodic_dihedral_k_n_d_values[2][1],
                    periodic_dihedral_k_n_d_values[3][1],
                    periodic_dihedral_k_n_d_values[4][1],
                    periodic_dihedral_k_n_d_values[5][1],
                    periodic_dihedral_k_n_d_values[0][2],
                    periodic_dihedral_k_n_d_values[1][2],
                    periodic_dihedral_k_n_d_values[2][2],
                    periodic_dihedral_k_n_d_values[3][2],
                    periodic_dihedral_k_n_d_values[4][2],
                    periodic_dihedral_k_n_d_values[5][2],
                )
            )

            # check periodic dihedral equation fit matches the OPLS value
            if not math.isclose(
                opls_dihedral_energy_check, periodic_dihedral_energy_check
            ):
                raise ValueError(
                    "ERROR: The OPLS and periodic/CHARMM style dihedral energies do not match. "
                    "The is likely a conversion error in the the constants between OPLS and the periodic/CHARMM style."
                )

            # *********************************
            # Periodic / CHARMM dihedrals form calculations and checks (transformed from OPLS)  (END)
            # *********************************

            # *********************************
            # RB torsions form calculations and checks (transformed from OPLS) (START)
            # NOTE: For RB torsions, the angle psi is used, which is psi = phi - Pi
            # *********************************
            # function = OPLS_to_RB(f0, f1, f2, f3, f4, error_tolerance=1e-4):

            # mdf_math.RB_torsion_n_1_2_3_4_5(phi_data, c_0, c_1, c_2, c_3, c_4, c_5):

            # Convert the OPLS fit to Periodic/CHARMM style dihedral
            # (function = OPLS_to_RB(f0, f1, f2, f3, f4, error_tolerance=1e-4)).
            RB_torsion_k_values = OPLS_to_RB(
                opls_fit_data_k0_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k1_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k2_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k3_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k4_kcal_per_mol_list[opls_fit_i],
            )

            RB_torsion_energy_check = mdf_math.RB_torsion_n_1_2_3_4_5(
                phi_values_for_check_degrees_list[opls_fit_i],
                RB_torsion_k_values[0],
                RB_torsion_k_values[1],
                RB_torsion_k_values[2],
                RB_torsion_k_values[3],
                RB_torsion_k_values[4],
                RB_torsion_k_values[5],
            )

            # check periodic dihedral equation fit matches the OPLS value
            if not math.isclose(
                opls_dihedral_energy_check, RB_torsion_energy_check
            ):
                raise ValueError(
                    "ERROR: The OPLS dihedral and RB torsion style energies do not match. "
                    "The is likely a conversion error in the the constants between OPLS dihedral and RB torsion style."
                )

        # *********************************
        # RB torsions form calculations and checks (transformed from OPLS)  (END)
        # NOTE: For RB torsions, the angle psi is used, which is psi = phi - Pi
        # *********************************

        # *********************************
        # Write out the periodic and RB torsions values to file (transformed from OPLS)  (START)
        # NOTE: For RB torsions, the angle psi is used, which is psi = phi - Pi
        # *********************************

        # Write the periodic dihedral constants to file.
        # Reused the OPLS R^2 fitted values as the energeries where
        # validated between energies periodic/CHARMM dihedral style and OPLS.
        # Added the '0_' to the 'non_zero_k_constants' as the k0 needed
        # for the OPLS conversion to the periodic/CHARMM style.
        periodic_dihedral_k_constants_fit_energy_kcal_per_mol_txt_file.write(
            f"\n{f'0_{opls_fit_data_non_zero_k_constants_list[opls_fit_i]}': <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[0][0], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[1][0], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[2][0], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[3][0], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[4][0], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[5][0], sig_figs=12): <25} "
            f"{int(periodic_dihedral_k_n_d_values[0][1]): <25} "
            f"{int(periodic_dihedral_k_n_d_values[1][1]): <25} "
            f"{int(periodic_dihedral_k_n_d_values[2][1]): <25} "
            f"{int(periodic_dihedral_k_n_d_values[3][1]): <25} "
            f"{int(periodic_dihedral_k_n_d_values[4][1]): <25} "
            f"{int(periodic_dihedral_k_n_d_values[5][1]): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[0][2], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[1][2], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[2][2], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[3][2], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[4][2], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(periodic_dihedral_k_n_d_values[5][2], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(opls_fit_data_r_squared_list[opls_fit_i], sig_figs=12): <25} "
        )

        # Write the RB torsion constants to file.
        # Reused the OPLS R^2 fitted values as the energeries where
        # validated between energies RB torsion style and OPLS.
        # Added the '0_' to the 'non_zero_k_constants' as the k0 needed
        # for the OPLS conversion to the RB torsion style.
        RB_torsion_k_constants_fit_energy_kcal_per_mol_txt_file.write(
            f"\n{f'0_{opls_fit_data_non_zero_k_constants_list[opls_fit_i]}': <25} "
            f"{mdf_math.round_to_sig_figs(RB_torsion_k_values[0], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(RB_torsion_k_values[1], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(RB_torsion_k_values[2], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(RB_torsion_k_values[3], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(RB_torsion_k_values[4], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(RB_torsion_k_values[5], sig_figs=12): <25} "
            f"{mdf_math.round_to_sig_figs(opls_fit_data_r_squared_list[opls_fit_i], sig_figs=12): <25} "
        )

        # *********************************
        # Write out the periodic and RB torsions values to file (transformed from OPLS)  (END)
        # NOTE: For RB torsions, the angle psi is used, which is psi = phi - Pi
        # *********************************

    # close the Periodic / CHARMM dihedrals file
    periodic_dihedral_k_constants_fit_energy_kcal_per_mol_txt_file.close()

    # close the RB torsions file
    RB_torsion_k_constants_fit_energy_kcal_per_mol_txt_file.close()

    # *********************************
    # Convert the OPLS dihedral to other forms (END)
    # *********************************

    # *********************************
    # Check the all the OPLS dihedral forms are correct
    # by running GOMC with the fitted values and comparing it to QM
    # (START)
    # *********************************
    opls_r_squared_fitted_data_via_gomc_list = []
    for opls_q, opls_fit_q in enumerate(
        opls_fit_data_non_zero_k_constants_list
    ):
        gomc_fitted_gaussian_kcal_per_mol_energy_filename = (
            f"all_normalized_energies_OPLS_fit_{opls_fit_q}_in_kcal_per_mol.txt"
        )
        # write the GOMC control files
        for scan_iter_q in range(1, len(gaussian_raw_degrees_list) + 1):
            read_gomc_fitted_restart_file_coor_dir_and_name = f"../{xyz_xsc_coor_files_directory}/dihedral_coords_position_{scan_iter_q}.coor"
            read_gomc_fitted_restart_file_xsc_dir_and_name = (
                f"../{xyz_xsc_coor_files_directory}/starting_point.xsc"
            )

            # The gomc raw energy filename in Kelvin Energies
            gomc_raw_energy_fitted_filename = (
                f"gomc_raw_OPLS_fit_{opls_fit_q}_energies_in_Kelvin.txt"
            )

            # The combined GOMC and Gaussian dihedral angles and energies in kcal/mol
            gomc_gaussian_kcal_per_mol_energy_fitted_filename = f"all_normalized_OPLS_fit_{opls_fit_q}_energies_in_kcal_per_mol.txt"

            control_file_name_fitted_str = (
                f"GOMC_OPLS_fit_{opls_fit_q}_dihedral_coords_{scan_iter_q}.conf"
            )
            output_name_control_fitted_file_name_str = f"output_GOMC_OPLS_fit_{opls_fit_q}_dihedral_coords_{scan_iter_q}.txt"

            # "Started: Writing NVT GOMC control file for the GOMC simulation with fitted dihedral k values."
            gomc_control.write_gomc_control_file(
                charmm,
                f"{gomc_runs_folder_name}/{control_file_name_fitted_str}",
                "NVT",
                MC_steps,
                temperature_unyt_units,
                ff_psf_pdb_file_directory=None,
                check_input_files_exist=False,
                Parameters=f"{output_gomc_pdb_psf_ff_file_name_str}_OPLS_fit_{opls_fit_q}_dihedral.inp",
                Restart=True,
                Checkpoint=False,
                ExpertMode=False,
                Coordinates_box_0=f"{output_gomc_pdb_psf_ff_file_name_str}.pdb",
                Structure_box_0=f"{output_gomc_pdb_psf_ff_file_name_str}.psf",
                binCoordinates_box_0=read_gomc_fitted_restart_file_coor_dir_and_name,
                extendedSystem_box_0=read_gomc_fitted_restart_file_xsc_dir_and_name,
                binVelocities_box_0=None,
                Coordinates_box_1=None,
                Structure_box_1=None,
                binCoordinates_box_1=None,
                extendedSystem_box_1=None,
                binVelocities_box_1=None,
                input_variables_dict={
                    "PRNG": seed_no,
                    "Pressure": None,
                    "Ewald": True,
                    "ElectroStatic": True,
                    "VDWGeometricSigma": VDWGeometricSigma,
                    "Rcut": Rcut,
                    "RcutLow": RcutLow,
                    "LRC": LRC,
                    "Exclude": Exclude,
                    "DisFreq": DisFreq,
                    "VolFreq": VolFreq,
                    "RotFreq": RotFreq,
                    "RegrowthFreq": RegrowthFreq,
                    "IntraSwapFreq": IntraSwapFreq,
                    "CrankShaftFreq": CrankShaftFreq,
                    "SwapFreq": SwapFreq,
                    "OutputName": output_name_control_fitted_file_name_str,
                    "EqSteps": EqSteps,
                    "AdjSteps": AdjSteps,
                    "PressureCalc": output_true_list_input,
                    "RestartFreq": output_false_list_input,
                    "CheckpointFreq": output_false_list_input,
                    "ConsoleFreq": output_true_list_input,
                    "BlockAverageFreq": output_false_list_input,
                    "HistogramFreq": output_false_list_input,
                    "CoordinatesFreq": output_false_list_input,
                    "DCDFreq": output_false_list_input,
                    "CBMC_First": 12,
                    "CBMC_Nth": 10,
                    "CBMC_Ang": 50,
                    "CBMC_Dih": 50,
                },
            )
            # "Completed: NVT GOMC control file written for the GOMC simulation with fitted dihedral k values."

            # **************************************************************
            # make the GOMC control file for each dihedral angle (END)
            # **************************************************************

            # Change the units of the dihedral angles constants (k's) to the correct units
            # based on their potential form
            if charmm.utilized_NB_expression in ["LJ"]:
                k_constant_units_str = "kcal/mol"
            elif charmm.utilized_NB_expression in ["Mie", "Exp6"]:
                k_constant_units_str = "K"
            else:
                raise ValueError(
                    "ERROR: The non-bonded equation type is not the LJ, Mie or Exp6 "
                    "potential, which are the only available non-bonded equation potentials."
                )

            opls_k_constant_fitted_q_list_kcal_per_mol = [
                opls_fit_data_k0_kcal_per_mol_list[opls_q],
                opls_fit_data_k1_kcal_per_mol_list[opls_q],
                opls_fit_data_k2_kcal_per_mol_list[opls_q],
                opls_fit_data_k3_kcal_per_mol_list[opls_q],
                opls_fit_data_k4_kcal_per_mol_list[opls_q],
            ]

            # Add the modified k values for the simulation in their correct units,
            # kcal/mol for LJ and K for Mie or Exp6
            opls_k_constant_fitted_q_list_for_ff_modifications = [
                u.unyt_quantity(ki, "kcal/mol").to_value(
                    k_constant_units_str, equivalence="thermal"
                )
                for ki in opls_k_constant_fitted_q_list_kcal_per_mol
            ]

            mdf_frw.change_gomc_ff_file_dihedral_values(
                f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_per_xml.inp",
                f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_OPLS_fit_{opls_fit_q}_dihedral.inp",
                fit_dihedral_atom_types,
                fit_dihedral_opls_k_0_1_2_3_4_values=opls_k_constant_fitted_q_list_for_ff_modifications,
                zero_dihedral_atom_types=zero_dihedral_atom_types,
            )

            # **************************************************************
            # Run GOMC and get the system/fragment energy (START)
            # **************************************************************

            # only run NVT as we only want/need the initial energy and want the box a constant size
            run_gomc_command = (
                f"cd {gomc_runs_folder_name} && {gomc_binary_path}/GOMC_CPU_NVT +p{gomc_cpu_cores} "
                f"{control_file_name_fitted_str} > {output_name_control_fitted_file_name_str}"
            )

            exec_gomc_run_command = subprocess.Popen(
                run_gomc_command, shell=True, stderr=subprocess.STDOUT
            )

            os.wait4(
                exec_gomc_run_command.pid, os.WSTOPPED
            )  # pauses python until exec_gomc_run_command sim done

            # **************************************************************
            # Run GOMC and get the system/fragment energy (END)
            # **************************************************************

            # open the gomc raw energy file so it can be writen over in the loop
            read_gomc_log_fitted_file = open(
                f"{gomc_runs_folder_name}/{output_name_control_fitted_file_name_str}",
                "r",
            ).readlines()

            get_e_titles = True
            for log_file_iter, log_file_line_iter in enumerate(
                read_gomc_log_fitted_file
            ):
                log_file_splitline_iter = log_file_line_iter.split()

                # scan_iter starts at 1
                dihedral_angle_degrees = gaussian_raw_degrees_list[
                    scan_iter_q - 1
                ]

                # only open the gomc raw energy file and write header for 1st iteration (1)
                if len(log_file_splitline_iter) >= 2:
                    if (
                        scan_iter_q == 1
                        and log_file_splitline_iter[0] == "ETITLE:"
                        and get_e_titles is True
                    ):
                        gomc_combined_raw_fitted_energy = open(
                            gomc_raw_energy_fitted_filename, "w"
                        )
                        # remove the wrongly entered 5 spaces in before "ETITLE:
                        # (This will be fixed in GOMC so it is not required)
                        extra_spaces_for_header_space_gomc_bug = 5
                        gomc_combined_raw_fitted_energy.write(
                            f'{"Dihedral_Position": <19} '
                            f'{"Dihedral_Degrees": <19} '
                            f"{log_file_line_iter[extra_spaces_for_header_space_gomc_bug:]}"
                        )
                        if get_e_titles is True:
                            get_e_titles = False

                    # get only the initial configuration energy lin (i.e., step 0 line)
                    if (
                        log_file_splitline_iter[0] == "ENER_0:"
                        and log_file_splitline_iter[1] == "0"
                    ):
                        gomc_combined_raw_fitted_energy.write(
                            f"{scan_iter_q: <19}"
                            f"{dihedral_angle_degrees: <19} "
                            f"{log_file_line_iter[5:]}"
                        )

        # This fails when the GOMC simulations did not run
        try:
            # close the gomc_combined_raw_energy file
            gomc_combined_raw_fitted_energy.close()

        except:
            raise ValueError(
                "ERROR: The GOMC fit test simulations did not run. "
                "There is likely an error in creating the "
                "required GOMC test simulations files or user inputs to the desired files."
            )

        # extract the raw data
        GOMC_data_fitted_df = pd.DataFrame(
            pd.read_csv(gomc_raw_energy_fitted_filename, sep="\s+")
        )
        GOMC_data_fitted_dihedral_degrees_list = GOMC_data_fitted_df.loc[
            :, "Dihedral_Degrees"
        ].tolist()
        GOMC_data_fitted_total_energy_K_list = GOMC_data_fitted_df.loc[
            :, "TOTAL"
        ].tolist()

        # convert from Kelvin to kcal/mol normalize so the min value is 0
        GOMC_data_fitted_total_energy_kcal_per_mol_list = [
            i * conversion_K_to_kcal_per_mol
            for i in GOMC_data_fitted_total_energy_K_list
        ]
        GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list = [
            i - min(GOMC_data_fitted_total_energy_kcal_per_mol_list)
            for i in GOMC_data_fitted_total_energy_kcal_per_mol_list
        ]

        # *********************************
        # get GOMC data (END)
        # *********************************

        # *********************************
        # get Gaussian data (START)
        # *********************************
        # extract the raw data
        gaussian_data_df = pd.DataFrame(
            pd.read_csv(qm_energy_file_dir_and_name, sep="\s+", header=3)
        )
        gaussian_data_total_energy_Hartree_list = gaussian_data_df.iloc[
            :, 1
        ].tolist()

        # convert from Hartree to kcal/mol energy units
        gaussian_data_total_energy_kcal_per_mol_list = [
            i * conversion_hartree_to_kcal_per_mol
            for i in gaussian_data_total_energy_Hartree_list
        ]

        # normalize so the min value is 0
        gaussian_data_total_energy_kcal_per_mol_normalize_list = [
            i - min(gaussian_data_total_energy_kcal_per_mol_list)
            for i in gaussian_data_total_energy_kcal_per_mol_list
        ]

        # get the Gaussian minus GOMC total energy and then it normalized
        Gaussian_minus_GOMC_data_fitted_dihedral_degrees_list = (
            GOMC_data_fitted_dihedral_degrees_list
        )
        Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_list = [
            gaussian_data_total_energy_kcal_per_mol_normalize_list[i]
            - GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list[i]
            for i in range(
                0, len(gaussian_data_total_energy_kcal_per_mol_normalize_list)
            )
        ]

        Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_normalized_list = [
            Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_list[i]
            - min(
                Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_list
            )
            for i in range(
                0,
                len(
                    Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_list
                ),
            )
        ]

        # get R**2 for the fit, running through GOMC to get the new energy of the
        # individual fit.
        opls_r_squared_fitted_data_via_gomc_iter = mdf_math.get_r_squared(
            gaussian_data_total_energy_kcal_per_mol_normalize_list,
            GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list,
        )
        opls_r_squared_fitted_data_via_gomc_list.append(
            opls_r_squared_fitted_data_via_gomc_iter
        )

        # Use R**2 (R-squared) to compare the data used to fit the dihedral(s)
        # via a single or multi-dihedral fit, to the individual fit entered in
        # GOMC and then recompared to Gaussian and write out the data to a file.
        for q_angle in range(0, len(GOMC_data_fitted_dihedral_degrees_list)):
            if q_angle == 0:
                # write out the GOMC and Gaussian data in a file
                gomc_fitted_gaussian_kcal_per_mol_energy_data_txt_file = open(
                    gomc_fitted_gaussian_kcal_per_mol_energy_filename, "w"
                )
                gomc_fitted_gaussian_kcal_per_mol_energy_data_txt_file.write(
                    f"{'Dihedral_Degrees': <30} "
                    f"{'GOMC_E_kcal_per_mol': <30} "
                    f"{'Gaussian_E_kcal_per_mol': <30} "
                    f"{'Gaussian_minus_GOMC_E_kcal_per_mol': <40} "
                    f"{'k0_OPLS_kcal_per_mol': <30} "
                    f"{'k1_OPLS_kcal_per_mol': <30} "
                    f"{'k2_OPLS_kcal_per_mol': <30} "
                    f"{'k3_OPLS_kcal_per_mol': <30} "
                    f"{'k4_OPLS_kcal_per_mol': <30} "
                )

            gomc_fitted_gaussian_kcal_per_mol_energy_data_txt_file.write(
                f"\n{Gaussian_minus_GOMC_data_fitted_dihedral_degrees_list[q_angle]: <30} "
                f"{GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list[q_angle]: <30} "
                f"{gaussian_data_total_energy_kcal_per_mol_normalize_list[q_angle]: <30} "
                f"{Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_normalized_list[q_angle]: <40} "
                f"{str(opls_k_constant_fitted_q_list_kcal_per_mol[0]): <30} "
                f"{str(opls_k_constant_fitted_q_list_kcal_per_mol[1]): <30} "
                f"{str(opls_k_constant_fitted_q_list_kcal_per_mol[2]): <30} "
                f"{str(opls_k_constant_fitted_q_list_kcal_per_mol[3]): <30} "
                f"{str(opls_k_constant_fitted_q_list_kcal_per_mol[4]): <30} "
            )

        # Compare original fit vs run through EACH GOMC as a validation test case
        if opls_q == 0:
            opls_fit_acceptable_r_squared_values_list = []
            opls_fit_acceptable_r_squared_values_not_in_rtol_list = []
            opls_fit_data_non_zero_k_constants_not_acceptable_r_squared_values_list = (
                []
            )
            opls_fit_data_r_squared_not_acceptable_r_squared_values_list = []
            opls_r_squared_fitted_data_via_gomc_list_not_acceptable_r_squared_values_list = (
                []
            )
            max_opls_fit_Rsquared = None
            max_opls_fit_Rsquared_list_index_number = None

        if opls_fit_data_r_squared_list[opls_q] >= r_squared_min:
            opls_fit_acceptable_r_squared_values_list.append(1)

            if np.isclose(
                opls_r_squared_fitted_data_via_gomc_list[opls_q],
                opls_fit_data_r_squared_list[opls_q],
                rtol=r_squared_rtol,
            ):
                opls_fit_acceptable_r_squared_values_not_in_rtol_list.append(1)

                # get maximum value for  fit vs run through EACH GOMC as a validation test case
                if (
                    max_opls_fit_Rsquared is None
                    and max_opls_fit_Rsquared_list_index_number is None
                ):
                    max_opls_fit_Rsquared = opls_fit_data_r_squared_list[opls_q]
                    max_opls_fit_Rsquared_list_index_number = opls_q

                elif (
                    opls_fit_data_r_squared_list[opls_q] > max_opls_fit_Rsquared
                ):
                    max_opls_fit_Rsquared = opls_fit_data_r_squared_list[opls_q]
                    max_opls_fit_Rsquared_list_index_number = opls_q

            else:
                opls_fit_acceptable_r_squared_values_not_in_rtol_list.append(0)

                opls_fit_data_non_zero_k_constants_not_acceptable_r_squared_values_list.append(
                    opls_fit_data_non_zero_k_constants_list[opls_q]
                )
                opls_fit_data_r_squared_not_acceptable_r_squared_values_list.append(
                    mdf_math.round_to_sig_figs(
                        opls_fit_data_r_squared_list[opls_q], sig_figs=8
                    )
                )
                opls_r_squared_fitted_data_via_gomc_list_not_acceptable_r_squared_values_list.append(
                    mdf_math.round_to_sig_figs(
                        opls_r_squared_fitted_data_via_gomc_list[opls_q],
                        sig_figs=8,
                    )
                )

    # round the 'opls_fit_data_r_squared_list' and 'opls_r_squared_fitted_data_via_gomc_list' for printing
    rounded_opls_fit_data_r_squared_list = [
        mdf_math.round_to_sig_figs(r_i, sig_figs=8)
        for r_i in opls_fit_data_r_squared_list
    ]
    rounded_opls_r_squared_fitted_data_via_gomc_list = [
        mdf_math.round_to_sig_figs(r_i, sig_figs=8)
        for r_i in opls_r_squared_fitted_data_via_gomc_list
    ]

    # combine in pairs for the dihedrals [constants_used, r_squared_fitted, r_squared_rerun_GOMC, tolerance]
    rounded_opls_combined_Rsq_tol_list = []
    # Get only ones that meet R^2 objective but not tolerance
    meet_r_sq_not_tol_rounded_opls_combined_Rsq_tol_list = []
    for r_n_tol_i in range(0, len(rounded_opls_fit_data_r_squared_list)):
        rounded_opls_combined_Rsq_tol_list.append(
            [
                opls_fit_data_non_zero_k_constants_list[r_n_tol_i],
                rounded_opls_fit_data_r_squared_list[r_n_tol_i],
                rounded_opls_r_squared_fitted_data_via_gomc_list[r_n_tol_i],
                mdf_math.round_to_sig_figs(
                    abs(
                        rounded_opls_fit_data_r_squared_list[r_n_tol_i]
                        - rounded_opls_r_squared_fitted_data_via_gomc_list[
                            r_n_tol_i
                        ]
                    ),
                    sig_figs=8,
                ),
            ]
        )

        # Get only ones that meet R^2 objective but not tolerance
        if opls_fit_data_r_squared_list[r_n_tol_i] >= r_squared_min:

            if not np.isclose(
                opls_r_squared_fitted_data_via_gomc_list[r_n_tol_i],
                opls_fit_data_r_squared_list[r_n_tol_i],
                rtol=r_squared_rtol,
            ):
                meet_r_sq_not_tol_rounded_opls_combined_Rsq_tol_list.append(
                    [
                        opls_fit_data_non_zero_k_constants_list[r_n_tol_i],
                        rounded_opls_fit_data_r_squared_list[r_n_tol_i],
                        rounded_opls_r_squared_fitted_data_via_gomc_list[
                            r_n_tol_i
                        ],
                        mdf_math.round_to_sig_figs(
                            abs(
                                rounded_opls_fit_data_r_squared_list[r_n_tol_i]
                                - rounded_opls_r_squared_fitted_data_via_gomc_list[
                                    r_n_tol_i
                                ]
                            ),
                            sig_figs=8,
                        ),
                    ]
                )

    # print some Info for the fit, which is the same for all the dihedral angles
    print("********************")
    print(f"charmm.combining_rule = {charmm.combining_rule}")
    print(f"charmm.utilized_NB_expression = {charmm.utilized_NB_expression}")
    print("********************")
    print(
        f"opls_r_squared_fitted_data_via_gomc_list = {opls_r_squared_fitted_data_via_gomc_list}"
    )
    print(f"opls_fit_data_r_squared_list = {opls_fit_data_r_squared_list}")
    print(f"user set 'r_squared_min' = {r_squared_min}")
    print(f"user set 'r_squared_rtol' = {r_squared_rtol}")
    print("********************")
    print(
        f"dihedral_degrees_list = {np.asarray(sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list)}"
    )
    print(f"normalized QM - MM total_energy (kcal_per_mol) = ")
    print(
        f"{np.asarray(sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)}"
    )
    print("********************")

    # Compare original fit vs run through ALL GOMC as a validation test case
    if np.sum(opls_fit_acceptable_r_squared_values_not_in_rtol_list) == 0:
        raise ValueError(
            f"ERROR: The calculated R-squared energy values from the fit types "
            f"do not match the validated case for 'r_squared_min' >= "
            f"{mdf_math.round_to_sig_figs(r_squared_min,sig_figs=8)}, "
            f"within the relative tolerance or 'r_squared_rtol' = "
            f"{mdf_math.round_to_sig_figs(r_squared_rtol,sig_figs=8)}. \n"
            f"Constants_used = The calculated R-squared energy values from the fit type constants_used."
            f"R-squared_fitted = Gaussian minus GOMC with the selected dihedral set to zero \n"
            f"R-squared_new_dihedral = Gaussian minus GOMC with the selected individual dihedral added in GOMC \n"
            f"Abs(delta) = absolute_value(R-squared_fitted - R-squared_new_dihedral) \n"
            f"- The fits for all are shown here for all the dihedral combinations \n"
            f"[opls_constants_used, R-squared_fitted, R-squared_new_dihedral, Abs(delta)] \n"
            f"{rounded_opls_combined_Rsq_tol_list}, \n"
            f"The 'r_squared_min' and 'r_squared_rtol' "
            f"variables may need to be adjusted, \n"
            f"there is likely something wrong with the fitting procedure, the "
            f"software parameters need tuned, or there is a bug in the software. \n\n "
            f"NOTE: Since the R-squared values are calculated via different parameters, \n"
            f"the compared R-squared values could be very different if they are not nearly \n"
            f"a perfect fit (R-squared --> ~0.98 to 0.99999999)."
        )

    elif np.sum(
        opls_fit_acceptable_r_squared_values_not_in_rtol_list
    ) != np.sum(opls_fit_acceptable_r_squared_values_list):
        warn(
            f"WARNING: The calculated R-squared energy values from the fit types "
            f"do match the validated case for 'r_squared_min' >= "
            f"{mdf_math.round_to_sig_figs(r_squared_min,sig_figs=8)}, "
            f"but do not fit within the relative tolerance of 'r_squared_rtol' = "
            f"{mdf_math.round_to_sig_figs(r_squared_rtol,sig_figs=8)}. \n"
            f"Constants_used = The calculated R-squared energy values from the fit type constants_used. \n"
            f"R-squared_fitted = Gaussian minus GOMC with the selected dihedral set to zero. \n"
            f"R-squared_new_dihedral = Gaussian minus GOMC with the selected individual dihedral added in GOMC. \n"
            f"Abs(delta) = Abs(R-squared_fitted - R-squared_new_dihedral). \n"
            f"- The fits for all are shown here for all the dihedral combinations \n"
            f"[opls_constants_used, R-squared_fitted, R-squared_new_dihedral, Abs(delta)] \n"
            f"{meet_r_sq_not_tol_rounded_opls_combined_Rsq_tol_list}. \n"
            f"The 'r_squared_min' and 'r_squared_rtol' "
            f"variables may need to be adjusted, \n"
            f"there is likely something wrong with the fitting procedure, the "
            f"software parameters need tuned, or there is a bug in the software. \n\n "
            f"NOTE: Since the R-squared values are calculated via different parameters, \n"
            f"the compared R-squared values could be very different if they are not nearly \n"
            f"a perfect fit (R-squared --> ~0.98 to 0.99999999)."
        )

    gomc_fitted_gaussian_kcal_per_mol_energy_data_txt_file.close()

    # *********************************
    # Check the all the OPLS dihedral forms are correct
    # by running GOMC with the fitted values and comparing it to QM
    # (END)
    # *********************************
    fig2, (ax2) = plt.subplots(1)
    axis_Label_font_size = 12
    legend_font_size = 12
    ax2.set_xlabel("phi (degrees)", size=axis_Label_font_size)
    ax2.set_ylabel("Dihedral Energy (kcal/mol)", size=axis_Label_font_size)

    plt.title(
        "OPLS: All viable individual dihedral fits. \n "
        "(Non-zero k's = 1 and 3 -> label k_non_0='1_3')"
    )

    # loop thru the different 'non_zero_k_constants' for the OPLS dihedral
    max_energy_list_iter = []
    for opls_fit_j in range(0, len(opls_fit_data_non_zero_k_constants_list)):
        # list the phi values to check
        delta_phi_degrees_iter = 1
        phi_degrees_iter = int(360 / delta_phi_degrees_iter)

        phi_degrees_list_iter = []
        cos_power_list_iter = []
        opls_dihedral_energy_list_iter = []
        for j in range(0, phi_degrees_iter):
            phi_degrees_list_iter.append(-180 + j * delta_phi_degrees_iter)

            cos_power_list_iter.append(
                opls_fit_data_non_zero_k_constants_list[opls_fit_j]
            )

            opls_dihedral_energy_iter = mdf_math.opls_dihedral(
                (
                    cos_power_list_iter[-1],
                    phi_degrees_list_iter[-1],
                    None,
                    None,
                    None,
                    None,
                    None,
                ),
                opls_fit_data_k0_kcal_per_mol_list[opls_fit_j],
                opls_fit_data_k1_kcal_per_mol_list[opls_fit_j],
                opls_fit_data_k2_kcal_per_mol_list[opls_fit_j],
                opls_fit_data_k3_kcal_per_mol_list[opls_fit_j],
                opls_fit_data_k4_kcal_per_mol_list[opls_fit_j],
            )

            opls_dihedral_energy_list_iter.append(opls_dihedral_energy_iter)

        max_energy_list_iter.append(max(opls_dihedral_energy_list_iter))

        plt.plot(
            np.asarray(phi_degrees_list_iter),
            np.asarray(opls_dihedral_energy_list_iter),
            "-",
            label=f" {opls_fit_data_non_zero_k_constants_list[opls_fit_j]} "
            f"| $R^{2}$={np.round(opls_r_squared_fitted_data_via_gomc_list[opls_fit_j], decimals=4)}",
        )

    ax2.set_xticks(major_xticks)
    ax2.set_xticks(minor_xticks, minor=True)
    ax2.set_xlim(-180, 180 + 0.001)

    plt.legend(
        ncol=2,
        loc="lower center",
        fontsize=legend_font_size,
        prop={"size": legend_font_size},
    )

    ax2.set_xticks(major_xticks)
    ax2.set_xticks(minor_xticks, minor=True)
    ax2.set_xlim(-180, 180 + 0.001)

    # get max energy
    plot_max = int(max(max_energy_list_iter) + 1.51)
    ax2.set_ylim(-int(plot_max * 2.5), plot_max)

    # plt.show()
    fig2.savefig(
        f"opls_{all_individual_fit_dihedral_k_constants_figure_filename}",
        dpi=300,
    )
    # **************************************************************
    # **************************************************************
    # **************************************************************
    # **************************************************************
    # Use the existing Gaussian (QM) file for their energy data.
    # Build the system using a mol2 file via mBuild and MoSDeF.
    # Use GOMC to run the MoSDeF build with the .coor and .xsc files
    # to yield the exact configuration of the QM equilibrium fragments.
    # Note: the MoSDeF XML file must set the target dihedral to zero (0) energy.
    # Subract the QM - GOMC eneries, normalize the dihedral to min of zero (0) energy,
    # fit the various dihedral combinations to standard OPLS form (no f0 constant),
    # then convert each fit to the periodic/CHARMM dihedral and RB torsions forms.
    # (END)
    # **************************************************************
    # **************************************************************
    # **************************************************************
    # **************************************************************
