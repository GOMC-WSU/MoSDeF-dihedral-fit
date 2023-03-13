import os
import pandas as pd
import numpy as np
import mbuild as mb
import subprocess
import unyt as u
import mosdef_gomc.formats.gmso_charmm_writer as mf_charmm
import mosdef_gomc.formats.gmso_gomc_conf_writer as gomc_control
import matplotlib.pyplot as plt
import math

from scipy.optimize import curve_fit
from mbuild.utils.conversion import OPLS_to_RB
from mosdef_gomc.utils.conversion import OPLS_to_periodic

import mosdef_dihedral_fit.utils.file_read_and_write as mdf_frw
import mosdef_dihedral_fit.utils.math_operations as mdf_math

import warnings
from warnings import warn
warnings.filterwarnings('ignore')


def fit_dihedral_with_gomc(
        fit_dihedral_atom_types,
        mol2_selection,
        forcefield_selection,
        temperature,
        gomc_binary_path,
        qm_log_files_and_entries_to_remove,
        zeroed_dihedral_atom_types=None,
        qm_engine="gaussian",
        override_VDWGeometricSigma=None,
        atom_type_naming_style='general',
        gomc_cpu_cores=1,
        fit_min_validated_r_squared=0.99,
        fit_validation_r_squared_rtol=1e-03
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

    NOTE: This dihedral fitting process can accomidate fitting more
    than 1 dihedral fit of the same type simultaneously.

    NOTE: Not all cos power terms are allowed to be utilized in the
    dihedral fits. This function tests if a cos power is able to be
    used for the fit, and only uses the valid power in the fitting
    process.

    NOTE: The 'extracted_guassian_data' and 'GOMC_simulations'
    folder are deleted at the beginning of this function,
    and recreated while running this function to ensure only the
    lasted data is in these folders.

    Parameters
    ----------
    fit_dihedral_atom_types: list of four (4) strings (Example: ['HC', 'CT, 'CT, 'HC'])
        The atom types/classes (strings in the list) of the dihedral which is
        being fitted with non-zero k-values.
    mol2_selection: str
        The mol2 file which matches the element, atom type, bonded connnections,
        the 'EXACT ATOM ORDER AND CONFIGURATION AS IN THE QM SIMULATION INPUT FILES'.
        This is required to know the MM bonding in the atoms, because QM simulations
        do not explictly specify the system bonds.
    forcefield_selection: str
        Apply a foyer or gmso forcefield to the output file by selecting a
        force field XML file with its path or by using the standard force
        field name provided the `foyer` package. This force field file must
        contain all the non-bonded, bonds, angles, dihedrals parameters for
        the system.

        * Example str for FF file: 'path_to file/trappe-ua.xml'

    temperature: unyt.unyt_quantity
        The temperature of the system that was performed for the Quantum Mechanics
        (QM) simulation.
    gomc_binary_path: str
        The path or directory of the GOMC binary file "GOMC_CPU_NVT", which is used to
        perform the Molecular Mechanics (MM) energy calculations. This does not include
        the "GOMC_CPU_NVT" in this variable.

        Example: '/home/brad/Programs/GOMC/GOMC_2_76/bin'

    qm_log_files_and_entries_to_remove_dict: dict, {str: [int, ..., int]}
        This is a dictionary comprised of a key (string) of the QM log file path and name,
        and a list of integers, which are the QM optimization parameters to remove from
        the written data, in order of reading from each file. These can be seen in the
        order of the dictionary file name (strings).  These removed parameters allow
        users to remove any bad or repeated data points for the QM log file when needed.

        Example 1: {'path/guassian_log_file.log': []}
        Uses all the optimized data points from the 'path/guassian_log_file.log' file.

        Example 2: {'path/guassian_log_file.log': [0, 23]}
        Uses all data points from the 'path/guassian_log_file.log' file, except points
        0 and 23.  NOTE: Python counting starts at 0.
    zeroed_dihedral_atom_types: nest list with lists of four (4) strings, default=None
        The nests list(s) of the other dihedrals, that need to have their k-values zeroed to
        properly fit the the 'fit_dihedral_atom_types' dihedral.

        Example: [['CT', 'CT, 'CT, 'HC'], ['NT', 'CT, 'CT, 'HC']]

    qm_engine: str (currently only 'guassian'), default='guassian'
        The Quantum Mechanics (QM) simulation engine utilized to produce the files listed
        in the 'qm_log_files_and_entries_to_remove' variable(s).
    override_VDWGeometricSigma: boolean, default = None
        Override the VDWGeometricSigma in the foyer or GMSO XML file.
        If this is None, it will use whatever is specified in the XML file, or the
        default foyer or GMSO values. BEWARE, if it is not specified XML file, it has a default.
        If this is None, it will use whatever is specified in the XML file,
        or the default foyer or GMSO values.

        The VDWGeometricSigma is the geometric mean, as required by OPLS force field,
        to combining Lennard-Jones sigma parameters for different atom types.
        If set to True, GOMC uses geometric mean to combine Lennard-Jones or VDW sigmas.
        Note: The default setting is pulled from the force field XML, if not present it
        defaults to geometric via MoSDeF's default setting of VDWGeometricSigma is True,
        which uses the arithmetic mean when combining Lennard-Jones or "
        VDW sigma parameters for different atom types.
    atom_type_naming_style: str, optional, default='all_unique', ('general' or 'all_unique')
        'general':

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

        Example of CHARMM style atom types in an all-atom ethane and ethanol system:
        * Ethane: alkane carbon = CT, alkane hydrogen = HC
        * Ethanol: alkane carbon = CT, alkane hydrogen = HC , oxygen in alcohol = OH, hydrogen in alcohol = OH

        This is only permitted when the following is true; otherwise it will default to the the 'all_unique':
        * All the MoSDeF force field XML's atom classes' non-bonded parameters
        (sigma, epsilon, mass, and Mie-n power constant) values are THE SAME.
        * If the general CHARMM style atom type in any residue/molecule's gomc_fix_bonds_angles,
        gomc_fix_bonds, or gomc_fix_angles NOT IN any other residue/molecule, the 'all_unique' type
        will be used.

        'all_unique':
        The 'all_unique' convention is the SAFE way to parameterize the system.
        The MoSDeF force field XML atom names within residue/molecule are all unique,
        where each one adds an alpha numberic value after the MoSDeF force field XML atom classes to
        ensure uniqueness within the within residue/molecule.
        The OPLS atom types/name do not require all the sigma, epsilon, mass values to be the same,
        but have less bonded class parameters.

        Example of CHARMM style atom types in an all-atom ethane and ethanol system:
        * Ethane: alkane carbon type 0 = CT0, alkane hydrogen type 0 = HC0
        * Ethanol: alkane carbon type 1 = CT1, alkane carbon type 2 = CT2,
        alkane hydrogen type 1 = HC1 , oxygen in alcohol type 0 = OH0, hydrogen in alcohol type 0 = OH0

        This is selected when auto-selected when:
        * All the MoSDeF force field XML's atom classes' non-bonded parameters
        (sigma, epsilon, mass, and Mie-n power constant) values are NOT THE SAME.
        * If the general CHARMM style atom type in any residue/molecule's gomc_fix_bonds_angles,
        gomc_fix_bonds, or gomc_fix_angles are IN any other residue/molecule.
    gomc_cpu_cores: int, default=1
        The number of CPU-cores that are used to perform the GOMC simulations, required
        for the Molecular Mechanics (MM) energy calulations.
    fit_min_validated_r_squared: float (0 <= fit_min_validated_r_squared <= 1), default=0.99
        The minimum R**2 (R-squared) value to test the validity of the fit with the
        new dihedral fitted constants, as fitted in the
        QM - MM energy data vs. the dihedral function fit, mentioned below.

        This mean that any R**2 (R-squared) value in the
        'opls_all_summed_dihedrals_k_constants_figure.pdf' plot less than this value
        will not be check for accuracy because the  R**2 (R-squared) values are compared
        differently as a check.  These are compared between the following calculations:

        * QM - MM energy data vs. the dihedral function fit:
        For the MM calculations, the 'fit_dihedral_atom_types' and
        'zeroed_dihedral_atom_types' are dihedral energies are set to zero, which is
        during the fitting process with 1 or more of the same dihedrals being fit simultaneously.

        * QM vs. the MM energy data:
        For the MM calculations, the 'fit_dihedral_atom_types' are set to the values which were
        fit for the specific cosine combinations during the fitting process with 1 or more of the
        same dihedrals being fit simultaneously, and the 'zeroed_dihedral_atom_types' are
        dihedral energies are set to zero.

    fit_validation_r_squared_rtol: float, default=1e-03
        Where the QM data is defined as the actual data; this is the fractional difference
        of the dihedral's calculated R-squared values between:
        * The QM-MM fitting process, where the fit MM dihedral k-values are zero (0).
        * The MM calculations where the fit k-value are entered in the MM data and
        compared to the QM data.

        fit_dihedral_atom_types,
        mol2_selection,
        forcefield_selection,
        temperature,
        gomc_binary_path,
        qm_log_files_and_entries_to_remove,
        zeroed_dihedral_atom_types=None,

    Returns
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
            ('zeroed_dihedral_atom_types') being zeroed out, where the XXX
            in the file name being the different  possibilities of cos power combinations
            of the dihedral fit. These may be 1 or more of these files/cos power
            combinations.
        GOMC_simulations/GOMC_zeroed_dihedral_coords_XXX.conf
            The GOMC control files for all the phi dihedral angles selected
            from the QM simulations where the 'fit_dihedral_atom_types' and the
            'zeroed_dihedral_atom_types' k-values are all set to zero.  The XXX is the
            integer number of the dihedrals starting at 1, adding 1 for every
            additional phi dihedral angles selected from the QM simulations.

            These are simulations are to get the total energy difference between
            QM and MM simulation with all the 'fit_dihedral_atom_types' and the
            'zeroed_dihedral_atom_types' k-values are all set to zero.
        GOMC_simulations/GOMC_OPLS_fit_YYY_dihedral_coords_XXX.conf
            The GOMC control files for all the phi dihedral angles selected
            from the QM simulations where the 'fit_dihedral_atom_types' are set to
            the solved k-values for the cos power equation combination, which is listed
            as the variable YYY, and the 'zeroed_dihedral_atom_types' k-values are all set
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
        extracted_guassian_data/dihedral.txt
            The QM data in a Gaussian-style output file which holds the scanned
            dihedral angles, in degrees, and the optimized energy value, in Hartree units,
            for the molecule/fragment.
        extracted_guassian_data/dihedral_coords_position_XXX.txt
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
            'zeroed_dihedral_atom_types' k-values are all set to zero. This also
            contains the sum of all the C1, C2, C3, and C4 values for every
            'fit_dihedral_atom_types' in the system, where there may be multiple of the
            same dihedral in the same molecule.

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

            OPLS_energy =   k1 * (1 + cos(1 * phi))
                          + k2 * (1 - cos(2 * phi))
                          + k3 * (1 + cos(3 * phi))
                          + k4 * (1 - cos(4 * phi))

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
            ('zeroed_dihedral_atom_types') are zeroed out.
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
    # write the qm data files data out
    mdf_frw.write_qm_data_files(
        qm_log_files_and_entries_to_remove,
        qm_engine=qm_engine
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
    gomc_runs_folder_name = 'GOMC_simulations'

    if os.path.isdir(gomc_runs_folder_name):
        os.rmdir(gomc_runs_folder_name)
    os.mkdir(gomc_runs_folder_name)

    # The gomc raw energy filename in Kelvin Energies
    gomc_raw_energy_filename = "gomc_raw_energies_in_Kelvin.txt"

    # The combined GOMC and Gaussian dihedral angles and energies in kcal/mol
    gomc_gaussian_kcal_mol_energy_filename = "all_normalized_energies_in_kcal_per_mol.txt"

    output_gomc_pdb_psf_ff_file_name_str = f'GOMC_pdb_psf_ff_files'
    seed_no = 12345

    # load bis(ethylhexylamido) with Gd
    fragment = mb.load(mol2_selection, smiles=False)
    fragment.name = 'TMP'

    residues_list = [fragment.name]
    fix_residues = None
    fix_residue_in_box = None

    # Build the methanol liquid box_1

    # liquid_box_0_length_nm must be a value <= 999.8 nm and an interger value in angstroms <= 9998 Ang
    liquid_box_0_length_nm = 999.8
    print("INFO: Started Building the fragment for the GOMC simulation with dihedral k values = 0")

    box_0_liq = mb.fill_box(
        compound=[fragment],
        n_compounds=[1],
        box=[liquid_box_0_length_nm, liquid_box_0_length_nm, liquid_box_0_length_nm],
        seed=seed_no,
    )

    print("INFO: Finished Building the fragment for the GOMC simulation with dihedral k values = 0")

    print("INFO: Started Building the Charmm Object for the GOMC simulation with dihedral k values = 0")
    # build the charmm object
    charmm = mf_charmm.Charmm(box_0_liq,
                              f'{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}',
                              structure_box_1=None,
                              filename_box_1=None,
                              ff_filename=f'{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_per_xml',
                              forcefield_selection=forcefield_selection,
                              residues=residues_list,
                              bead_to_atom_name_dict=None,
                              fix_residue=fix_residues,
                              fix_residue_in_box=fix_residue_in_box,
                              gomc_fix_bonds_angles=None,
                              atom_type_naming_style=atom_type_naming_style
                              )
    print("INFO: Finished Building the Charmm Object for the GOMC simulation with dihedral k values = 0")

    # Write the write the FF (.inp), psf, pdb, and GOMC control files
    print("INFO: Started Writing the PDB, PSF, and FF files for the GOMC simulation with dihedral k values = 0")
    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()
    print("INFO: Finished Writing the PDB, PSF, and FF files for the GOMC simulation with dihedral k values = 0")

    # **************************************************************
    # make the PDB, PSF and FF files for all the dihedral angles (END)
    # **************************************************************

    # **************************************************************
    # make the GOMC control file with the selected dihedral angles set to zero (START)
    # **************************************************************
    mdf_frw.change_gomc_ff_file_dihedral_values(
        f'{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_per_xml.inp',
        f'{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_zeroed.inp',
        fit_dihedral_atom_types,
        fit_dihedral_opls_k_0_1_2_3_4_values=[0, 0, 0, 0, 0],
        zeroed_dihedral_atom_types=zeroed_dihedral_atom_types,
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
    if os.isdir(xyz_xsc_coor_files_directory):
        os.rmdir(xyz_xsc_coor_files_directory)
    os.mkdir(xyz_xsc_coor_files_directory)

    # write all the xyz coordinate from the Guassian optimized coordinate file in the 'xyz_files' folder
    [atom_pdb_names_list, elementpdb_names_list] = mdf_frw.get_atom_names_and_elements_from_pdb(
        f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}.pdb"
    )

    qm_energy_file_dir_and_name = 'extracted_guassian_data/dihedral.txt'
    qm_parital_coordinate_file_starting_dir_and_name = 'extracted_guassian_data/dihedral_coords_position_'
    qm_coordinate_file_extension = 'txt'

    # check the gaussian file is correct
    mdf_frw.check_guassian_angle_energy_file_correct(
        qm_energy_file_dir_and_name)

    # Read the gaussian data and extract angles and number of scans (number of angles and degress analyzed)
    Guassian_raw_degrees_list = pd.DataFrame(pd.read_csv(
        qm_energy_file_dir_and_name, sep='\s+', header=3)).iloc[:, 0].tolist()
    total_qm_scans = len(Guassian_raw_degrees_list)

    mdf_frw.write_xyz_file_from_gaussian_coordinates(
        elementpdb_names_list,
        qm_parital_coordinate_file_starting_dir_and_name,
        qm_coordinate_file_extension,
        xyz_xsc_coor_files_directory,
        total_qm_scans
    )

    # Using vmd write the GOMC restart .coor files required for the
    mdf_frw.write_restart_coor_from_xyz_file(
        xyz_xsc_coor_files_directory, total_qm_scans)

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
    for scan_iter in range(1, len(Guassian_raw_degrees_list) + 1):

        read_gomc_restart_file_coor_dir_and_name = \
            f'../{xyz_xsc_coor_files_directory}/dihedral_coords_position_{scan_iter}.coor'
        read_gomc_restart_file_xsc_dir_and_name = f'../{xyz_xsc_coor_files_directory}/starting_point.xsc'

        control_file_name_str = f'GOMC_zeroed_dihedral_coords_{scan_iter}.conf'
        output_name_control_file_name_str = f'output_GOMC_zeroed_dihedral_coords_{scan_iter}.txt'

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

        print("#**********************")
        print("Started: NVT GOMC control file writing for the GOMC simulation with dihedral k values = 0")
        print("#**********************")
        # calc MC steps for gomc equilb
        MC_steps = 2
        EqSteps = 1
        AdjSteps = 1
        output_true_list_input = [True, 1]
        output_false_list_input = [False, 1]

        print(f"charmm.combining_rule = {charmm.combining_rule}")

        gomc_control.write_gomc_control_file(
            charmm,
            f'{gomc_runs_folder_name}/{control_file_name_str}',
            'NVT',
            MC_steps,
            temperature,
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
                "VDWGeometricSigma": override_VDWGeometricSigma,
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
                "Potential": "VDW",
                "CBMC_First": 12,
                "CBMC_Nth": 10,
                "CBMC_Ang": 50,
                "CBMC_Dih": 50,
            },
        )
        print("#**********************")
        print("Completed: NVT GOMC control file written for the GOMC simulation with dihedral k values = 0")
        print("#**********************")

        # **************************************************************
        # make the GOMC control file for each dihedral angle (END)
        # **************************************************************

        # *********************************
        # Write the restart .xsc file for GOMC (START)
        # *********************************
        gomc_restart_xsc_txt_file = open(
            f'{xyz_xsc_coor_files_directory}/starting_point.xsc', "w")
        gomc_restart_xsc_txt_file.write(
            f"# GOMC extended system configuration output file\n")
        gomc_restart_xsc_txt_file.write(
            f"#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w\n"
        )

        gomc_restart_xsc_txt_file.write(f"1 {Rcut.to_value('angstrom') * 2} 0 0 "
                                        f"0 {Rcut.to_value('angstrom') * 2} 0 "
                                        f"0 0 {Rcut.to_value('angstrom') * 2} "
                                        f"{Rcut.to_value('angstrom')} "
                                        f"{Rcut.to_value('angstrom')} "
                                        f"{Rcut.to_value('angstrom')} "
                                        f"0 0 0 "
                                        f"0 0 0"
                                        f"\n")
        gomc_restart_xsc_txt_file.close()
        # *********************************
        # Write the restart .xsc file for GOMC (END)
        # *********************************

        # **************************************************************
        # Run GOMC and get the system/fragment energy (START)
        # **************************************************************

        # only run NVT as we only want/need the initial energy and want the box a constant size
        run_gomc_command = \
            f"cd {gomc_runs_folder_name} && {gomc_binary_path}/GOMC_CPU_NVT +p{gomc_cpu_cores} " \
            f"{control_file_name_str} > {output_name_control_file_name_str}"

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
            f'{gomc_runs_folder_name}/{output_name_control_file_name_str}', "r").readlines()
        get_e_titles = True
        for log_file_iter, log_file_line_iter in enumerate(read_gomc_log_file):
            log_file_splitline_iter = log_file_line_iter.split()

            # scan_iter starts at 1
            dihedral_angle_degrees = Guassian_raw_degrees_list[scan_iter-1]

            # only open the gomc raw energy file and write header for 1st iteration (1)
            if len(log_file_splitline_iter) >= 2:
                if scan_iter == 1 and log_file_splitline_iter[0] == "ETITLE:" and get_e_titles is True:
                    gomc_combined_raw_energy = open(
                        gomc_raw_energy_filename, "w")
                    # remove the wrongly entered 5 spaces in before "ETITLE:
                    # (This will be fixed in GOMC so it is not required)
                    extra_spaces_for_header_space_gomc_bug = 5
                    gomc_combined_raw_energy.write(f'{"Dihedral_Position": <19} '
                                                   f'{"Dihedral_Degrees": <19} '
                                                   f'{log_file_line_iter[extra_spaces_for_header_space_gomc_bug:]}'
                                                   )
                    if get_e_titles is True:
                        get_e_titles = False

                # get only the initial configuration energy lin (i.e., step 0 line)
                if log_file_splitline_iter[0] == 'ENER_0:' and log_file_splitline_iter[1] == '0':
                    gomc_combined_raw_energy.write(f'{scan_iter: <19}'
                                                   f'{dihedral_angle_degrees: <19} '
                                                   f'{log_file_line_iter[5:]}'
                                                   )

    # This fails when the GOMC simulations did not run
    try:
        # close the gomc_combined_raw_energy file
        gomc_combined_raw_energy.close()

    except:
        raise ValueError("ERROR: The GOMC simulations did not run. There is likely an error in creating the "
                         "required GOMC files or user inputs to the desired files.")

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
    conversion_K_to_kcal_per_mol = (
        1 * u.Kelvin).to_value("kcal/mol", equivalence="thermal")

    # *********************************
    # get GOMC data (START)
    # *********************************
    # extract the raw data
    GOMC_data_df = pd.DataFrame(pd.read_csv(
        gomc_raw_energy_filename,  sep='\s+'))
    GOMC_data_dihedral_degrees_list = GOMC_data_df.loc[:, 'Dihedral_Degrees'].tolist(
    )
    GOMC_data_total_energy_K_list = GOMC_data_df.loc[:, 'TOTAL'].tolist()

    # convert from Kelvin to kcal/mol normalize so the min value is 0
    GOMC_data_total_energy_kcal_per_mol_list = \
        [i * conversion_K_to_kcal_per_mol for i in GOMC_data_total_energy_K_list]
    GOMC_data_total_energy_kcal_per_mol_normalize_list = \
        [i - min(GOMC_data_total_energy_kcal_per_mol_list)
         for i in GOMC_data_total_energy_kcal_per_mol_list]

    print(
        f"GOMC_data_dihedral_degrees_list = {GOMC_data_dihedral_degrees_list}")
    print(
        f"GOMC_data_total_energy_kcal_per_mol_list = {GOMC_data_total_energy_kcal_per_mol_list}")
    print(
        f"GOMC_data_total_energy_kcal_per_mol_normalize_list = {GOMC_data_total_energy_kcal_per_mol_normalize_list}")
    # *********************************
    # get GOMC data (END)
    # *********************************

    # *********************************
    # get Gaussian data (START)
    # *********************************

    # extract the raw data
    Guassian_data_df = pd.DataFrame(pd.read_csv(
        qm_energy_file_dir_and_name, sep='\s+', header=3))
    Guassian_data_dihedral_degrees_list = Guassian_data_df.iloc[:, 0].tolist()
    Guassian_data_total_energy_Hartree_list = Guassian_data_df.iloc[:, 1].tolist(
    )

    # convert from Hartree to kcal/mol energy units
    Guassian_data_total_energy_kcal_per_mol_list = \
        [i * conversion_hartree_to_kcal_per_mol for i in Guassian_data_total_energy_Hartree_list]

    # normalize so the min value is 0
    Guassian_data_total_energy_kcal_per_mol_normalize_list = \
        [i - min(Guassian_data_total_energy_kcal_per_mol_list)
         for i in Guassian_data_total_energy_kcal_per_mol_list]

    print(
        f"Guassian_data_dihedral_degrees_list = {Guassian_data_dihedral_degrees_list}")
    print(
        f"Guassian_data_total_energy_kcal_per_mol_list = {Guassian_data_total_energy_kcal_per_mol_list}")
    print(f"Guassian_data_total_energy_kcal_per_mol_normalize_list = "
          f"{Guassian_data_total_energy_kcal_per_mol_normalize_list}"
          )

    # get the Gaussian minus GOMC total energy and then it normalized
    Gaussian_minus_GOMC_data_dihedral_degrees_list = GOMC_data_dihedral_degrees_list
    Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list = [
        Guassian_data_total_energy_kcal_per_mol_normalize_list[i]
        - GOMC_data_total_energy_kcal_per_mol_normalize_list[i]
        for i in range(0, len(Guassian_data_total_energy_kcal_per_mol_normalize_list))
    ]

    Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list = [
        Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list[i] -
        min(Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list)
        for i in range(0, len(Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list))
    ]

    print(
        f'Gaussian_minus_GOMC_data_dihedral_degrees_list = {Gaussian_minus_GOMC_data_dihedral_degrees_list}')
    print(
        f'Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list = \
    {Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list}'
    )
    # *********************************
    # get Gaussian data (END)
    # *********************************

    # *********************************
    # get all other dihedral angles (phi) that match the atom type of the scanned dihedral,
    # and all sum 1/2*(k1 * scalar_i) = sum 1/2*(k1 * (1 +/- cos(n * phi)) values
    # (START)
    # *********************************

    [
        matching_dihedral_types_by_atom_numbers_list,
        matching_dihedral_types_by_atom_type_list,
        all_matching_dihedral_coordinates_angstroms_added_to_k_values_list,
        all_matching_dihedral_phi_degrees_added_to_k_values_list,
        all_sum_opls_const_1_plus_or_minus_cos_n_list
    ] = mdf_frw.get_matching_dihedral_info_and_opls_fitting_data(
        fit_dihedral_atom_types,
        f"{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}.psf",
        qm_log_files_and_entries_to_remove,
        qm_engine=qm_engine,
    )

    # get the individual sum_opls_const_1_plus_or_minus_cos_n_list ones for fitting and plotting
    const_1_minus_Cos_0_phi_data_lists = []
    const_1_plus_Cos_1_phi_data_lists = []
    const_1_minus_Cos_2_phi_data_lists = []
    const_1_plus_Cos_3_phi_data_lists = []
    const_1_minus_Cos_4_phi_data_lists = []
    for const_1_plus_or_minus_cos_i in all_sum_opls_const_1_plus_or_minus_cos_n_list:
        const_1_minus_Cos_0_phi_data_lists.append(
            const_1_plus_or_minus_cos_i[0])
        const_1_plus_Cos_1_phi_data_lists.append(
            const_1_plus_or_minus_cos_i[1])
        const_1_minus_Cos_2_phi_data_lists.append(
            const_1_plus_or_minus_cos_i[2])
        const_1_plus_Cos_3_phi_data_lists.append(
            const_1_plus_or_minus_cos_i[3])
        const_1_minus_Cos_4_phi_data_lists.append(
            const_1_plus_or_minus_cos_i[4])

    # *********************************
    # get all other dihedral angles (phi) that match the atom type of the scanned dihedral,
    # and all sum 1/2*(k1 * scalar_i) = sum 1/2*(k1 * (1 +/- cos(n * phi)) values
    # (END)
    # *********************************

    # Check if all the columns are the same length for GOMC and Gaussian data
    if not len(GOMC_data_dihedral_degrees_list) \
            == len(GOMC_data_total_energy_kcal_per_mol_list) \
            == len(GOMC_data_total_energy_kcal_per_mol_normalize_list) \
            == len(Guassian_data_dihedral_degrees_list) \
            == len(Guassian_data_total_energy_kcal_per_mol_list) \
            == len(Guassian_data_total_energy_kcal_per_mol_normalize_list) \
            == len(Gaussian_minus_GOMC_data_dihedral_degrees_list) \
            == len(Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_list) \
            == len(Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list) \
            == len(all_sum_opls_const_1_plus_or_minus_cos_n_list):
        raise ValueError(
            "ERROR: The GOMC and Guassian outputs do not match in length. "
            "This could mean something is changed and wrong in the code, "
            "or GOMC is outputting multiple Initial eneries in the log file "
            ", in this case use a new version of GOMC."
        )

    # Check if all the angles match between sorted GOMC and Gaussian data
    for j_angle in range(0, len(GOMC_data_dihedral_degrees_list)):
        if not len(GOMC_data_dihedral_degrees_list) == len(Guassian_data_dihedral_degrees_list):
            raise ValueError(
                "ERROR: The GOMC and Guassian output angles are not in the same angles in order.")

    # Check if all the angles match between sorted GOMC and Gaussian data
    for k_angle in range(0, len(GOMC_data_dihedral_degrees_list)):
        if not len(GOMC_data_dihedral_degrees_list) == len(Guassian_data_dihedral_degrees_list):
            raise ValueError(
                "ERROR: The GOMC and Guassian output angles are not in the same angles in order.")
        if k_angle == 0:
            # write out the GOMC and Gaussian data in a file
            gomc_gaussian_kcal_mol_energy_data_txt_file = open(
                gomc_gaussian_kcal_mol_energy_filename, "w")
            gomc_gaussian_kcal_mol_energy_data_txt_file.write(
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

        gomc_gaussian_kcal_mol_energy_data_txt_file.write(
            f"{Gaussian_minus_GOMC_data_dihedral_degrees_list[k_angle]: <30} "
            f"{GOMC_data_total_energy_kcal_per_mol_normalize_list[k_angle]: <30} "
            f"{Guassian_data_total_energy_kcal_per_mol_normalize_list[k_angle]: <30} "
            f"{Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list[k_angle]: <30} "
            f"{const_1_minus_Cos_0_phi_data_lists[k_angle]: <30} "
            f"{const_1_plus_Cos_1_phi_data_lists[k_angle]: <30} "
            f"{const_1_minus_Cos_2_phi_data_lists[k_angle]: <30} "
            f"{const_1_plus_Cos_3_phi_data_lists[k_angle]: <30} "
            f"{const_1_minus_Cos_4_phi_data_lists[k_angle]: <30} "
            f" \n"
        )
    gomc_gaussian_kcal_mol_energy_data_txt_file.close()

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
    ax1.set_xlabel('phi (degrees)', size=axis_Label_font_size)
    ax1.set_ylabel('Dihedral Energy (kcal/mol)', size=axis_Label_font_size)

    # sort the by Dihedral_Degrees at the same time for
    # Dihedral_Degrees, GOMC_E_kcal_per_mol, Gaussian_E_kcal_per_mol, and Gaussian_minus_GOMC_E_kcal_per_mol
    sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list, \
        sorted_GOMC_data_total_energy_kcal_per_mol_normalize_list, \
        sorted_Guassian_data_total_energy_kcal_per_mol_normalize_list, \
        sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list, \
        sorted_all_sum_opls_const_1_plus_or_minus_cos_n_list, \
        sorted_const_1_minus_Cos_0_phi_data_lists, \
        sorted_const_1_plus_Cos_1_phi_data_lists, \
        sorted_const_1_minus_Cos_2_phi_data_lists, \
        sorted_const_1_plus_Cos_3_phi_data_lists, \
        sorted_const_1_minus_Cos_4_phi_data_lists, \
        = zip(
            *sorted(zip(
                Gaussian_minus_GOMC_data_dihedral_degrees_list,
                GOMC_data_total_energy_kcal_per_mol_normalize_list,
                Guassian_data_total_energy_kcal_per_mol_normalize_list,
                Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list,
                all_sum_opls_const_1_plus_or_minus_cos_n_list,
                const_1_minus_Cos_0_phi_data_lists,
                const_1_plus_Cos_1_phi_data_lists,
                const_1_minus_Cos_2_phi_data_lists,
                const_1_plus_Cos_3_phi_data_lists,
                const_1_minus_Cos_4_phi_data_lists
            ))
        )

    plot_max = int(max(
        sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list) + 1.51)
    ax1.set_ylim(-int(plot_max * 1.5), plot_max)
    plt.title(
        "OPLS: All viable summed dihedral fits $\Longrightarrow$ $\Sigma$ matching dihedrals. \n"
        "(Non-zero k's = 1 and 3 -> label k_non_0='1_3')"
    )

    # loop thru dihderal_k_zeros_list_k0_k1_k2_k3_k4 list and fit all that are listed in here
    # add add the label designnator k's used (i.e., Non-zero k's = 1 and 3 -> label '1_3')
    # and write the constants out
    end_part_dihedral_k_constants_fit_energy_figure_filename = \
        "all_summed_dihedrals_k_constants_figure.pdf"

    all_individual_fit_dihedral_k_constants_figure_filename = \
        "all_single_fit_dihedral_k_constants_figure.pdf"

    end_part_dihedral_k_constants_fit_energy_filename = \
        "k_constants_fit_energy.txt"
    opls_dihedral_k_constants_fit_energy_kcal_mol_txt_file = open(
        f"{'opls_dihedral'}_{end_part_dihedral_k_constants_fit_energy_filename}", "w"
    )
    opls_dihedral_k_constants_fit_energy_kcal_mol_txt_file.write(
        f"{'non_zero_k_constants': <25} "
        f"{'k0_kcal_per_mol': <25} "
        f"{'k1_kcal_per_mol': <25} "
        f"{'k2_kcal_per_mol': <25} "
        f"{'k3_kcal_per_mol': <25} "
        f"{'k4_kcal_per_mol': <25} "
        f"{'r_squared': <25} "
        f" \n"
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
        '1',
        '2',
        '3',
        '4',
        '1_3',
        '2_4',
        '3_4',
        '1_2',
        '1_2_3',
        '1_2_3_4',
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
        25 / 100 * len(sorted_const_1_plus_Cos_1_phi_data_lists) + 1)
    for ck_const_i in range(0, len(sorted_const_1_plus_Cos_1_phi_data_lists)):
        # Check if can use cos power 1
        if bool(
                np.isclose(
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_plus_Cos_1_phi_data_lists[ck_const_i], sig_figs=check_sig_figs
                    ),
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_plus_Cos_1_phi_data_lists[0], sig_figs=check_sig_figs
                    )
                )
        ) is False:
            const_cos_count_powers_1_2_3_and_4_int_list[0] = 1 + \
                const_cos_count_powers_1_2_3_and_4_int_list[0]

        # Check if can use cos power 2
        if bool(
                np.isclose(
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_minus_Cos_2_phi_data_lists[ck_const_i], sig_figs=check_sig_figs
                    ),
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_minus_Cos_2_phi_data_lists[0], sig_figs=check_sig_figs
                    )
                )
        ) is False:
            const_cos_count_powers_1_2_3_and_4_int_list[1] = 1 + \
                const_cos_count_powers_1_2_3_and_4_int_list[1]

        # Check if can use cos power 3
        if bool(
                np.isclose(
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_plus_Cos_3_phi_data_lists[ck_const_i], sig_figs=check_sig_figs
                    ),
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_plus_Cos_3_phi_data_lists[0], sig_figs=check_sig_figs
                    )
                )
        ) is False:
            const_cos_count_powers_1_2_3_and_4_int_list[2] = 1 + \
                const_cos_count_powers_1_2_3_and_4_int_list[2]

        # Check if can use cos power 4
        if bool(
                np.isclose(
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_minus_Cos_4_phi_data_lists[ck_const_i], sig_figs=check_sig_figs
                    ),
                    mdf_math.round_to_sig_figs(
                        sorted_const_1_minus_Cos_4_phi_data_lists[0], sig_figs=check_sig_figs
                    )
                )
        ) is False:
            const_cos_count_powers_1_2_3_and_4_int_list[3] = 1 + \
                const_cos_count_powers_1_2_3_and_4_int_list[3]

        # add the cos power value if enough differences in the constants are seen
        if const_cos_count_powers_1_2_3_and_4_int_list[0] >= number_of_differnt_const_cos_needed_int:
            allowed_fitting_powers_1_2_3_and_4.append('1')

        if const_cos_count_powers_1_2_3_and_4_int_list[1] >= number_of_differnt_const_cos_needed_int:
            allowed_fitting_powers_1_2_3_and_4.append('2')

        if const_cos_count_powers_1_2_3_and_4_int_list[2] >= number_of_differnt_const_cos_needed_int:
            allowed_fitting_powers_1_2_3_and_4.append('3')

        if const_cos_count_powers_1_2_3_and_4_int_list[3] >= number_of_differnt_const_cos_needed_int:
            allowed_fitting_powers_1_2_3_and_4.append('4')

    allowed_fitting_powers_1_2_3_and_4.sort()

    # This selects only fit acceptable power based on dihedral symmetry of the molecule/fragment.
    # Not all the powers can be used when fitting a dihedral, so this outputs only the
    # powers that will produce the correct solution (i.e., removes all powers that don't
    # yeild the correct solution).
    fit_k_list_allowed = []
    for mod_i, mod_val in enumerate(fit_k_list):
        power_allowed_iter = True
        for str_i in mod_val:

            if str_i not in ['-', '_'] and str_i not in allowed_fitting_powers_1_2_3_and_4:
                power_allowed_iter = False
        if power_allowed_iter is True:
            fit_k_list_allowed.append(mod_val)

    # **********************************
    # Select only the OPLS fit cos powers that produce a unique solution
    # without any symmetry issues
    # (END)
    # **********************************

    # Run the fitting for only the allowed power types
    for k_iter_i, k_type_i in enumerate(fit_k_list_allowed):
        # make the list of k_type_i for fitting in the data
        k_type_list_i = []
        for v in range(0, len(sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list)):
            k_type_list_i.append(k_type_i)

        if k_type_i == '1':

            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )

            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[2] = 0
            parameters[3] = 0
            parameters[4] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '2':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[1] = 0
            parameters[3] = 0
            parameters[4] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '3':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[1] = 0
            parameters[2] = 0
            parameters[4] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '4':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[1] = 0
            parameters[2] = 0
            parameters[3] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '1_3':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[2] = 0
            parameters[4] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '2_4':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[1] = 0
            parameters[3] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '1_2':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[3] = 0
            parameters[4] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '3_4':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[1] = 0
            parameters[2] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '1_2_3':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0
            parameters[4] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        elif k_type_i == '1_2_3_4':
            parameters, covariance = curve_fit(
                mdf_math.opls_dihedral,
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                np.asarray(
                    sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list)
            )
            # fix parameters to zero for the unused values because we want the list length the same,
            # and the unused ones are auto-set to 1.
            parameters[0] = 0

            # fit the OPLS dihedral
            fit_opls_dihedral = mdf_math.opls_dihedral(
                (
                    np.asarray(k_type_list_i),
                    np.asarray(
                        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list),
                    np.asarray(sorted_const_1_minus_Cos_0_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_1_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_2_phi_data_lists),
                    np.asarray(sorted_const_1_plus_Cos_3_phi_data_lists),
                    np.asarray(sorted_const_1_minus_Cos_4_phi_data_lists)
                ),
                parameters[0],
                parameters[1],
                parameters[2],
                parameters[3],
                parameters[4],
            )

        else:
            raise ValueError(
                f"ERROR: The {k_type_i} selected in the 'fit_k_list' variable is not a valid selection")

        # calulate TSS, RSS and R**2
        r_squared = mdf_math.get_r_squared(
            sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list,
            fit_opls_dihedral
        )

        plt.plot(
            np.asarray(sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list), np.asarray(
                fit_opls_dihedral),
            '-', label=f" {k_type_i} | $R^{2}$={np.round(r_squared, decimals=4)}"
        )

        # wrie out the k constants and R^2
        opls_dihedral_k_constants_fit_energy_kcal_mol_txt_file.write(
            f"{k_type_i: <25} "
            f"{parameters[0]: <25} "
            f"{parameters[1]: <25} "
            f"{parameters[2]: <25} "
            f"{parameters[3]: <25} "
            f"{parameters[4]: <25} "
            f"{r_squared: <25} "
            f" \n"
        )

    # plot the data point that it is being fit too
    plt.plot(
        sorted_Gaussian_minus_GOMC_data_dihedral_degrees_list,
        sorted_Gaussian_minus_GOMC_data_total_energy_kcal_per_mol_normalized_list,
        'o', label='Data'
    )

    # close the file
    opls_dihedral_k_constants_fit_energy_kcal_mol_txt_file.close()

    major_xticks = np.arange(-180, 180 + 0.0001, 60)
    minor_xticks = np.arange(-180, 180 + 0.0001, 10)
    ax1.set_xticks(major_xticks)
    ax1.set_xticks(minor_xticks, minor=True)
    ax1.set_xlim(-180, 180 + 0.001)

    plt.legend(
        ncol=2,
        loc="lower center",
        fontsize=legend_font_size,
        prop={'family': 'Arial', 'size': legend_font_size}
    )

    # plt.show()
    fig1.savefig(
        f"opls_{end_part_dihedral_k_constants_fit_energy_figure_filename}", dpi=300)

    # *********************************
    # fit the Gaussian - GOMC dihedral (END)
    # *********************************

    # *********************************
    # Plot all OPLS dihedral to fitted forms together (START)
    # *********************************
    opls_fit_data_df = pd.DataFrame(
        pd.read_csv(f"{'opls_dihedral'}_{end_part_dihedral_k_constants_fit_energy_filename}",
                    sep='\s+',
                    header=0
                    )
    )
    opls_fit_data_non_zero_k_constants_list = opls_fit_data_df.loc[:, 'non_zero_k_constants'].tolist(
    )
    opls_fit_data_k0_kcal_per_mol_list = opls_fit_data_df.loc[:, 'k0_kcal_per_mol'].tolist(
    )
    opls_fit_data_k1_kcal_per_mol_list = opls_fit_data_df.loc[:, 'k1_kcal_per_mol'].tolist(
    )
    opls_fit_data_k2_kcal_per_mol_list = opls_fit_data_df.loc[:, 'k2_kcal_per_mol'].tolist(
    )
    opls_fit_data_k3_kcal_per_mol_list = opls_fit_data_df.loc[:, 'k3_kcal_per_mol'].tolist(
    )
    opls_fit_data_k4_kcal_per_mol_list = opls_fit_data_df.loc[:, 'k4_kcal_per_mol'].tolist(
    )
    opls_fit_data_r_squared_list = opls_fit_data_df.loc[:, 'r_squared'].tolist(
    )

    # *********************************
    # Plot all OPLS dihedral to fitted forms together  (END)
    # *********************************

    # *********************************
    # Convert the OPLS dihedral to other forms (START)
    # *********************************

    # create the Periodic / CHARMM dihedrals file
    periodic_dihedral_k_constants_fit_energy_kcal_mol_txt_file = open(
        f"{'periodic_dihedral'}_{end_part_dihedral_k_constants_fit_energy_filename}", "w"
    )
    periodic_dihedral_k_constants_fit_energy_kcal_mol_txt_file.write(
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
        f" \n"
    )

    # create the RB torsions file
    RB_torsion_k_constants_fit_energy_kcal_mol_txt_file = open(
        f"{'RB_torsion'}_{end_part_dihedral_k_constants_fit_energy_filename}", "w"
    )
    RB_torsion_k_constants_fit_energy_kcal_mol_txt_file.write(
        f"{'non_zero_k_constants': <25} "
        f"{'k0_kcal_per_mol': <25} "
        f"{'k1_kcal_per_mol': <25} "
        f"{'k2_kcal_per_mol': <25} "
        f"{'k3_kcal_per_mol': <25} "
        f"{'k4_kcal_per_mol': <25} "
        f"{'k5_kcal_per_mol': <25} "
        f"{'r_squared': <25} "
        f" \n"
    )

    # loop thru the different 'non_zero_k_constants' for the OPLS dihedral
    for opls_fit_i in range(0, len(opls_fit_data_non_zero_k_constants_list)):
        # list the phi values to check
        phi_check_degrees = 1
        phi_check_number_of_degree_values = int(360 / phi_check_degrees)
        phi_values_for_check_degrees_list = \
            [-180 + i *
                phi_check_degrees for i in range(0, phi_check_number_of_degree_values)]

        cos_power_list = [
            opls_fit_data_non_zero_k_constants_list[opls_fit_i]
            for i in range(0, len(phi_values_for_check_degrees_list))
        ]

        # check if the periodic and opls dihedral energies match all the values in the list
        for period_opls_ck_i in range(0, len(phi_values_for_check_degrees_list)):
            # calculate the dihedral energy for a check to the converted dihedrals
            opls_dihedral_energy_check = mdf_math.opls_dihedral(
                (
                    cos_power_list[opls_fit_i],
                    phi_values_for_check_degrees_list[opls_fit_i],
                    None,
                    None,
                    None,
                    None,
                    None
                ),
                opls_fit_data_k0_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k1_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k2_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k3_kcal_per_mol_list[opls_fit_i],
                opls_fit_data_k4_kcal_per_mol_list[opls_fit_i]
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
                opls_fit_data_k4_kcal_per_mol_list[opls_fit_i]
            )

            periodic_dihedral_energy_check = mdf_math.periodic_dihedral_n_1_2_3_4_5(
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

            # check periodic dihedral equation fit matches the OPLS value
            if not math.isclose(opls_dihedral_energy_check, periodic_dihedral_energy_check):
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
                opls_fit_data_k4_kcal_per_mol_list[opls_fit_i]
            )

            RB_torsion_energy_check = mdf_math.RB_torsion_n_1_2_3_4_5(
                phi_values_for_check_degrees_list[opls_fit_i],
                RB_torsion_k_values[0],
                RB_torsion_k_values[1],
                RB_torsion_k_values[2],
                RB_torsion_k_values[3],
                RB_torsion_k_values[4],
                RB_torsion_k_values[5]

            )

            # check periodic dihedral equation fit matches the OPLS value
            if not math.isclose(opls_dihedral_energy_check, RB_torsion_energy_check):
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
        periodic_dihedral_k_constants_fit_energy_kcal_mol_txt_file.write(
            f"{f'0_{opls_fit_data_non_zero_k_constants_list[opls_fit_i]}': <25} "
            f"{periodic_dihedral_k_n_d_values[0][0]: <25} "
            f"{periodic_dihedral_k_n_d_values[1][0]: <25} "
            f"{periodic_dihedral_k_n_d_values[2][0]: <25} "
            f"{periodic_dihedral_k_n_d_values[3][0]: <25} "
            f"{periodic_dihedral_k_n_d_values[4][0]: <25} "
            f"{periodic_dihedral_k_n_d_values[5][0]: <25} "
            f"{periodic_dihedral_k_n_d_values[0][1]: <25} "
            f"{periodic_dihedral_k_n_d_values[1][1]: <25} "
            f"{periodic_dihedral_k_n_d_values[2][1]: <25} "
            f"{periodic_dihedral_k_n_d_values[3][1]: <25} "
            f"{periodic_dihedral_k_n_d_values[4][1]: <25} "
            f"{periodic_dihedral_k_n_d_values[5][1]: <25} "
            f"{periodic_dihedral_k_n_d_values[0][2]: <25} "
            f"{periodic_dihedral_k_n_d_values[1][2]: <25} "
            f"{periodic_dihedral_k_n_d_values[2][2]: <25} "
            f"{periodic_dihedral_k_n_d_values[3][2]: <25} "
            f"{periodic_dihedral_k_n_d_values[4][2]: <25} "
            f"{periodic_dihedral_k_n_d_values[5][2]: <25} "
            f"{opls_fit_data_r_squared_list[opls_fit_i]: <25} "
            f" \n"
        )

        # Write the RB torsion constants to file.
        # Reused the OPLS R^2 fitted values as the energeries where
        # validated between energies RB torsion style and OPLS.
        # Added the '0_' to the 'non_zero_k_constants' as the k0 needed
        # for the OPLS conversion to the RB torsion style.
        RB_torsion_k_constants_fit_energy_kcal_mol_txt_file.write(
            f"{f'0_{opls_fit_data_non_zero_k_constants_list[opls_fit_i]}': <25} "
            f"{RB_torsion_k_values[0]: <25} "
            f"{RB_torsion_k_values[1]: <25} "
            f"{RB_torsion_k_values[2]: <25} "
            f"{RB_torsion_k_values[3]: <25} "
            f"{RB_torsion_k_values[4]: <25} "
            f"{RB_torsion_k_values[5]: <25} "
            f"{opls_fit_data_r_squared_list[opls_fit_i]: <25} "
            f" \n"
        )

        # *********************************
        # Write out the periodic and RB torsions values to file (transformed from OPLS)  (END)
        # NOTE: For RB torsions, the angle psi is used, which is psi = phi - Pi
        # *********************************

    # close the Periodic / CHARMM dihedrals file
    periodic_dihedral_k_constants_fit_energy_kcal_mol_txt_file.close()

    # close the RB torsions file
    RB_torsion_k_constants_fit_energy_kcal_mol_txt_file.close()

    # *********************************
    # Convert the OPLS dihedral to other forms (END)
    # *********************************

    # *********************************
    # Check the all the OPLS dihedral forms are correct
    # by running GOMC with the fitted values and comparing it to QM
    # (START)
    # *********************************
    opls_k_constant_fitted_q_list = []
    opls_r_squared_fitted_data_via_gomc_list = []
    for opls_q, opls_fit_q in enumerate(opls_fit_data_non_zero_k_constants_list):

        gomc_fitted_gaussian_kcal_mol_energy_filename = \
            f"all_normalized_energies_OPLS_fit_{opls_fit_q}_in_kcal_per_mol.txt"
        # write the GOMC control files
        for scan_iter_q in range(1, len(Guassian_raw_degrees_list) + 1):
            read_gomc_fitted_restart_file_coor_dir_and_name = \
                f'../{xyz_xsc_coor_files_directory}/dihedral_coords_position_{scan_iter_q}.coor'
            read_gomc_fitted_restart_file_xsc_dir_and_name = f'../{xyz_xsc_coor_files_directory}/starting_point.xsc'

            # The gomc raw energy filename in Kelvin Energies
            gomc_raw_energy_fitted_filename = f"gomc_raw_OPLS_fit_{opls_fit_q}_energies_in_Kelvin.txt"

            # The combined GOMC and Gaussian dihedral angles and energies in kcal/mol
            gomc_gaussian_kcal_mol_energy_fitted_filename = \
                f"all_normalized_OPLS_fit_{opls_fit_q}_energies_in_kcal_per_mol.txt"

            control_file_name_fitted_str = f'GOMC_OPLS_fit_{opls_fit_q}_dihedral_coords_{scan_iter_q}.conf'
            output_name_control_fitted_file_name_str = \
                f'output_GOMC_OPLS_fit_{opls_fit_q}_dihedral_coords_{scan_iter_q}.txt'

            print("#**********************")
            print(
                "Started: Writing NVT GOMC control file for the GOMC simulation with fitted dihedral k values.")
            print("#**********************")
            gomc_control.write_gomc_control_file(
                charmm,
                f'{gomc_runs_folder_name}/{control_file_name_fitted_str}',
                'NVT',
                MC_steps,
                temperature,
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
                    "VDWGeometricSigma": override_VDWGeometricSigma,
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
                    "Potential": "VDW",
                    "CBMC_First": 12,
                    "CBMC_Nth": 10,
                    "CBMC_Ang": 50,
                    "CBMC_Dih": 50,
                },
            )
            print("#**********************")
            print("Completed: NVT GOMC control file written for the GOMC simulation with fitted dihedral k values.")
            print("#**********************")

            # **************************************************************
            # make the GOMC control file for each dihedral angle (END)
            # **************************************************************

            # add the modified k values for the simulation
            opls_k_constant_fitted_q_list = [
                opls_fit_data_k0_kcal_per_mol_list[opls_q],
                opls_fit_data_k1_kcal_per_mol_list[opls_q],
                opls_fit_data_k2_kcal_per_mol_list[opls_q],
                opls_fit_data_k3_kcal_per_mol_list[opls_q],
                opls_fit_data_k4_kcal_per_mol_list[opls_q],
            ]

            mdf_frw.change_gomc_ff_file_dihedral_values(
                f'{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_dihedrals_per_xml.inp',
                f'{gomc_runs_folder_name}/{output_gomc_pdb_psf_ff_file_name_str}_OPLS_fit_{opls_fit_q}_dihedral.inp',
                fit_dihedral_atom_types,
                fit_dihedral_opls_k_0_1_2_3_4_values=opls_k_constant_fitted_q_list,
                zeroed_dihedral_atom_types=None,
            )

            # **************************************************************
            # Run GOMC and get the system/fragment energy (START)
            # **************************************************************

            # only run NVT as we only want/need the initial energy and want the box a constant size
            run_gomc_command = \
                f"cd {gomc_runs_folder_name} && {gomc_binary_path}/GOMC_CPU_NVT +p{gomc_cpu_cores} " \
                f"{control_file_name_fitted_str} > {output_name_control_fitted_file_name_str}"

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
                f'{gomc_runs_folder_name}/{output_name_control_fitted_file_name_str}', "r"
            ).readlines()

            get_e_titles = True
            for log_file_iter, log_file_line_iter in enumerate(read_gomc_log_fitted_file):
                log_file_splitline_iter = log_file_line_iter.split()

                # scan_iter starts at 1
                dihedral_angle_degrees = Guassian_raw_degrees_list[scan_iter_q - 1]

                # only open the gomc raw energy file and write header for 1st iteration (1)
                if len(log_file_splitline_iter) >= 2:
                    if scan_iter_q == 1 and log_file_splitline_iter[0] == "ETITLE:" and get_e_titles is True:
                        gomc_combined_raw_fitted_energy = open(
                            gomc_raw_energy_fitted_filename, "w")
                        # remove the wrongly entered 5 spaces in before "ETITLE:
                        # (This will be fixed in GOMC so it is not required)
                        extra_spaces_for_header_space_gomc_bug = 5
                        gomc_combined_raw_fitted_energy.write(f'{"Dihedral_Position": <19} '
                                                              f'{"Dihedral_Degrees": <19} '
                                                              f'{log_file_line_iter[extra_spaces_for_header_space_gomc_bug:]}'
                                                              )
                        if get_e_titles is True:
                            get_e_titles = False

                    # get only the initial configuration energy lin (i.e., step 0 line)
                    if log_file_splitline_iter[0] == 'ENER_0:' and log_file_splitline_iter[1] == '0':
                        gomc_combined_raw_fitted_energy.write(f'{scan_iter_q: <19}'
                                                              f'{dihedral_angle_degrees: <19} '
                                                              f'{log_file_line_iter[5:]}'
                                                              )

        # This fails when the GOMC simulations did not run
        try:
            # close the gomc_combined_raw_energy file
            gomc_combined_raw_fitted_energy.close()

        except:
            raise ValueError("ERROR: The GOMC fit test simulations did not run. "
                             "There is likely an error in creating the "
                             "required GOMC test simulations files or user inputs to the desired files.")

        # extract the raw data
        GOMC_data_fitted_df = pd.DataFrame(pd.read_csv(
            gomc_raw_energy_fitted_filename, sep='\s+'))
        GOMC_data_fitted_dihedral_degrees_list = GOMC_data_fitted_df.loc[:, 'Dihedral_Degrees'].tolist(
        )
        GOMC_data_fitted_total_energy_K_list = GOMC_data_fitted_df.loc[:, 'TOTAL'].tolist(
        )
        print(
            f'GOMC_data_fitted_dihedral_degrees_list = {GOMC_data_fitted_dihedral_degrees_list}')
        print(
            f'GOMC_data_fitted_total_energy_K_list = {GOMC_data_fitted_total_energy_K_list}')

        # convert from Kelvin to kcal/mol normalize so the min value is 0
        GOMC_data_fitted_total_energy_kcal_per_mol_list = \
            [i * conversion_K_to_kcal_per_mol for i in GOMC_data_fitted_total_energy_K_list]
        GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list = \
            [i - min(GOMC_data_fitted_total_energy_kcal_per_mol_list)
             for i in GOMC_data_fitted_total_energy_kcal_per_mol_list]

        print(
            f'GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list = {GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list}')

        # *********************************
        # get GOMC data (END)
        # *********************************

        # *********************************
        # get Gaussian data (START)
        # *********************************
        # extract the raw data
        Guassian_data_df = pd.DataFrame(pd.read_csv(
            qm_energy_file_dir_and_name, sep='\s+', header=3))
        Guassian_data_fitted_dihedral_degrees_list = Guassian_data_df.iloc[:, 0].tolist(
        )
        Guassian_data_total_energy_Hartree_list = Guassian_data_df.iloc[:, 1].tolist(
        )

        # convert from Hartree to kcal/mol energy units
        Guassian_data_total_energy_kcal_per_mol_list = \
            [i * conversion_hartree_to_kcal_per_mol for i in Guassian_data_total_energy_Hartree_list]

        # normalize so the min value is 0
        Guassian_data_total_energy_kcal_per_mol_normalize_list = \
            [i - min(Guassian_data_total_energy_kcal_per_mol_list) for i in
             Guassian_data_total_energy_kcal_per_mol_list]

        # get the Gaussian minus GOMC total energy and then it normalized
        Gaussian_minus_GOMC_data_fitted_dihedral_degrees_list = GOMC_data_fitted_dihedral_degrees_list
        Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_list = [
            Guassian_data_total_energy_kcal_per_mol_normalize_list[i]
            - GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list[i]
            for i in range(0, len(Guassian_data_total_energy_kcal_per_mol_normalize_list))
        ]

        Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_normalized_list = [
            Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_list[i] -
            min(Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_list)
            for i in range(0, len(Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_list))
        ]

        # get R**2 for the fit, running through GOMC to get the new energy of the
        # individual fit.
        opls_r_squared_fitted_data_via_gomc_iter = mdf_math.get_r_squared(
            Guassian_data_total_energy_kcal_per_mol_normalize_list,
            GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list
        )
        opls_r_squared_fitted_data_via_gomc_list.append(
            opls_r_squared_fitted_data_via_gomc_iter)

        # Use R**2 (R-squared) to compare the data used to fit the dihedral(s)
        # via a single or multi-dihedral fit, to the individual fit entered in
        # GOMC and then recompared to Gaussian and write out the data to a file.
        for q_angle in range(0, len(GOMC_data_fitted_dihedral_degrees_list)):
            if q_angle == 0:
                # write out the GOMC and Gaussian data in a file
                gomc_fitted_gaussian_kcal_mol_energy_data_txt_file = \
                    open(gomc_fitted_gaussian_kcal_mol_energy_filename, "w")
                gomc_fitted_gaussian_kcal_mol_energy_data_txt_file.write(
                    f"{'Dihedral_Degrees': <30} "
                    f"{'GOMC_E_kcal_per_mol': <30} "
                    f"{'Gaussian_E_kcal_per_mol': <30} "
                    f"{'Gaussian_minus_GOMC_E_kcal_per_mol': <40} "
                    f"{'k0_OPLS_kcal_mol': <30} "
                    f"{'k1_OPLS_kcal_mol': <30} "
                    f"{'k2_OPLS_kcal_mol': <30} "
                    f"{'k3_OPLS_kcal_mol': <30} "
                    f"{'k4_OPLS_kcal_mol': <30} "
                    f" \n"
                )

            gomc_fitted_gaussian_kcal_mol_energy_data_txt_file.write(
                f"{Gaussian_minus_GOMC_data_fitted_dihedral_degrees_list[q_angle]: <30} "
                f"{GOMC_data_fitted_total_energy_kcal_per_mol_normalize_list[q_angle]: <30} "
                f"{Guassian_data_total_energy_kcal_per_mol_normalize_list[q_angle]: <30} "
                f"{Gaussian_minus_GOMC_data_fitted_total_energy_kcal_per_mol_normalized_list[q_angle]: <30} "
                f"{str(opls_k_constant_fitted_q_list[0]): <30} "
                f"{str(opls_k_constant_fitted_q_list[1]): <30} "
                f"{str(opls_k_constant_fitted_q_list[2]): <30} "
                f"{str(opls_k_constant_fitted_q_list[3]): <30} "
                f"{str(opls_k_constant_fitted_q_list[4]): <30} "
                f" \n"
            )

        # Compare original fit vs run through GOMC as a validation test case
        if opls_fit_data_r_squared_list[opls_q] >= fit_min_validated_r_squared \
                and not np.isclose(
                opls_r_squared_fitted_data_via_gomc_list[opls_q],
                opls_fit_data_r_squared_list[opls_q],
                rtol=fit_validation_r_squared_rtol
        ):
            raise ValueError(
                f"ERROR: The calculated R-squared energy values from the fit type "
                f"{opls_fit_data_non_zero_k_constants_list[opls_q]} does not match "
                f"the validated case for "
                f"'fit_min_validated_r_squared' >= {fit_min_validated_r_squared}, "
                f"within the relative tolerance or "
                f"'fit_validation_r_squared_rtol' = {fit_validation_r_squared_rtol}. \n"
                f"- Fit via the individual or multi-dihedral fit, when "
                f"Gaussian minus GOMC with the selected dihedral set to zero \n"
                f"--> R-squared = {opls_fit_data_r_squared_list[opls_q]} \n"
                f"- Fit via the validation test case, when "
                f"Gaussian minus GOMC with the selected individual dihedral added in GOMC \n"
                f"-- >R-squared = {opls_r_squared_fitted_data_via_gomc_list[opls_q]} \n"
                f"The 'fit_min_validated_r_squared' and 'fit_validation_r_squared_rtol' "
                f"variables may need to be adjusted, \n"
                f"there is likely something wrong with the fitting procedure, the "
                f"software parameters need tuned, or there is a bug in the software. \n\n "
                f"NOTE: Since the R-squared values are calculated via different parameters, \n"
                f"the compared R-squared values could be very different if they are not nearly \n"
                f"a perfect fit (R-squared --> 0.98 to 1)."
                f""
            )

    gomc_fitted_gaussian_kcal_mol_energy_data_txt_file.close()

    # *********************************
    # Check the all the OPLS dihedral forms are correct
    # by running GOMC with the fitted values and comparing it to QM
    # (END)
    # *********************************
    fig2, (ax2) = plt.subplots(1)
    axis_Label_font_size = 12
    legend_font_size = 12
    ax2.set_xlabel('phi (degrees)', size=axis_Label_font_size)
    ax2.set_ylabel('Dihedral Energy (kcal/mol)', size=axis_Label_font_size)

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
                opls_fit_data_non_zero_k_constants_list[opls_fit_j])

            opls_dihedral_energy_iter = mdf_math.opls_dihedral(
                (
                    cos_power_list_iter[-1],
                    phi_degrees_list_iter[-1],
                    None,
                    None,
                    None,
                    None,
                    None
                ),
                opls_fit_data_k0_kcal_per_mol_list[opls_fit_j],
                opls_fit_data_k1_kcal_per_mol_list[opls_fit_j],
                opls_fit_data_k2_kcal_per_mol_list[opls_fit_j],
                opls_fit_data_k3_kcal_per_mol_list[opls_fit_j],
                opls_fit_data_k4_kcal_per_mol_list[opls_fit_j]
            )

            opls_dihedral_energy_list_iter.append(opls_dihedral_energy_iter)

        max_energy_list_iter.append(max(opls_dihedral_energy_list_iter))

        plt.plot(
            np.asarray(phi_degrees_list_iter),
            np.asarray(opls_dihedral_energy_list_iter),
            '-',
            label=f" {opls_fit_data_non_zero_k_constants_list[opls_fit_j]} "
                  f"| $R^{2}$={np.round(opls_r_squared_fitted_data_via_gomc_list[opls_fit_j], decimals=4)}"
        )

    ax2.set_xticks(major_xticks)
    ax2.set_xticks(minor_xticks, minor=True)
    ax2.set_xlim(-180, 180 + 0.001)

    plt.legend(
        ncol=2,
        loc="lower center",
        fontsize=legend_font_size,
        prop={'family': 'Arial', 'size': legend_font_size}
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
        dpi=300
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
