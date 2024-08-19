import os
import shutil

import numpy as np
import vmd
from mosdef_gomc.utils.conversion import OPLS_to_periodic
from scipy.optimize import curve_fit

import mosdef_dihedral_fit.utils.math_operations as mdf_math


def get_atom_names_and_elements_from_mol2(mol2_directory_and_filename):
    """Get the atom names and element names from the mol2 file.

    This gets the atom names and elements from the mol2 file,
    outputting them as lists.

    Parameters
    ----------
    mol2_directory_and_filename: str
        The mol2 path/directory and the filename

    Returns
    -------
    atom_name_list: list, Example = [atom_name_0, atom_name_1, ..., atom_name_n]
        A list of atom names in the mol2 file, in order.
    element_name_list: list, Example = [element_name_0, element_name_1, ..., element_name_n]
        A list of element names in the mol2 file, in order.
    """
    atom_name_list = []
    element_name_list = []
    get_atom_type_bool = False
    with open(mol2_directory_and_filename, "r") as fp:
        readlines_mol2_file = fp.readlines()
        for m, line_m in enumerate(readlines_mol2_file):
            split_line_m = line_m.split()

            if len(split_line_m) > 0 and str(split_line_m[0]) in [
                "@<TRIPOS>BOND"
            ]:
                get_atom_type_bool = False

            if len(split_line_m) in [8, 9] and get_atom_type_bool is True:
                atom_name_list.append(split_line_m[1])
                element_name_list.append(split_line_m[5])

            if len(split_line_m) > 0 and str(split_line_m[0]) in [
                "@<TRIPOS>ATOM"
            ]:
                get_atom_type_bool = True

    if len(atom_name_list) == 0 or len(element_name_list) == 0:
        error_message = (
            "ERROR: The provided mol2 format is not the required VMD TRIPOS format "
            "or the mol2 file has zero (0) atoms or beads in it."
        )
        raise TypeError(error_message)

    return [atom_name_list, element_name_list]


def get_atom_names_and_elements_from_pdb(pdb_directory_and_filename):
    """Get the atom names and element names from the pdb file.

    This gets the atom names and elements from the pdb file,
    outputting them as lists.

    NOTE: This only works if the first column for all atoms/beads is
    labeled 'ATOM' or 'HETATM', which is currently the only output for the PDB
    files via MoSDeF-GOMC.

    Parameters
    ----------
    pdb_directory_and_filename: str
        The PDB path/directory and the filename

    Returns
    -------
    atom_name_list: list, Example = [atom_name_0, atom_name_1, ..., atom_name_n]
        A list of atom names in the PDB file, in order.
    element_name_list: list, Example = [element_name_0, element_name_1, ..., element_name_n]
        A list of element names in the PDB file, in order.
    """
    atom_name_list = []
    element_name_list = []
    get_atom_type_bool = False
    with open(pdb_directory_and_filename, "r") as fp:
        readlines_pdb_file = fp.readlines()
        for m, line_m in enumerate(readlines_pdb_file):
            split_line_m = line_m.split()

            if len(split_line_m) > 0 and str(split_line_m[0]) in ["END"]:
                get_atom_type_bool = False

            # and len(split_line_m) in [10, 11, 12]
            if (
                split_line_m[0] in ["ATOM", "HETATM"]
                and get_atom_type_bool is True
            ):
                atom_name_j = ""
                element_name_j = ""
                for fix_space_j in range(12, 12 + len(line_m[12:16])):
                    atom_name_str_space_iter = line_m[fix_space_j]
                    if atom_name_str_space_iter != " ":
                        atom_name_j += atom_name_str_space_iter

                for fix_space_k in range(76, 76 + len(line_m[76:79])):
                    element_name_str_space_iter = line_m[fix_space_k]
                    if element_name_str_space_iter != " ":
                        element_name_j += element_name_str_space_iter

                atom_name_list.append(atom_name_j)
                element_name_list.append(element_name_j)

            if len(split_line_m) > 0 and str(split_line_m[0]) in ["CRYST1"]:
                get_atom_type_bool = True

    if len(atom_name_list) == 0:
        raise TypeError(
            "ERROR: The provided pdb format is not the required PDB format "
            "or the pdb file has zero (0) atoms or beads in it."
        )

    return [atom_name_list, element_name_list]


def write_xyz_file_from_gaussian_coordinates(
    atom_names_list,
    qm_parital_coordinate_file_starting_dir_and_name,
    qm_coordinate_file_extension,
    xyz_files_directory,
    total_qm_scans,
):
    """Write the xyz files from the extracted Gaussian coordinate file.

    This write the xyz files from the extracted Gaussian coordinate file.
    This Gaussian coordinate file maintains the GausView coordinate file format.

    Parameters
    ----------
    atom_name_list: list, Example = [atom_name_0, atom_name_1, ..., atom_name_n]
        A list of atom names in the PDB file, in order.
    qm_parital_coordinate_file_starting_dir_and_name: str
        The first part of the QM directory and file name. There are many files or
        dihedral angles to evaluate so each file will have an numeric number starting
        at 1 to n. In order to read these from 1 to n, this is the first part of the
        file name, where the full file name is listed as follows:
        f'{qm_parital_coordinate_file_starting_dir_and_name}{i_iter}.{qm_coordinate_file_extension}'
    qm_coordinate_file_extension: str
        The last part of the QM directory and file name. There are many files or
        dihedral angles to evaluate so each file will have a numeric number starting
        at 1 to n. In order to read these from 1 to n, this is the last part of the
        file name, where the full file name is listed as follows:
        f'{qm_parital_coordinate_file_starting_dir_and_name}{i_iter}.{qm_coordinate_file_extension}'
    xyz_files_directory: str
        The directory where the xyz files will be written, with the file name
        'dihedral_coords_position_n.xyz'
    total_qm_scans: int
        The number of QM scans in the QM coordinate file.

    Outputs
    -------
    Writes the xyz file in the file in the selected 'xyz_files_directory' directory,
    with the file name 'dihedral_coords_position_n.xyz'.
    """
    # check if file extension starts with a "." or not remove for correct insertion later
    if (
        len(qm_coordinate_file_extension) > 0
        and qm_coordinate_file_extension[0] == "."
    ):
        qm_coordinate_file_extension = qm_coordinate_file_extension[1:]

    elif len(qm_coordinate_file_extension) > 0:
        qm_coordinate_file_extension = qm_coordinate_file_extension

    else:
        raise ValueError(
            "ERROR: The 'qm_coordinate_file_extension' variable extension is not listed or a an empty string."
        )

    # add 1 spaces infront of atom_names_list
    atom_names_with_2_spaces_in_front = []
    for j_iter in range(0, len(atom_names_list)):
        atom_names_with_2_spaces_in_front.append(f"  {atom_names_list[j_iter]}")

    for i_iter in range(1, total_qm_scans + 1):
        read_gausian_file_dir_name = f"{qm_parital_coordinate_file_starting_dir_and_name}{i_iter}.{qm_coordinate_file_extension}"

        new_file_name_each_point = (
            f"{xyz_files_directory}/dihedral_coords_position_{i_iter}.xyz"
        )

        number_of_atoms = str(int(len(atom_names_with_2_spaces_in_front)))
        comment_1_space_in_front = (
            f" generated by Python script from Gaussian data"
        )

        output_file_xyz_file = open(new_file_name_each_point, "w")
        output_file_xyz_file.write(f"{number_of_atoms}\n")
        output_file_xyz_file.write(f"{comment_1_space_in_front}\n")

        # The Gaussian file need to have these columns "Row	Highlight	Display	Tag	Symbol	X	Y	Z"
        check_gaussian_optimized_coordinate_file_correct(
            read_gausian_file_dir_name
        )
        with open(read_gausian_file_dir_name, "r") as fp:
            readlines_gausian_file = fp.readlines()
            for i, line in enumerate(readlines_gausian_file):
                if i == 0 and not (
                    "Row	Highlight	Display	Tag	Symbol	X	Y	Z" in line
                ):
                    raise ValueError(
                        "# The Gaussian file need to have these columns "
                        "'Row	Highlight	Display	Tag	Symbol	X	Y	Z'"
                    )

                if i != 0:
                    output_file_xyz_file.write(
                        f"{atom_names_with_2_spaces_in_front[i-1]: <20}\t"
                        f"{readlines_gausian_file[i].split()[5]: <20}\t"
                        f"{readlines_gausian_file[i].split()[6]: <20}\t"
                        f"{readlines_gausian_file[i].split()[7]: <20}\n"
                    )

    output_file_xyz_file.close()


def write_restart_coor_from_xyz_file(coor_files_directory, total_qm_scans):
    """Write the restart coor file from the xyz file.

    This function utilized VMD to write the restart coor files from the xyz file.

    NOTE: NAMD or GOMC need to be restarted from the restart coor file to have
    more coordinate precision than the PDB file format provides.

    Parameters
    ----------
    coor_files_directory: str
        The directory where the restart coor files will be written, with the file name
        'dihedral_coords_position_n.coor'
    total_qm_scans: int
        The number of QM scans in the QM coordinate file.

    Outputs
    -------
    Writes and the restart coor file in the file in the selected
    'coor_files_directory' directory, with the file name
    'dihedral_coords_position_n.coor'.
    """
    # *******************************************
    # write the "write_restart_files.tcl" file that will be run using VMD
    # (START)
    # *******************************************
    write_coor_vmd_source_file = f"write_restart_files.tcl"
    output_file_restart_coor = open(
        f"{coor_files_directory}/{write_coor_vmd_source_file}", "w"
    )

    output_file_restart_coor.write(f"package require topotools\n\n")
    output_file_restart_coor.write(f"set NumberPositions {total_qm_scans}\n")
    output_file_restart_coor.write(
        f'set xyzAndCoorDirectoryBaseFilename "dihedral_coords_position_"\n'
    )
    output_file_restart_coor.write(f'set xyzExtension ".xyz"\n')
    output_file_restart_coor.write(f'set NAMDbinExtension ".coor"\n\n')

    output_file_restart_coor.write(
        "for {set x 1} {$x <= $NumberPositions} {incr x} {\n"
        "\tset xyzFilename $xyzAndCoorDirectoryBaseFilename$x$xyzExtension\n"
        "\ttopo readvarxyz $xyzFilename\n\n"
        "\tset NAMDBinFilename $xyzAndCoorDirectoryBaseFilename$x$NAMDbinExtension\n"
        "\tanimate write namdbin $NAMDBinFilename\n\n"
        "}\n\n"
    )

    output_file_restart_coor.close()
    # *******************************************
    # write the "write_restart_files.tcl" file that will be run using VMD
    # (END)
    # *******************************************

    # *******************************************
    # change to the 'coor_files_directory' directory and write the .coor restart files via VMD
    # (START)
    # *******************************************

    current_dir = os.getcwd()
    os.chdir(f"{coor_files_directory}")
    vmd.evaltcl(f"source {write_coor_vmd_source_file}")
    os.chdir(current_dir)  # go back to starting directory

    # *******************************************
    # change to the 'coor_files_directory' directory and write the .coor restart files via VMD
    # (END)
    # *******************************************


def check_gaussian_angle_energy_file_correct(gaussian_energy_file_dir_and_name):
    """Check that the gaussian/Gausview file containing the angle and energy is formatted correctly.

    This checks the gaussian/Gausview containing the angle and energy is formatted correctly.

    The proper header format for the GausView/Gaussian output is as follows:

    | # Scan of Total Energy
    | # X-Axis:  Scan Coordinate
    | # Y-Axis:  Total Energy (Hartree)
    | #                  X                   Y

    Parameters
    ----------
    gaussian_energy_file_dir_and_name: str
        The directory and filename of the Gaussian/Gausview angle and energy file.

    Returns
    -------
    gaussian_angle_energy_file_correct: True or TypeError
        True; if the file is formatted correctly
        TypeError; if the file is not formatted correctly
    """
    with open(gaussian_energy_file_dir_and_name, "r") as fp:
        gaussian_dihedral_header_line_correct_bool_list = [
            False,
            False,
            False,
            False,
            False,
        ]
        gaussian_dihedral_header_file = fp.readlines()
        for m, line_m in enumerate(gaussian_dihedral_header_file):
            split_line_m = line_m.split()

            if m == 0 and len(split_line_m) == 5:
                if (
                    split_line_m[0] == "#"
                    and split_line_m[1] == "Scan"
                    and split_line_m[2] == "of"
                    and split_line_m[3] == "Total"
                    and split_line_m[4] == "Energy"
                ):
                    gaussian_dihedral_header_line_correct_bool_list[0] = True

            if m == 1 and len(split_line_m) == 4:
                if (
                    split_line_m[0] == "#"
                    and split_line_m[1] == "X-Axis:"
                    and split_line_m[2] == "Scan"
                    and split_line_m[3] == "Coordinate"
                ):
                    gaussian_dihedral_header_line_correct_bool_list[1] = True

            if m == 2 and len(split_line_m) == 5:
                if (
                    split_line_m[0] == "#"
                    and split_line_m[1] == "Y-Axis:"
                    and split_line_m[2] == "Total"
                    and split_line_m[3] == "Energy"
                    and split_line_m[4] == "(Hartree)"
                ):
                    gaussian_dihedral_header_line_correct_bool_list[2] = True

            if m == 3 and len(split_line_m) == 3:
                if (
                    split_line_m[0] == "#"
                    and split_line_m[1] == "X"
                    and split_line_m[2] == "Y"
                ):
                    gaussian_dihedral_header_line_correct_bool_list[3] = True

            # this is the 1st data line and should look like '0.0000000000     -266.8384100090'
            if m == 4 and len(split_line_m) == 2:
                gaussian_dihedral_header_line_correct_bool_list[4] = True

        if False in gaussian_dihedral_header_line_correct_bool_list:
            raise TypeError(
                "ERROR: The Gaussian dihedral scan file or the 'Scan of Total Energy' file' "
                "does not have the proper header as follows:\n"
                "# Scan of Total Energy\n"
                "# X-Axis:  Scan Coordinate\n"
                "# Y-Axis:  Total Energy (Hartree)\n"
                "#                  X                   Y\n"
            )
            gaussian_angle_energy_file_correct = False

        else:
            gaussian_angle_energy_file_correct = True

    return gaussian_angle_energy_file_correct


def check_gaussian_optimized_coordinate_file_correct(
    gaussian_optimized_coordinate_path_and_name,
):
    """Check that the gaussian/Gausview file containing the optimized coordinates is formatted correctly.

    This checks the gaussian/Gausview containing the optimized coordinates is formatted correctly.

    The proper header format for the GausView/Gaussian output is as follows:

    | Row	Highlight	Display	Tag	Symbol	X	Y	Z

    Parameters
    ----------
    gaussian_optimized_coordinate_path_and_name: str
        The directory and filename of the Gaussian/Gausview optimized coordinates file.

    Returns
    -------
    gaussian_optimized_coordinate_header_line_correct_bool: True or TypeError
        True; if the file is formatted correctly
        TypeError; if the file is not formatted correctly
    """
    with open(gaussian_optimized_coordinate_path_and_name, "r") as fp:
        gaussian_optimized_coordinate_header_line_correct_bool = False
        gaussian_dihedral_header_file = fp.readlines()
        for m, line_m in enumerate(gaussian_dihedral_header_file):
            split_line_m = line_m.split()

            if m == 0 and len(split_line_m) == 8:
                if (
                    split_line_m[0] == "Row"
                    and split_line_m[1] == "Highlight"
                    and split_line_m[2] == "Display"
                    and split_line_m[3] == "Tag"
                    and split_line_m[4] == "Symbol"
                    and split_line_m[5] == "X"
                    and split_line_m[6] == "Y"
                    and split_line_m[7] == "Z"
                ):
                    gaussian_optimized_coordinate_header_line_correct_bool = (
                        True
                    )

        if gaussian_optimized_coordinate_header_line_correct_bool is False:
            raise TypeError(
                "ERROR: The Gaussian optimized coordinates file for each dihedral angle "
                "does not have the proper header as follows:\n"
                "Row	Highlight	Display	Tag	Symbol	X	Y	Z\n"
            )

    return gaussian_optimized_coordinate_header_line_correct_bool


def get_final_gaussian_output_file_data(
    qm_log_file_dict, manual_dihedral_atom_numbers_list
):
    """Get the gaussian/Gausview file data from the existing gaussian/Gausview files.

    This gets the gaussian/Gausview log file data for all the optimized configurations,
    moving it to the folder that will be analyzed.

    The proper header format for the GausView/Gaussian output is as follows:

    | Row	Highlight	Display	Tag	Symbol	X	Y	Z

    Parameters
    ----------
    qm_log_file_dict: dict, {str: [int, ..., int]}
        This is a dictionary comprised of a key (string) of the QM gaussian/Gausview file
        data path and name,
        and a list of integers, which are the QM optimization parameters to remove from
        the written data, in order of reading from each file. These can be seen in the
        order of the dictionary file name (strings).  These removed parameters allow
        users to remove any bad or repeated data points for the QM log file when needed.

        Example 1: {'path/gaussian_log_file.log': []}
        Uses all the optimized data points from the 'path/gaussian_log_file_data_path' file.

        Example 2: {'path/gaussian_log_file.log': [0, 23]}
        Uses all data points from the 'path/gaussian_log_file_data_path' file, except points
        0 and 23.  NOTE: Python counting starts at 0.
    manual_dihedral_atom_numbers_list: list, list of four (4) int (example: [3,2,1,5])
        This is a list of the dihedral atom numbers in order that were used for the dihedral
        fit. This information needs to be correct and in order to produce correct results.
        The values must be the same in all the combined files.

    Returns
    -------
    list of:
        all_dihedral_angle_degrees_list: list (nested list)
            This is the list of the optimized Gaussian dihedral angles (degrees) with the
            specific angles removed per the 'qm_log_file_dict' list
            (value) input. This is a nested list, with an inner list for every
            Gaussian log file or 'qm_log_file_dict' (key).
        all_energy_hartree_list: list (nested list)
            This is the list of the optimized Gaussian energies with the specific energies
            (hartree) removed per the 'qm_log_file_dict' list
            (value) input. This is a nested list, with an inner list for every
            Gaussian log file or 'qm_log_file_dict' (key).
        all_coordinates_ang_list: list (nested list)
            This is the list of the optimized Gaussian energies with the specific energies
            (hartree) removed per the 'qm_log_file_dict' list
            (value) input. This is a nested list, with an inner list for every
            Gaussian log file or 'qm_log_file_dict' (key).
        element_names_list: list
            This is the list of the Gaussian element names per the
            'qm_log_file_dict' list (value) input.
            The list length is dependant on the number of elements in the Gaussian log file
            These values are confirmed to be the same for all entered Gaussian log files.
        number_of_atoms: int
            This is the number of atoms in the Gaussian
            'qm_log_file_dict' list (value) input.
            These values are confirmed to be the same for all entered Gaussian log files.
        manual_dihedral_atom_numbers_list: list, list of four (4) int (example: [3,2,1,5])
        This is a list of the dihedral atom numbers in order that were used for the dihedral
        fit. This information needs to be correct and in order to produce correct results.
        The values must be the same in all the combined files.
    """
    if not isinstance(qm_log_file_dict, dict):
        raise TypeError(
            f"ERROR: In the 'get_final_gaussian_output_file_data' function, "
            f"the Gaussian 'log_files_and_entries_to_remove_dict' variable "
            f"is a {type(qm_log_file_dict)} not a dict."
        )

    else:
        for key_j, value_j in qm_log_file_dict.items():
            if not isinstance(key_j, str):
                raise ValueError(
                    f"ERROR: In the 'get_final_gaussian_output_file_data' function, "
                    f"the 'qm_log_file_dict' key "
                    f"'{key_j}' is a {type(key_j)} not a string."
                )
            if isinstance(value_j, list):
                for r_i in value_j:
                    if not isinstance(r_i, (int, np.int64)) or r_i < 0:
                        raise TypeError(
                            f"ERROR: In the 'get_final_gaussian_output_file_data' function, "
                            f"the 'qm_log_file_dict' values '{value_j}' "
                            f"are all not integers and >=0."
                        )
            else:
                raise TypeError(
                    f"ERROR: In the 'get_final_gaussian_output_file_data' function, "
                    f"the 'qm_log_file_dict' "
                    f"'{value_j}' is a {type(value_j)} not a list."
                )

    if (
        not isinstance(manual_dihedral_atom_numbers_list, list)
        or len(manual_dihedral_atom_numbers_list) != 4
    ):
        raise TypeError(
            "ERROR: The 'manual_dihedral_atom_numbers_list' is not a list of length 4."
        )

    elif (
        isinstance(manual_dihedral_atom_numbers_list, list)
        and len(manual_dihedral_atom_numbers_list) == 4
    ):
        for x_i in manual_dihedral_atom_numbers_list:
            if not isinstance(x_i, int):
                raise TypeError(
                    "ERROR: The 'manual_dihedral_atom_numbers_list' values are not integers."
                )

    all_coordinates_ang_list = []
    all_energy_hartree_list = []
    all_dihedral_angle_degrees_list = []

    all_dihedral_atom_numbers_list = []
    all_number_of_atoms_list = []
    all_element_names_list = []

    dihedral_counter = 0  # start at 1 and add a +1 initially
    for (
        direct_gaussian_folder_iter,
        entries_to_remove_list_iter,
    ) in qm_log_file_dict.items():
        # reset the dihedral used per file
        all_used_and_unused_dihedral_angle_degrees_list = []

        # run each log file
        # check the file is correctly formated
        direct_gaussian_angles_energy_formated_file_name_iter = (
            f"{direct_gaussian_folder_iter}/dihedral.txt"
        )
        check_gaussian_angle_energy_file_correct(
            direct_gaussian_angles_energy_formated_file_name_iter
        )

        with open(
            direct_gaussian_angles_energy_formated_file_name_iter, "r"
        ) as fp1:
            first_enerery_dihedral_file_data_header_lines = 4
            first_coord_file_data_header_lines = 1

            # get the 1st QM gaussian/Gausview file angles and dihedrals
            direct_gaussian_angles_energy_iter = fp1.readlines()
            for m, line_m in enumerate(direct_gaussian_angles_energy_iter):
                m_less_spacers = int(
                    m - first_enerery_dihedral_file_data_header_lines - 0
                )
                split_line_m = line_m.split()

                # get all used and unused dihedrals for later pulling coordinates in order
                if m_less_spacers >= 0 and len(split_line_m) == 2:
                    all_used_and_unused_dihedral_angle_degrees_list.append(
                        split_line_m[0]
                    )

                if (
                    m_less_spacers >= 0
                    and len(split_line_m) == 2
                    and m_less_spacers not in entries_to_remove_list_iter
                ):
                    all_dihedral_angle_degrees_list.append(split_line_m[0])
                    all_energy_hartree_list.append(split_line_m[1])

                elif m_less_spacers >= 0 and len(split_line_m) != 2:
                    raise ValueError(
                        f"ERROR: The directly input file {direct_gaussian_angles_energy_formated_file_name_iter} "
                        f"is not in the correct gaussian sytle format."
                    )

            total_used_and_unused_dihedrals_per_file = [
                len(all_used_and_unused_dihedral_angle_degrees_list)
            ]

            # iterate through the coord files (added +1 as file names start with 1)
            for dih_per_file_i in range(
                1, total_used_and_unused_dihedrals_per_file[-1] + 1
            ):
                direct_coordinates_ang_list_iter = []
                direct_number_of_atoms_list_iter = []
                direct_element_names_list_iter = []
                direct_dihedral_atom_numbers_list_iter = []
                if int(dih_per_file_i - 1) not in entries_to_remove_list_iter:
                    dihedral_counter += 1
                    direct_gaussian_coord_formated_file_name_iter = f"{direct_gaussian_folder_iter}/dihedral_coords_position_{dih_per_file_i}.txt"
                    check_gaussian_optimized_coordinate_file_correct(
                        direct_gaussian_coord_formated_file_name_iter
                    )

                    with open(
                        direct_gaussian_coord_formated_file_name_iter, "r"
                    ) as fp2:
                        # get the 1st QM gaussian/Gausview file angles and dihedrals
                        direct_gaussian_coord_iter = fp2.readlines()
                        for n, line_n in enumerate(direct_gaussian_coord_iter):
                            split_line_n = line_n.split()
                            n_less_spacers = int(
                                n - first_coord_file_data_header_lines - 0
                            )
                            if n_less_spacers >= 0 and len(split_line_n) == 8:
                                direct_coordinates_ang_list_iter.append(
                                    [
                                        float(split_line_n[5]),
                                        float(split_line_n[6]),
                                        float(split_line_n[7]),
                                    ]
                                )
                                direct_number_of_atoms_list_iter.append(
                                    float(split_line_n[3])
                                )
                                direct_element_names_list_iter.append(
                                    split_line_n[4]
                                )
                                direct_dihedral_atom_numbers_list_iter.append(
                                    manual_dihedral_atom_numbers_list
                                )

                            elif n_less_spacers >= 0 and len(split_line_n) != 8:
                                raise ValueError(
                                    f"ERROR: The directly input file {direct_gaussian_coord_formated_file_name_iter} "
                                    f"is not in the correct gaussian sytle format"
                                )

                        # check the values against the past ones
                        if (
                            len(all_coordinates_ang_list) == 0
                            and len(all_number_of_atoms_list) == 0
                            and len(all_element_names_list) == 0
                        ):
                            all_coordinates_ang_list.append(
                                direct_coordinates_ang_list_iter
                            )
                            all_number_of_atoms_list.append(
                                direct_number_of_atoms_list_iter
                            )
                            all_element_names_list.append(
                                direct_element_names_list_iter
                            )
                            all_dihedral_atom_numbers_list.append(
                                direct_dihedral_atom_numbers_list_iter
                            )

                        else:
                            all_coordinates_ang_list.append(
                                direct_coordinates_ang_list_iter
                            )

                            # check if the number_of_atoms are the same for all Gaussian files
                            all_number_of_atoms_list = mdf_math.check_previous_qm_values_match(
                                all_number_of_atoms_list,
                                direct_number_of_atoms_list_iter,
                                "number of atoms",
                                "Direct gaussian output file",
                                direct_gaussian_coord_formated_file_name_iter,
                            )

                            # check if the nelement_names are the same for all Gaussian files
                            all_element_names_list = mdf_math.check_previous_qm_values_match(
                                all_element_names_list,
                                direct_element_names_list_iter,
                                "element names",
                                "Direct gaussian output file",
                                direct_gaussian_coord_formated_file_name_iter,
                            )

                            # check if the nelement_names are the same for all Gaussian files
                            all_dihedral_atom_numbers_list = mdf_math.check_previous_qm_values_match(
                                all_dihedral_atom_numbers_list,
                                direct_dihedral_atom_numbers_list_iter,
                                "dihedral atom numbers",
                                "Direct gaussian output file",
                                direct_gaussian_coord_formated_file_name_iter,
                            )

    element_names_list = all_element_names_list[0]
    number_of_atoms = len(element_names_list)

    return [
        all_dihedral_angle_degrees_list,
        all_energy_hartree_list,
        all_coordinates_ang_list,
        element_names_list,
        number_of_atoms,
        manual_dihedral_atom_numbers_list,
    ]


def get_gaussian_log_file_data(
    qm_log_file_dict,
):
    """Get the gaussian/Gausview file data from the log filefor all the optimized configurations.

    This gets the gaussian/Gausview log file data for all the optimized configurations,
    allowing the data to be analyzed further.

    The proper header format for the GausView/Gaussian output is as follows:

    | Row	Highlight	Display	Tag	Symbol	X	Y	Z

    Parameters
    ----------
    qm_log_file_dict: dict, {str: [int, ..., int]}
        This is a dictionary comprised of a key (string) of the QM log file path and name,
        and a list of integers, which are the QM optimization parameters to remove from
        the written data, in order of reading from each file. These can be seen in the
        order of the dictionary file name (strings).  These removed parameters allow
        users to remove any bad or repeated data points for the QM log file when needed.

        Example 1: {'path/gaussian_log_file.log': []}
        Uses all the optimized data points from the 'path/gaussian_log_file.log' file.

        Example 2: {'path/gaussian_log_file.log': [0, 23]}
        Uses all data points from the 'path/gaussian_log_file.log' file, except points
        0 and 23.  NOTE: Python counting starts at 0.

    Returns
    -------
    list of:
        all_dihedral_angle_degrees_list: list (nested list)
            This is the list of the optimized Gaussian dihedral angles (degrees) with the
            specific angles removed per the 'qm_log_file_dict' list
            (value) input. This is a nested list, with an inner list for every
            Gaussian log file or 'qm_log_file_dict' (key).
        all_energy_hartree_list: list (nested list)
            This is the list of the optimized Gaussian energies with the specific energies
            (hartree) removed per the 'qm_log_file_dict' list
            (value) input. This is a nested list, with an inner list for every
            Gaussian log file or 'qm_log_file_dict' (key).
        all_coordinates_ang_list: list (nested list)
            This is the list of the optimized Gaussian energies with the specific energies
            (hartree) removed per the 'qm_log_file_dict' list
            (value) input. This is a nested list, with an inner list for every
            Gaussian log file or 'qm_log_file_dict' (key).
        element_names_list: list
            This is the list of the Gaussian element names per the
            'qm_log_file_dict' list (value) input.
            The list length is dependant on the number of elements in the Gaussian log file
            These values are confirmed to be the same for all entered Gaussian log files.
        number_of_atoms: int
            This is the number of atoms in the Gaussian
            'qm_log_file_dict' list (value) input.
            These values are confirmed to be the same for all entered Gaussian log files.
        dihedral_atom_numbers_list: list of 4 integers
            This is the list of the dihedral atom numbers used in the Gaussian
            dihedral scan, which are taken from the Gaussian
            'qm_log_file_dict' list (value) input.
            These values are confirmed to be the same for all entered Gaussian log files.

    """
    if not isinstance(qm_log_file_dict, dict):
        raise TypeError(
            f"ERROR: In the 'get_gaussian_log_file_data' function, "
            f"the Gaussian 'log_files_and_entries_to_remove_dic' variable "
            f"is a {type(qm_log_file_dict)} not a dict."
        )

    else:
        for key_j, value_j in qm_log_file_dict.items():
            if not isinstance(key_j, str):
                raise ValueError(
                    f"ERROR: In the 'get_gaussian_log_file_data' function, "
                    f"the 'qm_log_file_dict' key "
                    f"'{key_j}' is a {type(key_j)} not a string."
                )
            if isinstance(value_j, list):
                for r_i in value_j:
                    if not isinstance(r_i, int) or r_i < 0:
                        raise TypeError(
                            f"ERROR: In the 'get_gaussian_log_file_data' function, "
                            f"the 'qm_log_file_dict' values '{value_j}' "
                            f"are all not integers and >=0."
                        )
            else:
                raise TypeError(
                    f"ERROR: In the 'get_gaussian_log_file_data' function, "
                    f"the 'qm_log_file_dict' "
                    f"'{value_j}' is a {type(value_j)} not a list."
                )

    all_coordinates_ang_list = []
    all_energy_hartree_list = []
    all_dihedral_angle_degrees_list = []

    all_dihedral_atom_numbers_list = []
    all_molecule_charge_list = []
    all_molecule_multiplicity_list = []
    all_number_of_atoms_list = []
    all_element_names_list = []

    for (
        log_file_iter,
        entries_to_remove_list_iter,
    ) in qm_log_file_dict.items():
        # run each log file
        with open(log_file_iter, "r") as fp:
            dihedral_scan_line = None
            dihedral_atom_numbers_list = None
            number_of_optimized_scans = None
            scan_degrees_per_scan = None
            molecule_charge = None
            molecule_multiplicity = None
            number_of_atoms = None

            energy_hartree_list_iter = []

            optimized_energy_hartree_list = []
            optimized_coordinates_ang_list = []
            optimized_dihedral_angle_degrees_list = []

            element_names_list = []
            get_element_names_bool = True
            # need to store the last coordinates and only add to the list once the 'Optimized Parameters' is seen
            gaussian_log_file = fp.readlines()
            for m, line_m in enumerate(gaussian_log_file):
                split_line_m = line_m.split()

                # **********************************
                # get the starting dihedral scan location, elements, charge, multiplicity and number of atoms (START)
                # **********************************
                # get and mark the initial system setup for a dihedral scan
                if (
                    len(split_line_m) >= 2
                    and split_line_m[0] == "Dihedral"
                    and split_line_m[1] == "Scan"
                ):
                    dihedral_scan_line = m

                # get the gaussian molecule charge and multiplicity
                if (
                    len(split_line_m) >= 4
                    and split_line_m[0] == "Charge"
                    and split_line_m[1] == "="
                    and split_line_m[3] == "Multiplicity"
                    and split_line_m[4] == "="
                ):
                    molecule_charge = split_line_m[2]
                    molecule_multiplicity = split_line_m[5]

                    # check if the molecule charges are the same for all Gaussian files
                    all_molecule_charge_list = (
                        mdf_math.check_previous_qm_values_match(
                            all_molecule_charge_list,
                            molecule_charge,
                            "molecule charge",
                            "Gaussian",
                            log_file_iter,
                        )
                    )
                    # check if the molecule multiplicity are the same for all Gaussian files
                    all_molecule_multiplicity_list = (
                        mdf_math.check_previous_qm_values_match(
                            all_molecule_multiplicity_list,
                            molecule_multiplicity,
                            "molecule multiplicity",
                            "Gaussian",
                            log_file_iter,
                        )
                    )

                # use the NAtoms= to get number of atoms
                if (
                    len(split_line_m) > 2
                    and number_of_atoms is None
                    and split_line_m[0] == "NAtoms="
                ):
                    number_of_atoms = int(split_line_m[1])

                # get the element names
                spaces_to_element = 2
                if (
                    len(split_line_m) == 5
                    and get_element_names_bool is True
                    and number_of_atoms is not None
                    and split_line_m[0] == "Condensed"
                    and split_line_m[1] == "to"
                    and split_line_m[2] == "atoms"
                    and split_line_m[3] == "(all"
                    and split_line_m[4] == "electrons):"
                    and len(gaussian_log_file[m + spaces_to_element].split())
                    > 2
                ):
                    for spaces_to_element_i in range(
                        spaces_to_element, spaces_to_element + number_of_atoms
                    ):
                        element_names_list.append(
                            gaussian_log_file[m + spaces_to_element_i].split()[
                                1
                            ]
                        )

                    get_element_names_bool = False

                # **********************************
                # get the starting dihedral scan location, elements, charge, multiplicity and number of atoms (END)
                # **********************************

                # **********************************
                # get the dihedral atom numbers  (START)
                # **********************************

                # get the dihedral atom numbers (in Gaussian numbering starts at 1)
                # and the dihedral degrees per scan
                if len(split_line_m) == 8 and dihedral_scan_line is not None:
                    if (
                        split_line_m[0] == "D"
                        and split_line_m[5] == "S"
                        and gaussian_log_file[m - 1].split()[0] == "The"
                        and gaussian_log_file[m - 1].split()[1] == "following"
                    ):
                        dihedral_atom_numbers_list = [
                            int(split_line_m[1]),
                            int(split_line_m[2]),
                            int(split_line_m[3]),
                            int(split_line_m[4]),
                        ]
                        scan_degrees_per_scan = float(split_line_m[7])

                        # check if the all the dihedral atom numbers  are the same for all Gaussian files
                        all_dihedral_atom_numbers_list = (
                            mdf_math.check_previous_qm_values_match(
                                all_dihedral_atom_numbers_list,
                                dihedral_atom_numbers_list,
                                "dihedral atom numbers list",
                                "Gaussian",
                                log_file_iter,
                            )
                        )

                # **********************************
                # get the dihedral atom numbers  (END)
                # **********************************

                # **********************************
                # get the number of optimized scans  (START)
                # **********************************
                # get the number of optimized scans
                if (
                    len(split_line_m) >= 6
                    and split_line_m[0] == "Number"
                    and split_line_m[1] == "of"
                    and split_line_m[2] == "optimizations"
                    and split_line_m[3] == "in"
                    and split_line_m[4] == "scan="
                ):
                    number_of_optimized_scans = int(split_line_m[5])
                # **********************************
                # get the number of optimized scans  (END)
                # **********************************

                # **********************************
                # get the optimized coordinates  (START)
                # **********************************

                # get the last coordinate point
                if (
                    len(split_line_m) == 2
                    and split_line_m[0] == "Standard"
                    and split_line_m[1] == "orientation:"
                    and gaussian_log_file[m + 2].split()[0] == "Center"
                    and gaussian_log_file[m + 2].split()[1] == "Atomic"
                    and gaussian_log_file[m + 2].split()[2] == "Atomic"
                    and gaussian_log_file[m + 2].split()[3] == "Coordinates"
                    and gaussian_log_file[m + 2].split()[4] == "(Angstroms)"
                    and gaussian_log_file[m + 3].split()[0] == "Number"
                    and gaussian_log_file[m + 3].split()[1] == "Number"
                    and gaussian_log_file[m + 3].split()[2] == "Type"
                    and gaussian_log_file[m + 3].split()[3] == "X"
                    and gaussian_log_file[m + 3].split()[4] == "Y"
                    and gaussian_log_file[m + 3].split()[5] == "Z"
                ):
                    if isinstance(number_of_atoms, int):
                        coord_list_iter = []
                        for atom_i in range(0, number_of_atoms):
                            spaces_to_1st_atom_coord = 5
                            coor_line_m = gaussian_log_file[
                                m + spaces_to_1st_atom_coord + atom_i
                            ].split()

                            if len(coor_line_m) == 6:
                                coord_list_iter.append(
                                    [
                                        float(coor_line_m[3]),
                                        float(coor_line_m[4]),
                                        float(coor_line_m[5]),
                                    ]
                                )

                        last_coordinates_ang_list = coord_list_iter

                if len(optimized_dihedral_angle_degrees_list) > len(
                    optimized_coordinates_ang_list
                ):
                    optimized_coordinates_ang_list.append(
                        last_coordinates_ang_list
                    )

                # **********************************
                # get the optimized coordinates  (END)
                # **********************************

                # **********************************
                # get energy data in Hartree (START)
                # **********************************

                if (
                    len(split_line_m) >= 5
                    and split_line_m[0] == "SCF"
                    and split_line_m[1] == "Done:"
                    and split_line_m[3] == "="
                ):
                    energy_hartree_list_iter.append(split_line_m[4])

                # only collect parameters (optimized_energy_hartree_list) after the system is optimized
                if len(optimized_dihedral_angle_degrees_list) > len(
                    optimized_energy_hartree_list
                ):
                    optimized_energy_hartree_list.append(
                        energy_hartree_list_iter[-1]
                    )

                # **********************************
                # get energy data in Hartree (START)
                # **********************************

                # **********************************
                # get actual degree from optimization rounded to 3 decimal places (START)
                # **********************************
                if (
                    len(split_line_m) == 8
                    and split_line_m[0] == "!"
                    and split_line_m[2] == f"D("
                    f"{dihedral_atom_numbers_list[0]},"
                    f"{dihedral_atom_numbers_list[1]},"
                    f"{dihedral_atom_numbers_list[2]},"
                    f"{dihedral_atom_numbers_list[3]}"
                    f")"
                    and split_line_m[4] == "-DE/DX"
                    and split_line_m[5] == "="
                    and split_line_m[7] == "!"
                ):
                    optimized_dihedral_angle_degrees_list.append(
                        str(np.round(float(split_line_m[3]), decimals=3))
                    )
                # **********************************
                # get actual degree from optimization  rounded to 3 decimal places (End)
                # **********************************

        if (
            number_of_optimized_scans != len(optimized_coordinates_ang_list)
            or number_of_optimized_scans != len(optimized_energy_hartree_list)
            or number_of_optimized_scans
            != len(optimized_dihedral_angle_degrees_list)
        ):
            raise ValueError(
                f"ERROR: The lengths of "
                f"the optimized coordinate list 'len(optimized_coordinates_ang_list)' "
                f"= {len(optimized_coordinates_ang_list)}, "
                f"or the optimized energy list 'len(optimized_energy_hartree_list)' "
                f"= {len(optimized_energy_hartree_list)}, "
                f"or the optimized dihedral list 'len(optimized_dihedral_angle_degrees_list)' "
                f"= {len(optimized_dihedral_angle_degrees_list)}, "
                f"is not the same as "
                f"the Gaussian file's number of optimized scans 'len(number_of_optimized_scans)' "
                f"= {len(number_of_optimized_scans)}."
            )

        if (
            isinstance(entries_to_remove_list_iter, list)
            and len(entries_to_remove_list_iter) != 0
        ):
            if max(entries_to_remove_list_iter) > number_of_optimized_scans - 1:
                raise ValueError(
                    f"ERROR: For the '{log_file_iter}' log file, "
                    f"the 'entries_to_remove_list_iter' list '{entries_to_remove_list_iter}' "
                    f"has values greater than the number of optimized scans in the dihedral. "
                    f"The only allowed values are from 0 to {number_of_optimized_scans-1}."
                )

        for entry_no_i in range(0, len(optimized_coordinates_ang_list)):
            if entry_no_i not in entries_to_remove_list_iter:
                all_coordinates_ang_list.append(
                    optimized_coordinates_ang_list[entry_no_i]
                )
                all_energy_hartree_list.append(
                    optimized_energy_hartree_list[entry_no_i]
                )
                all_dihedral_angle_degrees_list.append(
                    optimized_dihedral_angle_degrees_list[entry_no_i]
                )

                # check if theelement names  are the same for all Gaussian files
                all_element_names_list = (
                    mdf_math.check_previous_qm_values_match(
                        all_element_names_list,
                        element_names_list,
                        "element names ",
                        "Gaussian",
                        log_file_iter,
                    )
                )

                all_number_of_atoms_list = (
                    mdf_math.check_previous_qm_values_match(
                        all_number_of_atoms_list,
                        number_of_atoms,
                        "number of atoms",
                        "Gaussian",
                        log_file_iter,
                    )
                )

    return [
        all_dihedral_angle_degrees_list,
        all_energy_hartree_list,
        all_coordinates_ang_list,
        element_names_list,
        number_of_atoms,
        dihedral_atom_numbers_list,
    ]


def write_qm_data_files(
    qm_log_file_dict,
    manual_dihedral_atom_numbers_list=None,
    qm_engine="gaussian",
):
    """Write out the optimized QM simulation data from the QM log files.

    This extracts and writes out the ptimized QM simulation data from the QM log files.

    Parameters
    ----------
    qm_log_file_dict: dict, {str: [int, ..., int]}
        * qm_engine="gaussian"
            This is a dictionary comprised of a key (string) of the QM log file path and name,
            and a list of integers, which are the QM optimization parameters to remove from
            the written data, in order of reading from each file. These can be seen in the
            order of the dictionary file name (strings).  These removed parameters allow
            users to remove any bad or repeated data points for the QM log file when needed.

            Example 1: {'path/gaussian_log_file.log': []}

            Uses all the optimized data points from the 'path/gaussian_log_file.log' file.

            Example 2: {'path/gaussian_log_file.log': [0, 23]}
            Uses all data points from the 'path/gaussian_log_file.log' file, except points
            0 and 23.  NOTE: Python counting starts at 0.

        * qm_engine="gaussian_style_final_files"
            This is a dictionary comprised of a key (string) of the  file paths to the
            Gaussian style final formatted files, and a list of integers, which are the
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

    qm_engine: str, default='gaussian' (options = 'gaussian')
        The QM simulation engine that was utilized also tells the log file readers
        what QM log file read to use for the analysis.

    Notes
    -------
    files are written to the created 'extracted_gaussian_data' folder:
        - 'dihedral.txt' file is in the standard Gaussian/Gausview format and
        contains the optimized scanned/rotated dihedral angle and energy in Hartree
        energy units.
        - 'dihedral_coords_position_n.txt' file is in the standard Gaussian/Gausview
        format. There are n coordinate files written (n starts at 1), one for each
        optimized dihedral angle/energy, wwhose output is numbered in the same
        order as the 'dihedral.txt' file data. The coordinate file distances are
        in Angstrom units.
    """
    # delete any existing directories and make a new one
    gaussian_directory_name = "extracted_gaussian_data"
    if os.path.isdir(gaussian_directory_name):
        shutil.rmtree(gaussian_directory_name)
    os.mkdir(gaussian_directory_name)

    if qm_engine == "gaussian":
        # extract the required data to write the gaussian style formatted output files
        [
            all_dihedral_angle_degrees_list,
            all_energy_hartree_list,
            all_coordinates_ang_list,
            element_names_list,
            number_of_atoms,
            dihedral_atom_numbers_list,
        ] = get_gaussian_log_file_data(qm_log_file_dict)

    elif qm_engine == "gaussian_style_final_files":
        if manual_dihedral_atom_numbers_list == None:
            raise TypeError(
                "ERROR: The 'dihedral_atom_numbers_list' is not a list of length 4, "
                "and needs entered when using 'gaussian_style_final_files'."
            )
        [
            all_dihedral_angle_degrees_list,
            all_energy_hartree_list,
            all_coordinates_ang_list,
            element_names_list,
            number_of_atoms,
            dihedral_atom_numbers_list,
        ] = get_final_gaussian_output_file_data(
            qm_log_file_dict,
            manual_dihedral_atom_numbers_list=manual_dihedral_atom_numbers_list,
        )

    else:
        raise ValueError(
            f"ERROR: In the 'write_qm_data_files' function, "
            f"the 'qm_engine' variable = {qm_engine}, which is not "
            f"any of the available options."
        )
    # write the gaussian style formatted angle (degrees) and energy output files
    output_file_dihedral_energy = open(
        f"extracted_gaussian_data/dihedral.txt", "w"
    )

    output_file_dihedral_energy.write(
        f"# Scan of Total Energy\n"
        f"# X-Axis:  Scan Coordinate\n"
        f"# Y-Axis:  Total Energy (Hartree)\n"
        f"#                  X                   Y\n"
    )

    for qm_i in range(0, len(all_dihedral_angle_degrees_list)):
        output_file_dihedral_energy.write(
            f"{all_dihedral_angle_degrees_list[qm_i]: >20}"
            f"{all_energy_hartree_list[qm_i]: >20}\n"
        )

    output_file_dihedral_energy.close()

    # write the gaussian style formatted coordinates file
    for qm_j in range(0, len(all_dihedral_angle_degrees_list)):
        output_file_dihedral_coordinates = open(
            f"extracted_gaussian_data/dihedral_coords_position_{int(qm_j+1)}.txt",
            "w",
        )

        output_file_dihedral_coordinates.write(
            f"Row	Highlight	Display	Tag	Symbol	X	Y	Z\n"
        )

        for qm_k in range(0, number_of_atoms):
            # atom_number_i starts at 1
            atom_number_i = int(qm_k + 1)
            highlight = "No"
            display = "Show"
            output_file_dihedral_coordinates.write(
                "{: <8}"
                "{: <8}"
                "{: <8}"
                "{: <8}"
                "{: <8}"
                "{: <16}"
                "{: <16}"
                "{: <16}\n".format(
                    atom_number_i,
                    highlight,
                    display,
                    atom_number_i,
                    element_names_list[qm_k],
                    all_coordinates_ang_list[qm_j][qm_k][0],
                    all_coordinates_ang_list[qm_j][qm_k][1],
                    all_coordinates_ang_list[qm_j][qm_k][2],
                )
            )

        output_file_dihedral_coordinates.close()


def get_matching_dihedral_info_and_opls_fitting_data(
    fit_dihedral_atom_types,
    psf_path_and_filename,
    qm_log_file_dict,
    mol2_file,
    qm_engine="gaussian",
    manual_dihedral_atom_numbers_list=None,
):
    """Get all dihedral angles from the dihedrals which match the fitted dihedrals atom types/classes.

    This gets all dihedral angles from the dihedrals which match the fitted
    dihedrals atom types/classes. This is required to properly fit the rotated QM
    dihedral, if there are any other matching dihedrals which need to be accounted for.
    If they are not accounted for, then all dihedrals will add the forces to the single
    dihedral fit, which will overfit the energy of a single dihedral.

    NOTE: In the PSF file and Gaussian, all atom numbers start at one (1).


    Parameters
    ----------
    fit_dihedral_atom_types: list of 4 str
        The atom types/classes of the dihedral's atoms that is being fit, in order,
        like a proper dihedral is written. The code will also check the reverse order
        of the dihedral.

        Example 1: Ethane dihedral = ['HC', 'CT', 'CT', 'CH']

        Example 2: Ethane dihedral with wildcards = ['X', 'CT', 'CT', 'X']
                   or ['*', CT, 'CT', '*'] or ['', CT, 'CT', '']

    psf_path_and_filename: str
        The path and filename of the PSF file, which shall be used to extract and map the
        dihedrals atom types/classes and all their associated atom number.
    mol2_file: str
        The mol2 file which matches the element, atom type, bonded connnections,
        the 'EXACT ATOM ORDER AND CONFIGURATION AS IN THE QM SIMULATION INPUT FILES'.
        This is required to know the MM bonding in the atoms, because QM simulations
        do not explictly specify the system bonds.
    qm_log_file_dict: dict, {str: [int, ..., int]}
        * qm_engine="gaussian"
            This is a dictionary comprised of a key (string) of the QM log file path and name,
            and a list of integers, which are the QM optimization parameters to remove from
            the written data, in order of reading from each file. These can be seen in the
            order of the dictionary file name (strings).  These removed parameters allow
            users to remove any bad or repeated data points for the QM log file when needed.

            Example 1: {'path/gaussian_log_file.log': []}

            Uses all the optimized data points from the 'path/gaussian_log_file.log' file.

            Example 2: {'path/gaussian_log_file.log': [0, 23]}
            Uses all data points from the 'path/gaussian_log_file.log' file, except points
            0 and 23.  NOTE: Python counting starts at 0.

        * qm_engine="gaussian_style_final_files"
            This is a dictionary comprised of a key (string) of the  file paths to the
            Gaussian style final formatted files, and a list of integers, which are the
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

    manual_dihedral_atom_numbers_list: list of 4 integers, default=None
        NOTE: Only needed for qm_engine="gaussian_style_final_files"

        This is a list of the dihedral atom numbers in order that were used for the dihedral
        fit. This information needs to be correct and in order to produce correct results.
        The values must be the same in all the combined files.

    Returns
    -------
    matching_dihedral_types_by_atom_numbers_list: nested list
        A list of the four atom numbers (integers) that match the selected
        'fit_dihedral_atom_types', extracted from the PSF file, including the dihedral being fit.

        The list length is the number of dihedrals which match the QM scanned atom/bead types,
        with the nested list has four (4) integers

    matching_dihedral_types_by_atom_type_list: nested list
        A list of the four atom numbers that match the selected
        'fit_dihedral_atom_types', extracted from the PSF file, including the dihedral being fit.

        The list length is the number of dihedrals which match the QM scanned atom/bead types,
        with the nested list has four (4) strings.

    all_matching_dihedral_coordinates_angstroms_added_to_k_values_list: nested lists of 3
        The coordinates from all dihedrals, which match the fitted dihedral,
        including the dihedral being fit.

        The list length is the number of dihedrals which match the QM scanned atom/bead types,
        with the next (2nd) nested list consisting of (4) lists, being one (1) atom/bead for,
        with a list length the same as the 'matching_dihedral_types_by_atom_numbers_list' and
        'matching_dihedral_types_by_atom_type_list' lists.
        The last (3rd) nested list is the x, y, and z coordinates [x, y, z], with a list
        length of three (3)

        Example:
        [
           [
               [qm_scan_1_dih_1_coor_x, qm_scan_1_dih_1_coor_y, qm_scan_1_dih_1_coor_z],
               [qm_scan_1_dih_2_coor_x, qm_scan_1_dih_2_coor_y, qm_scan_2_dih_1_coor_z],
           ]
        ...
           [
               [qm_scan_n_dih_1_coor_x, qm_scan_n_dih_1_coor_y, qm_scan_n_dih_1_coor_z],
               [qm_scan_n_dih_2_coor_x, qm_scan_n_dih_2_coor_y, qm_scan_n_dih_1_coor_z],
           ]
         ]

    all_matching_dihedral_phi_degrees_added_to_k_values_list: nested list
        The angles/phis, in degrees, from all dihedrals, which match the fitted dihedral,
        including the dihedral being fit.

        The list length is the number of dihedrals which match the QM scanned atom/bead types,
        with the next nested list consisting of the all the dihedral angles/phis that have
        the same atom type/class as the fitted dihedral, with a list length the same as the
        'matching_dihedral_types_by_atom_numbers_list' and
        'matching_dihedral_types_by_atom_type_list' lists.

    all_sum_opls_const_1_plus_or_minus_cos_n_list: nested list, nested list
        A list of the OPLS (1 +/- cos(n * phi)) values for all k-values, which
        is required to fit the data, especially when multiple dihedrals of the same
        atom types/classes exist in the molecule that the dihedral is fit too.

        all_sum_opls_const_1_plus_or_minus_cos_n_list = [
        const_1_minus_Cos_0_phi,
        const_1_plus_Cos_1_phi,
        const_1_minus_Cos_2_phi,
        const_1_plus_Cos_3_phi,
        const_1_minus_Cos_4_phi
        ]

        Example 0:
        const_1_minus_Cos_0_phi  = 0 , since k0 is not used in this form

        Example 1:
        const_1_plus_Cos_1_phi  = sum of all phis in the list [(1 + cos(1 * phi))] in k1 * (1 + cos(1 * phi))

        Example 2:
        const_1_minus_Cos_2_phi = sum of all phis in the list [(1 - cos(2 * phi))] in k2 * (1 - cos(2 * phi))

        Example 3:
        const_1_plus_Cos_3_phi  = sum of all phis in the list [(1 + cos(3 * phi))] in k3 * (1 + cos(3 * phi))

        Example 4:
        const_1_minus_Cos_4_phi = sum of all phis in the list [(1 - cos(4 * phi))] in k4 * (1 - cos(4 * phi))

        The list length is the number of dihedrals which match the QM scanned atom/bead types,
        with the next nested list being the 'all_sum_opls_const_1_plus_or_minus_cos_n_list'
        as described.
    """

    if isinstance(fit_dihedral_atom_types, list):
        if len(fit_dihedral_atom_types) != 4:
            raise TypeError(
                f"ERROR: The 'fit_dihedral_atom_types' length = {len(fit_dihedral_atom_types)}, not 4."
            )

        for fit_i in fit_dihedral_atom_types:
            if not isinstance(fit_i, str):
                raise TypeError(
                    f"ERROR: The 'fit_dihedral_atom_types' variable = {fit_dihedral_atom_types} does not "
                    f"only contain strings. The value {type(fit_i)} is not a string."
                )
    else:
        raise TypeError(
            f"ERROR: The 'fit_dihedral_atom_types' variable is a {type(fit_dihedral_atom_types)}, not a list"
        )

    read_psf_file = open(f"{psf_path_and_filename}", "r").readlines()
    atom_number_to_atom_type_dict = {}
    all_dihedral_numbers_as_read_list = []
    for psf_iter, psf_line_iter in enumerate(read_psf_file):
        psf_splitline_iter = psf_line_iter.split()

        # get the PSF files atom number to atom type/class map dictionary
        if len(psf_splitline_iter) == 2 and psf_splitline_iter[1] == "!NATOM":
            number_of_atoms_iter = int(psf_splitline_iter[0])

            if len(read_psf_file[psf_iter + 1].split()) == 8:
                for atom_i in range(0, number_of_atoms_iter):
                    atom_number_to_atom_type_dict.update(
                        {
                            int(
                                read_psf_file[psf_iter + atom_i + 1].split()[0]
                            ): read_psf_file[psf_iter + atom_i + 1].split()[5]
                        }
                    )

        # get the PSF files dihedrals
        if (
            len(psf_splitline_iter) == 3
            and psf_splitline_iter[1] == "!NPHI:"
            and psf_splitline_iter[2] == "dihedrals"
        ):
            # Get the number of dihedral lines with 2 dihedrals per line (full_lines).
            # See if the last line has 1 or 2 (half_full_lines), which are only present in the last line.
            full_lines = int(int(psf_splitline_iter[0]) / 2)
            half_full_lines = int(int(psf_splitline_iter[0]) % 2)

            if len(read_psf_file[psf_iter + 1].split()) == 8:
                for dihedral_i in range(0, full_lines):
                    # get the first dihedral in the full line
                    all_dihedral_numbers_as_read_list.append(
                        [
                            int(
                                read_psf_file[
                                    psf_iter + dihedral_i + 1
                                ].split()[0]
                            ),
                            int(
                                read_psf_file[
                                    psf_iter + dihedral_i + 1
                                ].split()[1]
                            ),
                            int(
                                read_psf_file[
                                    psf_iter + dihedral_i + 1
                                ].split()[2]
                            ),
                            int(
                                read_psf_file[
                                    psf_iter + dihedral_i + 1
                                ].split()[3]
                            ),
                        ]
                    )
                    # get the second dihedral in the full line
                    all_dihedral_numbers_as_read_list.append(
                        [
                            int(
                                read_psf_file[
                                    psf_iter + dihedral_i + 1
                                ].split()[4]
                            ),
                            int(
                                read_psf_file[
                                    psf_iter + dihedral_i + 1
                                ].split()[5]
                            ),
                            int(
                                read_psf_file[
                                    psf_iter + dihedral_i + 1
                                ].split()[6]
                            ),
                            int(
                                read_psf_file[
                                    psf_iter + dihedral_i + 1
                                ].split()[7]
                            ),
                        ]
                    )

                    # get the first dihedral in the line in the half full line
                    if half_full_lines == 1 and full_lines == dihedral_i + 1:
                        all_dihedral_numbers_as_read_list.append(
                            [
                                int(
                                    read_psf_file[
                                        psf_iter + dihedral_i + 2
                                    ].split()[0]
                                ),
                                int(
                                    read_psf_file[
                                        psf_iter + dihedral_i + 2
                                    ].split()[1]
                                ),
                                int(
                                    read_psf_file[
                                        psf_iter + dihedral_i + 2
                                    ].split()[2]
                                ),
                                int(
                                    read_psf_file[
                                        psf_iter + dihedral_i + 2
                                    ].split()[3]
                                ),
                            ]
                        )

    # get the reversed order of the atom types/classes that you want to compare/exact from the PSF file
    fit_dihedral_atom_types_reversed = [
        fit_dihedral_atom_types[3],
        fit_dihedral_atom_types[2],
        fit_dihedral_atom_types[1],
        fit_dihedral_atom_types[0],
    ]

    # check for a forward and revers match, including wildcard combinations, for the
    # fitted dihedral parameter types ('fit_dihedral_atom_types')
    matching_dihedral_types_by_atom_numbers_list = []
    matching_dihedral_types_by_atom_type_list = []
    for iter_p, dih_num_p in enumerate(all_dihedral_numbers_as_read_list):
        # convert from atom numbers to atom types/classes
        dih_atom_type_p = []
        for iter_q in range(0, len(dih_num_p)):
            dih_atom_type_p.append(
                atom_number_to_atom_type_dict[dih_num_p[iter_q]]
            )

        # check for the input dihedral fit and reversed dihedral fit
        match_dihedral_atom_types_4_bool_list = [False, False, False, False]
        match_reverse_dihedral_atom_types_4_bool_list = [
            False,
            False,
            False,
            False,
        ]
        for iter_r, dih_atom_type_r in enumerate(dih_atom_type_p):
            if dih_atom_type_r in ["x", "X", "*", ""]:
                match_dihedral_atom_types_4_bool_list[iter_r] = True
                match_reverse_dihedral_atom_types_4_bool_list[iter_r] = True

            elif dih_atom_type_r == fit_dihedral_atom_types[iter_r]:
                match_dihedral_atom_types_4_bool_list[iter_r] = True

            elif dih_atom_type_r == fit_dihedral_atom_types_reversed[iter_r]:
                match_reverse_dihedral_atom_types_4_bool_list[iter_r] = True

        # check if either the forward or reversed dihedral matches the read dihedral type/class.
        # if so save the atom numbers
        # matching_dihedral_types_by_atom_numbers_list = []
        # matching_dihedral_types_by_atom_type_list = []
        if (
            False not in match_dihedral_atom_types_4_bool_list
            or False not in match_reverse_dihedral_atom_types_4_bool_list
        ):
            matching_dihedral_types_by_atom_numbers_list.append(dih_num_p)
            matching_dihedral_types_by_atom_type_list.append(dih_atom_type_p)

        # check if the QM and PSF file match in atom numbers
        # NOTE: gaussian and the PSF file start atom counting at 1.

        # extract the data from the QM log file
        if qm_engine == "gaussian":
            # NOTE: gaussian start atom counting at 1.
            [
                dihedral_angle_degrees_list,
                energy_hartree_list,
                coordinates_ang_list,
                all_element_names_list,
                number_of_atoms_list,
                dihedral_atom_numbers_list,
            ] = get_gaussian_log_file_data(qm_log_file_dict)

        elif qm_engine == "gaussian_style_final_files":
            if manual_dihedral_atom_numbers_list == None:
                raise TypeError(
                    "ERROR: The 'manual_dihedral_atom_numbers_list' is not a list of length 4, "
                    "and needs entered when using 'direct_gaussian_final_files'."
                )
            [
                dihedral_angle_degrees_list,
                energy_hartree_list,
                coordinates_ang_list,
                all_element_names_list,
                number_of_atoms_list,
                dihedral_atom_numbers_list,
            ] = get_final_gaussian_output_file_data(
                qm_log_file_dict,
                manual_dihedral_atom_numbers_list=manual_dihedral_atom_numbers_list,
            )
        else:
            raise ValueError(
                f"ERROR: The entered qm_engine = {qm_engine} and the only valid options are "
                f"{['gaussian', 'gaussian_style_final_files']}"
            )

        # check if QM and mol2 file elements match
        [
            atom_name_mol2_list,
            element_name_mol2_list,
        ] = get_atom_names_and_elements_from_mol2(mol2_file)
        if all_element_names_list != element_name_mol2_list:
            raise ValueError(
                f"ERROR: The QM elements do not match the mol2 file elements, in order. \n"
                f"This does not guarantee that the element postions are correct. \n"
                f"mol2 file element names = {element_name_mol2_list} \n"
                f"QM file element names = {all_element_names_list}"
            )

        qm_atom_types_from_psf_map = [
            atom_number_to_atom_type_dict[dihedral_atom_numbers_list[0]],
            atom_number_to_atom_type_dict[dihedral_atom_numbers_list[1]],
            atom_number_to_atom_type_dict[dihedral_atom_numbers_list[2]],
            atom_number_to_atom_type_dict[dihedral_atom_numbers_list[3]],
        ]

        for iter_u, dih_num_u in enumerate(
            matching_dihedral_types_by_atom_type_list
        ):
            if dih_num_u not in [
                fit_dihedral_atom_types,
                fit_dihedral_atom_types_reversed,
            ]:
                raise ValueError(
                    f"ERROR: The extracted identical atom dihedral types = {dih_num_u} "
                    f"are not equal to the fitted dihedral type "
                    f"('fit_dihedral_atom_types') = {fit_dihedral_atom_types}, "
                    f"or the reverse order of the fitted dihedral type "
                    f"('fit_dihedral_atom_types') = {fit_dihedral_atom_types_reversed}."
                )

        # check the atom numbers of the PSF against QM
        if qm_atom_types_from_psf_map not in [
            fit_dihedral_atom_types,
            fit_dihedral_atom_types_reversed,
        ]:
            raise ValueError(
                f"ERROR: When the QM ({qm_engine}) dihedral atoms = {qm_atom_types_from_psf_map}, "
                f"when they are mapped to atom types/classes via the PSF data, which "
                f"do not match fitted dihedral type "
                f"('fit_dihedral_atom_types') = {fit_dihedral_atom_types}, "
                f"or the reverse order of the fitted dihedral type "
                f"('fit_dihedral_atom_types') = {fit_dihedral_atom_types_reversed}."
            )

        # get the dihedral coordinates for each qm scan  via the atom numbers and the
        all_matching_dihedral_coordinates_angstroms_added_to_k_values_list = []
        all_matching_dihedral_phi_degrees_added_to_k_values_list = []
        all_sum_opls_const_1_plus_or_minus_cos_n_list = []
        for qm_scan_i, coords_qm_scan_i in enumerate(coordinates_ang_list):
            dih_coor_iter_list = []
            dih_phi_iter_list = []
            for dih_q, dih_list_4_atom_number_q in enumerate(
                matching_dihedral_types_by_atom_numbers_list
            ):
                atom_coor_iter_list = []
                for dih_atom_number_q in dih_list_4_atom_number_q:
                    # added - 1 to below because atom numbering starts at 1 and Python lists at 0.
                    atom_coor_iter_list.append(
                        coords_qm_scan_i[dih_atom_number_q - 1]
                    )

                # get the dihedral coordinates lists for each QM scan
                dih_coor_iter_list.append(atom_coor_iter_list)

                # get the dihedral phi agles lists for each QM scan

                dih_phi_iter_list.append(
                    mdf_math.dihedral_angle(
                        atom_coor_iter_list[0],
                        atom_coor_iter_list[1],
                        atom_coor_iter_list[2],
                        atom_coor_iter_list[3],
                    )
                )

            all_matching_dihedral_coordinates_angstroms_added_to_k_values_list.append(
                dih_coor_iter_list
            )
            all_matching_dihedral_phi_degrees_added_to_k_values_list.append(
                dih_phi_iter_list
            )
            all_sum_opls_const_1_plus_or_minus_cos_n_list.append(
                mdf_math.sum_opls_const_1_plus_or_minus_cos_n_values(
                    dih_phi_iter_list
                )
            )

    return [
        matching_dihedral_types_by_atom_numbers_list,
        matching_dihedral_types_by_atom_type_list,
        all_matching_dihedral_coordinates_angstroms_added_to_k_values_list,
        all_matching_dihedral_phi_degrees_added_to_k_values_list,
        all_sum_opls_const_1_plus_or_minus_cos_n_list,
    ]


def change_gomc_ff_file_dihedral_values(
    read_gomc_ff_filename,
    new_gomc_ff_filename,
    fit_dihedral_atom_types,
    fit_dihedral_opls_k_0_1_2_3_4_values=[0, 0, 0, 0, 0],
    zero_dihedral_atom_types=None,
):
    """Rewrite the GOMC/CHARMM style force field file with modified fit dihedral values and other zeroed dihedrals.

    This rewrites the the GOMC/CHARMM style force field file with modified/entered values
    for the dihedral being fit ('fit_dihedral_opls_k_0_1_2_3_4_values'). It also allows
    other dihedrals to be set to zero (via 'zero_dihedral_atom_types'), allowing
    proper dihedral fitting procedures, such as in COOH and amide groups.  This
    function replaces the fit dihedral ('fit_dihedral_atom_types') values in
    with their new set values ('fit_dihedral_opls_k_0_1_2_3_4_values'), and
    zeros all the other selected dihedrals ('zero_dihedral_atom_types').


    Parameters
    ----------
    read_gomc_ff_filename: str
        The GOMC FF file name that will be modified with the new settable dihedral values
        (via 'fit_dihedral_atom_types' and 'fit_dihedral_opls_k_0_1_2_3_4_values'), and other
        dihedrals which will be zeroed ('zero_dihedral_atom_types').
    fit_dihedral_atom_types: list of four (4) strings (Example: ['HC', 'CT, 'CT, 'HC'])
        The atom types/classes (strings in the list) of the dihedral which is
        being fitted with non-zero k-values.
    fit_dihedral_opls_k_0_1_2_3_4_values: list of five (5) floats or int (Example: [k0, k1, k2, k3, k4])
        The OPLS k-values (k0, k1, k2, k3, k4) that the 'fit_dihedral_atom_types' are changed too.
        NOTE: These need input in their correct unit form (i.e., kcal/mol or K) for the
        given GOMC force file outputs (i.e., kcal/mol for LJ and K for Mie or Exp6).
    zero_dihedral_atom_types: nest list with lists of four (4) strings, default=None
        The nests list(s) of the other dihedrals, that need to have their k-values zeroed to
        properly fit the the 'fit_dihedral_atom_types' dihedral.

        Example: [['CT', 'CT, 'CT, 'HC'], ['NT', 'CT, 'CT, 'HC']]

    Notes
    -------
    Write a modified GOMC/CHARMM style force field file
        Force files are written copied from the existing force field file
        (Example: MoSDeF-GOMC force file (.inp) file), and rewritten with the
        fit new dihedral values 'fit_dihedral_opls_k_0_1_2_3_4_values', and any other
        dihedrals that need to be zeroed out ('zero_dihedral_atom_types') for
        the 'fit_dihedral_atom_types' dihedral to be fit.
    """

    # check if the fitted dihedral is input correctly
    if (
        not isinstance(fit_dihedral_atom_types, list)
        or len(fit_dihedral_atom_types) != 4
        or not isinstance(fit_dihedral_atom_types[0], str)
        or not isinstance(fit_dihedral_atom_types[1], str)
        or not isinstance(fit_dihedral_atom_types[2], str)
        or not isinstance(fit_dihedral_atom_types[3], str)
    ):
        raise TypeError(
            f"ERROR: The input 'fit_dihedral_atom_types' variable = {fit_dihedral_atom_types}, "
            f"but it needs to be a list of 4 strings, "
            f"where the strings are the atom types/classes. Example: ['HC', 'CT', 'CT', 'HC']."
        )

    # check if the other dihedral which need zeroed are input correctly
    zero_dihedral_atom_types_error = (
        f"ERROR: The zero_dihedral_atom_types' variable need to be 'None, "
        f"a list of strings, "
        f"or a nested list containing 1 or more list(s) of 4 strings, "
        f"where the strings are the atom types/classes."
    )
    if isinstance(zero_dihedral_atom_types, list):
        if (
            len(zero_dihedral_atom_types) == 4
            and isinstance(zero_dihedral_atom_types[0], str)
            and isinstance(zero_dihedral_atom_types[1], str)
            and isinstance(zero_dihedral_atom_types[2], str)
            and isinstance(zero_dihedral_atom_types[3], str)
        ):
            # make it a nested list for the other files input
            zero_dihedral_atom_types = [zero_dihedral_atom_types]

        elif len(zero_dihedral_atom_types) >= 1 and isinstance(
            zero_dihedral_atom_types, list
        ):
            for list_i in zero_dihedral_atom_types:
                if not len(list_i) == 4:
                    raise TypeError(zero_dihedral_atom_types_error)

                    if (
                        isinstance(list_i[0], str)
                        and isinstance(list_i[1], str)
                        and isinstance(list_i[2], str)
                        and isinstance(list_i[3], str)
                    ):
                        raise TypeError(zero_dihedral_atom_types_error)

    elif not isinstance(zero_dihedral_atom_types, type(None)):
        raise TypeError(zero_dihedral_atom_types_error)

    # check fit_dihedral_opls_k_0_1_2_3_4_values if all are zero
    all_k_values_zero_fit_dihedral_atom_types_bool = True
    for k_i in fit_dihedral_opls_k_0_1_2_3_4_values:
        if k_i != 0:
            all_k_values_zero_fit_dihedral_atom_types_bool = False

    # ensure 'fit_dihedral_atom_types' are not in the 'zero_dihedral_atom_types'
    if zero_dihedral_atom_types is not None:
        if fit_dihedral_atom_types in zero_dihedral_atom_types:
            raise ValueError(
                "ERROR: The 'fit_dihedral_atom_types' can not also be in the 'zero_dihedral_atom_types'."
            )

    # create dictionary for fitted dihedral types and its status of being written to the ff file
    status_written_fit_dihedral_atom_types_dict = {
        str(fit_dihedral_atom_types): False
    }

    # create dictionary for fitted dihedral types to be zeroed and its status of being written to the ff file
    status_written_zero_dihedral_atom_types_dict = {}
    if zero_dihedral_atom_types is not None:
        for other_i in zero_dihedral_atom_types:
            status_written_zero_dihedral_atom_types_dict.update(
                {str(other_i): False}
            )

    # write the new GOMC force field file with the selected dihedrals zeroed out
    gomc_modified_kvalues_ff_file = open(f"{new_gomc_ff_filename}", "w")
    if zero_dihedral_atom_types is None:
        gomc_modified_kvalues_ff_file.write(
            f"* This file was modified from the original GOMC FF file, by zeroing out the "
            f"fitted dihedral = {fit_dihedral_atom_types}.\n"
        )
    else:
        gomc_modified_kvalues_ff_file.write(
            f"* This file was modified from the original GOMC FF file, by zeroing out the "
            f"fitted dihedral = {fit_dihedral_atom_types} "
            f"and the other selected dihedrals = {zero_dihedral_atom_types}.\n"
        )
    gomc_modified_kvalues_ff_file.write(
        f"* NOTE: The selected dihedrals may have been zeroed out in the original force "
        f"field file (XML file), but are rezeroed here also.\n"
    )

    dih_spacing = "{:10s} {:10s} {:10s} {:10s} {:15s} {:10s} {:15s} ! {:20s} {:20s} {:20s} {:20s}\n"
    number_sig_fig_for_ff_file = 10
    ff_inp_file = open(read_gomc_ff_filename, "r").readlines()
    for m, line_m in enumerate(ff_inp_file):
        split_line_m = line_m.split()

        # if the correct line(s) added the changed dihedral value to the set k-value (Kchi)
        if len(split_line_m) == 12:
            # get the forward and reversed dihedral, as it could be printed either way
            dih_m = [
                split_line_m[0],
                split_line_m[1],
                split_line_m[2],
                split_line_m[3],
            ]
            dih_m_reverse = [
                split_line_m[3],
                split_line_m[2],
                split_line_m[1],
                split_line_m[0],
            ]

            # set the fitted dihedral types to the 'fit_dihedral_opls_k_0_1_2_3_4_values'
            # in the GOMC/CHARMM style FF file
            if (
                (str(dih_m) or str(dih_m_reverse))
                in list(status_written_fit_dihedral_atom_types_dict.keys())
                and status_written_fit_dihedral_atom_types_dict[str(dih_m)]
                is False
                and (
                    dih_m == fit_dihedral_atom_types
                    or dih_m_reverse == fit_dihedral_atom_types
                )
            ):
                # get the opls dihedrals and convert them to CHARMM style periodic
                # periodic_dihedral_k_n_d_values =
                #   [[K0, n0, d0],
                #   [K1, n1, d1],
                #   [K2, n2, d2],
                #   [K3, n3, d3],
                #   [K4, n4, d4],
                #   [K5, n5, d5]]
                periodic_dihedral_k_n_d_values = OPLS_to_periodic(
                    fit_dihedral_opls_k_0_1_2_3_4_values[0],
                    fit_dihedral_opls_k_0_1_2_3_4_values[1],
                    fit_dihedral_opls_k_0_1_2_3_4_values[2],
                    fit_dihedral_opls_k_0_1_2_3_4_values[3],
                    fit_dihedral_opls_k_0_1_2_3_4_values[4],
                )

                if all_k_values_zero_fit_dihedral_atom_types_bool is True:
                    gomc_modified_kvalues_ff_file.write(
                        dih_spacing.format(
                            split_line_m[0],
                            split_line_m[1],
                            split_line_m[2],
                            split_line_m[3],
                            str(
                                mdf_math.round_to_sig_figs(
                                    periodic_dihedral_k_n_d_values[1][0],
                                    sig_figs=number_sig_fig_for_ff_file,
                                )
                            ),
                            str(int(periodic_dihedral_k_n_d_values[1][1])),
                            str(
                                mdf_math.round_to_sig_figs(
                                    periodic_dihedral_k_n_d_values[1][2],
                                    sig_figs=number_sig_fig_for_ff_file,
                                )
                            ),
                            split_line_m[8],
                            split_line_m[9],
                            split_line_m[10],
                            split_line_m[11],
                        )
                    )

                else:
                    if periodic_dihedral_k_n_d_values[0][0] != 0:
                        gomc_modified_kvalues_ff_file.write(
                            dih_spacing.format(
                                split_line_m[0],
                                split_line_m[1],
                                split_line_m[2],
                                split_line_m[3],
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[0][0],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                str(int(periodic_dihedral_k_n_d_values[0][1])),
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[0][2],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                split_line_m[8],
                                split_line_m[9],
                                split_line_m[10],
                                split_line_m[11],
                            )
                        )

                    if periodic_dihedral_k_n_d_values[1][0] != 0:
                        gomc_modified_kvalues_ff_file.write(
                            dih_spacing.format(
                                split_line_m[0],
                                split_line_m[1],
                                split_line_m[2],
                                split_line_m[3],
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[1][0],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                str(int(periodic_dihedral_k_n_d_values[1][1])),
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[1][2],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                split_line_m[8],
                                split_line_m[9],
                                split_line_m[10],
                                split_line_m[11],
                            )
                        )

                    if periodic_dihedral_k_n_d_values[2][0] != 0:
                        gomc_modified_kvalues_ff_file.write(
                            dih_spacing.format(
                                split_line_m[0],
                                split_line_m[1],
                                split_line_m[2],
                                split_line_m[3],
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[2][0],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                str(int(periodic_dihedral_k_n_d_values[2][1])),
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[2][2],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                split_line_m[8],
                                split_line_m[9],
                                split_line_m[10],
                                split_line_m[11],
                            )
                        )

                    if periodic_dihedral_k_n_d_values[3][0] != 0:
                        gomc_modified_kvalues_ff_file.write(
                            dih_spacing.format(
                                split_line_m[0],
                                split_line_m[1],
                                split_line_m[2],
                                split_line_m[3],
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[3][0],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                str(int(periodic_dihedral_k_n_d_values[3][1])),
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[3][2],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                split_line_m[8],
                                split_line_m[9],
                                split_line_m[10],
                                split_line_m[11],
                            )
                        )

                    if periodic_dihedral_k_n_d_values[4][0] != 0:
                        gomc_modified_kvalues_ff_file.write(
                            dih_spacing.format(
                                split_line_m[0],
                                split_line_m[1],
                                split_line_m[2],
                                split_line_m[3],
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[4][0],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                str(int(periodic_dihedral_k_n_d_values[4][1])),
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[4][2],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                split_line_m[8],
                                split_line_m[9],
                                split_line_m[10],
                                split_line_m[11],
                            )
                        )

                    if periodic_dihedral_k_n_d_values[5][0] != 0:
                        gomc_modified_kvalues_ff_file.write(
                            dih_spacing.format(
                                split_line_m[0],
                                split_line_m[1],
                                split_line_m[2],
                                split_line_m[3],
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[5][0],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                str(int(periodic_dihedral_k_n_d_values[5][1])),
                                str(
                                    mdf_math.round_to_sig_figs(
                                        periodic_dihedral_k_n_d_values[5][2],
                                        sig_figs=number_sig_fig_for_ff_file,
                                    )
                                ),
                                split_line_m[8],
                                split_line_m[9],
                                split_line_m[10],
                                split_line_m[11],
                            )
                        )

                    status_written_fit_dihedral_atom_types_dict[str(dih_m)] = (
                        True
                    )

            # set the zeroed dihedral types to zero in the GOMC/CHARMM style FF file
            elif (
                (str(dih_m) or str(dih_m_reverse))
                in list(status_written_zero_dihedral_atom_types_dict.keys())
                and zero_dihedral_atom_types is not None
                and status_written_zero_dihedral_atom_types_dict[str(dih_m)]
                is False
                and (
                    dih_m in zero_dihedral_atom_types
                    or dih_m_reverse in zero_dihedral_atom_types
                )
            ):
                # get the opls dihedrals and convert them to CHARMM style periodic
                # periodic_dihedral_k_n_d_values =
                #   [[K0, n0, d0],
                #   [K1, n1, d1],
                #   [K2, n2, d2],
                #   [K3, n3, d3],
                #   [K4, n4, d4],
                #   [K5, n5, d5]]
                gomc_modified_kvalues_ff_file.write(
                    dih_spacing.format(
                        split_line_m[0],
                        split_line_m[1],
                        split_line_m[2],
                        split_line_m[3],
                        str(
                            mdf_math.round_to_sig_figs(
                                0, sig_figs=number_sig_fig_for_ff_file
                            )
                        ),
                        str(int(1)),
                        str(
                            mdf_math.round_to_sig_figs(
                                180.0, sig_figs=number_sig_fig_for_ff_file
                            )
                        ),
                        split_line_m[8],
                        split_line_m[9],
                        split_line_m[10],
                        split_line_m[11],
                    )
                )

                status_written_zero_dihedral_atom_types_dict[str(dih_m)] = True

            elif (
                str(dih_m)
                not in list(status_written_fit_dihedral_atom_types_dict.keys())
                and str(dih_m_reverse)
                not in list(status_written_fit_dihedral_atom_types_dict.keys())
            ) and (
                str(dih_m)
                not in list(status_written_zero_dihedral_atom_types_dict.keys())
                and str(dih_m_reverse)
                not in list(status_written_zero_dihedral_atom_types_dict.keys())
            ):
                gomc_modified_kvalues_ff_file.write(f"{line_m}")

        else:
            gomc_modified_kvalues_ff_file.write(f"{line_m}")

    gomc_modified_kvalues_ff_file.close()
