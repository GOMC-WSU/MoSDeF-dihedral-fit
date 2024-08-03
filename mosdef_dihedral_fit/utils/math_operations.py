import math

import numpy as np


def round_to_sig_figs(value, sig_figs=3):
    """Round a number to the selected number of significant figures.

    Round a number to the selected number of significant figures.

    Parameters
    ----------
    value: int or float
        The number that will be rounded to the selected number of
        significant figures.
    sig_figs: int, default=3
        The number significant figures that the 'value' variable
        will be rounded too. If sig_fig=0, it will return 0.0.

    Returns
    -------
    value_rounded_to_sig_figs: float
        The input 'value' variable rounded to the selected number
        significant figures. If sig_fig=0, it will return 0.0.
    """
    if value == 0:
        value_rounded_to_sig_figs = 0

    else:
        value_rounded_to_sig_figs = float(
            round(value, sig_figs - int(math.floor(math.log10(abs(value)))) - 1)
        )

    return value_rounded_to_sig_figs


# **************************************************************
# dihedral angles calculations
# (START)
# **************************************************************


def normalize_vector(vector):
    """Generates the normalized vector from a given vector.

    Generates the normalized vector from a given vector.

    Parameters
    ----------
    vector: list or np.array, or floats or integers
        A vector in a list of np.array of any dimension.

    Returns
    -------
    normal_vector: ndarray, or floats or integers
        The normalized vector from the provided 'vector', with the
        same dimensions.
    """
    vector = np.array(vector)
    normalized_length = np.linalg.norm(vector)

    if normalized_length == 0:
        raise ValueError(
            f"ERROR: The normal vector = 0, indicating that these lines or planes "
            f"lay on top each other or are perpendicular. \n"
            f"The input vector = {vector}"
        )

    else:
        normal_vector = vector / normalized_length

    return np.array(normal_vector)


def angle_between_2_vectors(vector_1, vector_2):
    """Gets the angle between 2 vectors.

    Gets the angle between 2 vectors, which can be used to obtain the
    dihedral angles given their normal vectors to the planes.

    Parameters
    ----------
    vector_1: list or np.array, or floats or integers
        A vector in a list of np.array of any dimension.

    Returns
    -------
    angle_degrees: floats (in degrees)
        The angle (degrees), between the two (2) provided vectors (vector_1 and vector_2).
    """
    vector_1 = np.array(vector_1)
    vector_2 = np.array(vector_2)

    normal_vector_1 = normalize_vector(vector_1)
    normal_vector_2 = normalize_vector(vector_2)

    # get the required divisor 'normal_1_times_normal_2' and ensure it it not zero (0)
    if np.linalg.norm(vector_1) == 0 or np.linalg.norm(vector_2) == 0:
        raise ValueError(
            f"ERROR: The 'vector_1' or 'vector_2', and  |vector_1||vector_2| == 0, which means the \n"
            f"angle can not be calculated due to the divisor being zero in the formula; \n"
            f"angle = arccos[(vector_1 dot vector_2)/(|vector_1||vector_2|)]"
        )

    arc_cos_value = np.dot(normal_vector_1, normal_vector_2) / (
        np.linalg.norm(vector_1) * np.linalg.norm(vector_2)
    )

    # ensure arc_cos_value  is -1<= arc_cos_value <= 1
    # ... not 1.000000001 or -1.000000001 or it will yield 'nan'
    if np.isclose(arc_cos_value, 1.0):
        arc_cos_value = 1

    elif np.isclose(arc_cos_value, -1.0):
        arc_cos_value = -1

    angle_degrees = np.arccos(arc_cos_value) * 180 / np.pi

    return float(angle_degrees)


def dihedral_angle(
    atom_xyz_coord_1,
    atom_xyz_coord_2,
    atom_xyz_coord_3,
    atom_xyz_coord_4,
):
    """Gets the dihedral angle between the four (4) atom coordinates given in cartesian coordinates.

    Gets the dihedral angle between the four (4) atom coordinates given in cartesian coordinates.

    NOTE: the atoms must be entered in the proper order or the dihedral to obtain the
    correct results.

    Parameters
    ----------
    atom_xyz_coord_1: list of three (3) floats or ints
        The first (1st) atom/bead in the dihedral, in order.
    atom_xyz_coord_2: list of three (3) floats or ints
        The second (2nd) atom/bead in the dihedral, in order.
    atom_xyz_coord_3: list of three (3) floats or ints
        The third (3rd) atom/bead in the dihedral, in order.
    atom_xyz_coord_4: list of three (3) floats or ints
        The fourth/last (4th) atom/bead in the dihedral, in order.


    Returns
    -------
    dihedral_angle_degrees: float (in degrees)
        The dihedral angle, in degrees, between the four (4) atoms.
    """
    # check if any atom coordinates are the same
    if (
        list(atom_xyz_coord_1) == list(atom_xyz_coord_2)
        or list(atom_xyz_coord_1) == list(atom_xyz_coord_3)
        or list(atom_xyz_coord_1) == list(atom_xyz_coord_4)
        or list(atom_xyz_coord_2) == list(atom_xyz_coord_3)
        or list(atom_xyz_coord_2) == list(atom_xyz_coord_4)
        or list(atom_xyz_coord_3) == list(atom_xyz_coord_4)
    ):
        raise ValueError(
            f"ERROR: The one or more of the atom coordinates are the same when entered into \n"
            f"the 'get_dihedral_angle' functions.  In order, the entered atom coordinates: \n"
            f"atom_xyz_coord_1 = {atom_xyz_coord_1}; \n"
            f"atom_xyz_coord_2 = {atom_xyz_coord_2}; \n"
            f"atom_xyz_coord_3 = {atom_xyz_coord_3}; \n"
            f"atom_xyz_coord_4 = {atom_xyz_coord_4}. "
        )

    # push to array
    atom_xyz_coord_1 = np.array(atom_xyz_coord_1)
    atom_xyz_coord_2 = np.array(atom_xyz_coord_2)
    atom_xyz_coord_3 = np.array(atom_xyz_coord_3)
    atom_xyz_coord_4 = np.array(atom_xyz_coord_4)

    # get the vector between the bonded atoms/beads
    vector_between_atoms_1_2 = atom_xyz_coord_2 - atom_xyz_coord_1
    vector_between_atoms_2_3 = atom_xyz_coord_3 - atom_xyz_coord_2
    vector_between_atoms_3_4 = atom_xyz_coord_4 - atom_xyz_coord_3

    # get the normal vectors from the from the 3 bonded atoms/beads at both ends of the dihedral
    # (i.e., normal vector from the 1-2-3 atoms/beads and 2-3-4 atoms/bead planes
    normal_vector_for_atoms_1_2_3_plane = (
        # np.cross(vector_between_atoms_1_2, vector_between_atoms_2_3) /
        normalize_vector(
            np.cross(vector_between_atoms_1_2, vector_between_atoms_2_3)
        )
    )
    normal_vector_for_atoms_2_3_4_plane = (
        # np.cross(vector_between_atoms_2_3, vector_between_atoms_3_4) /
        normalize_vector(
            np.cross(vector_between_atoms_2_3, vector_between_atoms_3_4)
        )
    )

    dihedral_angle_degrees = float(
        angle_between_2_vectors(
            normal_vector_for_atoms_1_2_3_plane,
            normal_vector_for_atoms_2_3_4_plane,
        )
    )

    # correct the dihedral angle ('dihedral_angle_degrees') for location, rotation direction
    if (
        not np.dot(
            normal_vector_for_atoms_1_2_3_plane, vector_between_atoms_3_4
        )
        >= 0
    ):
        dihedral_angle_degrees = -dihedral_angle_degrees

    if dihedral_angle_degrees == 180:
        dihedral_angle_degrees = -180

    return dihedral_angle_degrees


# **************************************************************
# dihedral angles calculations
# (End)
# **************************************************************


def check_previous_qm_values_match(
    all_value_list, current_value, value_name, qm_engine, log_file_name
):
    """Checks if the QM log file read values match the last value in the appended list.

    Checks if the QM log file read values match the last value in the appended list,
    which were taken from previously read QM log files.

    Parameters
    ----------
    all_value_list: list of any type
        A list of the previously read QM log files values.
    current_value: any type
        The currently QM log file property value.
    value_name: str
        The name the property, which will be printed in the error output
        (Example: Number of atoms).
    qm_engine: str
        The name of the QM engine (Example: Gaussian)
    log_file_name: str
        The name of the log file being read.

    Returns
    -------
    list or ValueError:
        all_value_list; list
            An appended 'all_value_list' with the current value matches
            the last value in the list 'all_value_list'.
        ValueError:
            If the current value does not match the last value in the list 'all_value_list'.
    """
    if len(all_value_list) > 0:
        if current_value not in all_value_list:
            raise ValueError(
                f"ERROR: The {qm_engine} log file '{log_file_name}' does not have the same "
                f"{value_name} = {current_value}, "
                f"as other previous entries = {all_value_list[-1]}. "
                f"The molecule or property may not be the same in the multiple {qm_engine} log files, "
                f"but this may be be desired by the user to obtain a more general dihedral fit."
            )
    all_value_list.append(current_value)

    return all_value_list


def sum_opls_const_1_plus_or_minus_cos_n_values(phi_list):
    """Get OPLS (1 +/- cos(n * phi)) values for all k-values.

    This gets the OPLS (1 +/- cos(n * phi)) values for all k-values,
    which is required to fit the data, especially when multiple
    dihedrals of the same atom types/classes exist in the molecule
    that the dihedral is fit too.

    Example 0:
    const_1_minus_Cos_0_phi  = 1 , since k0 is a constant

    Example 1:
    const_1_plus_Cos_1_phi  = sum of the list using all phis [(1 + cos(1 * phi))] in k1 * (1 + cos(1 * phi))

    Example 2:
    const_1_minus_Cos_2_phi = sum of the list using all phis [(1 - cos(2 * phi))] in k2 * (1 - cos(2 * phi))

    Example 3:
    const_1_plus_Cos_3_phi  = sum of the list using all phis [(1 + cos(3 * phi))] in k3 * (1 + cos(3 * phi))

    Example 4:
    const_1_minus_Cos_4_phi = sum of the list using all phis [(1 - cos(4 * phi))] in k4 * (1 - cos(4 * phi))

    NOTE: These values are used to simplify the diheral fitting process,
    allowing the user to deal with four (4) OPLS constants (see below).


    Standard OPLS dihedral form
    .. math::
    opls_dihedral &= 1/2 *(
        &= k0
        &= + k1 * (1 + cos(1 * phi))
        &= + k2 * (1 - cos(2 * phi))
        &= + k3 * (1 + cos(3 * phi))
        &= + k4 * (1 - cos(4 * phi))
        &= )

    Modified OPLS dihedral form for all dihedrals with the same atom types/classes exist in the molecule.
    opls_dihedral_n_1 &= 1/2 *(
        &= k0
        &= + k1 * const_1_plus_Cos_1_phi
        &= + k2 * const_1_minus_Cos_2_phi
        &= + k3 * const_1_plus_Cos_3_phi
        &= + k4 * const_1_minus_Cos_4_phi
        &= )

    Parameters
    ----------
    phi_list: list
        A list of the dihedral angle phi values (degrees).

    Returns
    -------
    List of:
        const_1_minus_Cos_0_phi: float (unitless)
            This values is always zero, since k0 is not used this standard OPLS form.
        const_1_plus_Cos_1_phi: float (unitless)
            The sum of all phi values in the list [(1 - cos(2 * phi))] in k2 * (1 - cos(2 * phi))
        const_1_minus_Cos_2_phi: float (unitless)
            The sum of all phi values in the list [(1 - cos(2 * phi))] in k2 * (1 - cos(2 * phi))
        const_1_plus_Cos_3_phi: float (unitless)
            The sum of all phi values in the list [(1 + cos(3 * phi))] in k3 * (1 + cos(3 * phi))
        const_1_minus_Cos_4_phi: float (unitless)
            The sum of all phi values in the list [(1 - cos(4 * phi))] in k4 * (1 - cos(4 * phi))
    """
    const_1_minus_Cos_0_phi = 1
    const_1_plus_Cos_1_phi = 0
    const_1_minus_Cos_2_phi = 0
    const_1_plus_Cos_3_phi = 0
    const_1_minus_Cos_4_phi = 0

    for phi_i, phi_value_i in enumerate(phi_list):
        const_1_plus_Cos_1_phi += 1 + np.cos(1 * phi_value_i * np.pi / 180)
        const_1_minus_Cos_2_phi += 1 - np.cos(2 * phi_value_i * np.pi / 180)
        const_1_plus_Cos_3_phi += 1 + np.cos(3 * phi_value_i * np.pi / 180)
        const_1_minus_Cos_4_phi += 1 - np.cos(4 * phi_value_i * np.pi / 180)

    const_1_plus_or_minus_Cos_n_phi = [
        const_1_minus_Cos_0_phi,
        const_1_plus_Cos_1_phi,
        const_1_minus_Cos_2_phi,
        const_1_plus_Cos_3_phi,
        const_1_minus_Cos_4_phi,
    ]

    return const_1_plus_or_minus_Cos_n_phi


# opls dihedral function using for all combinations


def opls_dihedral(cos_powers_phi_and_constants_data, k0, k1, k2, k3, k4):
    """OPLS dihedral energy calculation with only the selected k-values.

    This is the OPLS dihedral energy calculation done with only the selected k-values.
    However, all k-values are to be input as this makes entering the
    k-values in a standard method for all the different OPLS dihderal
    configurations, and more importantly allows it to be fit to the data
    properly, which is not possible if all the k-values are included.
    The k-values are in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...),
    and the output ' dihedral_energy' energy are output in the same energy units.


    NOTE: THIS WILL FIT THE FORM IT IS FED, REGARDLESS IF IT IS THE CORRECT
    ANALYTICAL SOLUTION.  MEANING SOMETIMES YOU CAN NOT USE A K-VALUE BECAUSE
    OF MOLECULE SYMMETREY. THEREFORE, THIS IS MUST BE ACCOUNTED FOR OUTSIDE
    OF THIS FUNCTION.

    NOTE: ALL THE K-VALUE ENERGY UNITS MUST BE THE SAME.

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
    cos_powers_phi_and_constants_data: list or tuple of (cos_powers, scanned_phi, all_sum_opls_const_1_plus_or_minus_cos_n_list)
        cos_powers: str, int, or float,
            The str options are:   '1',  '1_3',  '2', '2_4', '1_2', '1_2_3', and '1_2_3_4'
            The int options are:   '1',   '13',  '2',  '24',  '12',   '123', and    '1234'
            The float options are: '1.', '13.', '2.', '24.', '12.',  '123.', and   '1234.'
            The type type of powers in the cosine series to use.

            Example 1:  '1' only set the k1 to a non-zero value
            Example 2:  '1_2_3' only set the k1, k2, and k3 to a non-zero value

            Example 3:  '1' only set the k1 to a non-zero value
            Example 4:  '123' only set the k1, k2, and k3 to a non-zero value

            Example 5:  '1' only set the k1 to a non-zero value
            Example 6:  '123.' only set the k1, k2, and k3 to a non-zero value

        scanned_phi: floats, int, nest list or tuple of floats/int
            The 'phi_data' angle of the dihedral is in degrees

            The floats, int options: the phi_data from a single point

            The nest list of floats/int: This nested list include the dihedral angle
                being rotated in QM, and all the other identical dihedral angles from
                the same atom typed/classed atoms in the test molecule/fragment.

        all_sum_opls_const_1_plus_or_minus_cos_n_list: list or tuple, default=None

            all_sum_opls_const_1_plus_or_minus_cos_n_list = [
            const_1_minus_Cos_0_phi,
            const_1_plus_Cos_1_phi,
            const_1_minus_Cos_2_phi,
            const_1_plus_Cos_3_phi,
            const_1_minus_Cos_4_phi
            ]

            A list of the OPLS k0 data, which is always 0 or 1 in this case.
            If all 0 values --> (0, 0, .., 0), then k0=0.
            If all 1 values --> (1, 1, .., 1), then 0=constant

            # It is critical that 'const_1_minus_Cos_0_phi' be all (0, 0, .., 0) for k0=0
            # and all (1, 1, .., 1) 'const_1_minus_Cos_0_phi' for k0=constant

            const_1_plus_Cos_1_phi:  values for all k-values, which
            is required to fit the data, especially when multiple dihedrals of the same
            atom types/classes exist in the molecule that the dihedral is fit too.

            Example 0:
            const_1_minus_Cos_0_phi  = k0

            Example 1:
            const_1_plus_Cos_1_phi  = sum of all phi values in the list
            [(1 + cos(1 * phi))] in k1 * (1 + cos(1 * phi))

            Example 2:
            const_1_minus_Cos_2_phi = sum of all phi values in the list
            [(1 - cos(2 * phi))] in k2 * (1 - cos(2 * phi))

            Example 3:
            const_1_plus_Cos_3_phi  = sum of all phi values in the list
            [(1 + cos(3 * phi))] in k3 * (1 + cos(3 * phi))

            Example 4:
            const_1_minus_Cos_4_phi = sum of all phi values in the list
            [(1 - cos(4 * phi))] in k4 * (1 - cos(4 * phi))

            The list length is the number of dihedrals which match the QM scanned atom/bead types,
            with the next nested list being the 'all_sum_opls_const_1_plus_or_minus_cos_n_list'
            as described.

    k0: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k0' value is the k-value for the opls dihedral where n=0,
        or the constant without a cosine multiple.
        NOTE: In this case, it is set to zero (0) regardless of the
        user entered value, because it is not in the equation form.
    k1: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k1' value is the k-value for the opls dihedral where n=1
        in the cosine multiple.
    k2: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k2' value is the k-value for the opls dihedral where n=2
        in the cosine multiple.
        NOTE: In this case, it is set to zero (0) regardless of the
        user entered value, because it is not in the equation form.
    k3: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k3' value is the k-value for the opls dihedral where n=3
        in the cosine multiple.
        NOTE: In this case, it is set to zero (0) regardless of the
        user entered value, because it is not in the equation form.
    k4: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k4' value is the k-value for the opls dihedral where n=4
        in the cosine multiple.
        NOTE: In this case, it is set to zero (0) regardless of the
        user entered value, because it is not in the equation form.

    Returns
    -------
    dihedral_energy; float or list of floats, same as the k-value energy units
        The OPLS dihedral energy from the valid phi_data and k-values
        (i.e., k1 in this case).
    """
    (
        cos_powers,
        scanned_phi,
        const_1_minus_Cos_0_phi_data,
        const_1_plus_Cos_1_phi_data,
        const_1_minus_Cos_2_phi_data,
        const_1_plus_Cos_3_phi_data,
        const_1_minus_Cos_4_phi_data,
    ) = cos_powers_phi_and_constants_data

    # check if all_sum_opls_const_1_plus_or_minus_cos_n_list is correct size:
    if (
        not isinstance(
            const_1_minus_Cos_0_phi_data,
            (list, np.ndarray, int, float, type(None)),
        )
        or not isinstance(
            const_1_minus_Cos_0_phi_data,
            (list, np.ndarray, int, float, type(None)),
        )
        or not isinstance(
            const_1_plus_Cos_1_phi_data,
            (list, np.ndarray, int, float, type(None)),
        )
        or not isinstance(
            const_1_minus_Cos_2_phi_data,
            (list, np.ndarray, int, float, type(None)),
        )
        or not isinstance(
            const_1_plus_Cos_3_phi_data,
            (list, np.ndarray, int, float, type(None)),
        )
        or not isinstance(
            const_1_minus_Cos_4_phi_data,
            (list, np.ndarray, int, float, type(None)),
        )
    ):
        raise TypeError(
            "ERROR: the 'cos_powers_phi_and_constants_data' values ("
            "const_1_minus_Cos_0_phi_data,"
            "const_1_plus_Cos_1_phi_data,"
            " const_1_minus_Cos_2_phi_data,"
            "const_1_plus_Cos_3_phi_data, and const_1_minus_Cos_4_phi_data"
            "), "
            f"are not all a (list, np.ndarray, int, float, or None) -> "
            f"const_1_minus_Cos_0_phi_data = {const_1_minus_Cos_0_phi_data}, \n"
            f"const_1_plus_Cos_1_phi_data = {const_1_plus_Cos_1_phi_data}, \n"
            f"const_1_minus_Cos_2_phi_data = {const_1_minus_Cos_2_phi_data}, \n"
            f"const_1_plus_Cos_3_phi_data = {const_1_plus_Cos_3_phi_data}, \n"
            f"const_1_minus_Cos_4_phi_data = {const_1_minus_Cos_4_phi_data}."
        )

    if (
        const_1_minus_Cos_0_phi_data is None
        or const_1_plus_Cos_1_phi_data is None
        or const_1_minus_Cos_2_phi_data is None
        or const_1_plus_Cos_3_phi_data is None
        or const_1_minus_Cos_4_phi_data is None
    ) and (
        const_1_minus_Cos_0_phi_data
        != const_1_plus_Cos_1_phi_data
        != const_1_minus_Cos_2_phi_data
        != const_1_plus_Cos_3_phi_data
        != const_1_minus_Cos_4_phi_data
    ):
        raise TypeError(
            "ERROR: If any of these values, "
            "const_1_minus_Cos_0_phi_data, "
            "const_1_plus_Cos_1_phi_data, "
            "const_1_minus_Cos_2_phi_data, "
            "const_1_plus_Cos_3_phi_data, "
            "const_1_minus_Cos_4_phi_data, "
            "are 'None', they all must be None"
        )

    else:
        if const_1_minus_Cos_0_phi_data is None:
            use_const_1_plus_minus_Cos_x_values_bool = False

        else:
            use_const_1_plus_minus_Cos_x_values_bool = True

    # check for acceptable cos_powers
    if isinstance(cos_powers, (str, int, float)):
        cos_powers_idential_value = cos_powers

        if not isinstance(scanned_phi, (int, float)):
            raise TypeError(
                f"ERROR: the 'cos_powers_phi_and_constants_data' values (cos_powers, scanned_phi), "
                f"are not the both a single value (str, int, or floats) -> "
                f"cos_powers = {cos_powers}, \n"
                f"scanned_phi = {scanned_phi}."
            )

    elif isinstance(cos_powers, (list, np.ndarray)):
        if not isinstance(scanned_phi, (list, np.ndarray)):
            raise TypeError(
                "ERROR: the 'cos_powers_phi_and_constants_data' values (cos_powers, scanned_phi), "
                f"are not the both a not a  (list or np.ndarray) -> "
                f"cos_powers = {cos_powers}, \n"
                f"scanned_phi = {scanned_phi}."
            )
        else:
            if (
                len(cos_powers)
                != len(scanned_phi)
                != len(const_1_minus_Cos_0_phi_data)
                != len(const_1_plus_Cos_1_phi_data)
                != len(const_1_minus_Cos_2_phi_data)
                != len(const_1_plus_Cos_3_phi_data)
                != len(const_1_minus_Cos_4_phi_data)
            ):
                raise ValueError(
                    f"ERROR: the 'cos_powers_phi_and_constants_data' values (cos_powers, phi_data), "
                    f"are not the same length -> "
                    f"len(cos_powers) = {len(cos_powers)}, \n"
                    f"len(scanned_phi) = {len(scanned_phi)}, \n"
                    f"len(const_1_minus_Cos_0_phi_data) = {len(const_1_minus_Cos_0_phi_data)}, \n"
                    f"len(const_1_plus_Cos_1_phi_data)  = {len(const_1_plus_Cos_1_phi_data) }, \n"
                    f"len(const_1_minus_Cos_2_phi_data) = {len(const_1_minus_Cos_2_phi_data)}, \n"
                    f"len(onst_1_plus_Cos_3_phi_data) = {len(const_1_plus_Cos_3_phi_data)}, \n"
                    f"len(const_1_minus_Cos_4_phi_data) = {len(const_1_minus_Cos_4_phi_data)}."
                )

            else:
                for z in range(0, len(cos_powers)):
                    if cos_powers[0] != cos_powers[z]:
                        raise ValueError(
                            "ERROR: the 'cos_powers' are not all the same for a given fit."
                        )

                cos_powers_idential_value = cos_powers[0]

    else:
        raise TypeError(
            f"ERROR: The cos_powers variable is a {type(cos_powers)}, but needs to be a "
            f"str, int, float, list, or np.ndarray."
        )

    if isinstance(cos_powers_idential_value, str):
        cos_powers_modified = ""
        for p in range(0, len(cos_powers_idential_value)):
            if cos_powers[p] == "-":
                cos_powers_modified += "_"
            else:
                cos_powers_modified += cos_powers_idential_value[p]

    else:
        cos_powers_modified = cos_powers_idential_value

    # check for acceptable cos_powers
    if cos_powers_modified not in [
        "1",
        "2",
        "3",
        "4",
        "1_3",
        "2_4",
        "1_2",
        "3_4",
        "1_2_3",
        "1_2_3_4",
        "1",
        "2",
        "3",
        "4",
        "1_3",
        "2_4",
        "1_2",
        "3_4",
        "1_2_3",
        "1_2_3_4",
        1,
        2,
        3,
        4,
        13,
        24,
        12,
        34,
        123,
        1234,
        1.0,
        2.0,
        3.0,
        4.0,
        13.0,
        24.0,
        12.0,
        34.0,
        123.0,
        1234.0,
    ]:
        raise ValueError(
            f"ERROR: {cos_powers} was entered for the cos_powers variable, but the only "
            f"available options are "
            f"["
            f"{'1', '2', '3', '4' '1_3', '2_4', '1_2', '3_4', '1_2_3', '1_2_3_4', } "
            f"'1, 2, 3, 4, 13, 24, 12, 34, 123, 1234, "
            f" 1., 2., 3., 4., 13., 24., 12., 34., 123., 1234.]"
        )

    # check if using k0=0 ('opls_force_k0_zero'=True) by summing
    # 'sorted_const_1_minus_Cos_0_phi_data_lists'
    # It is critical that 'const_1_minus_Cos_0_phi_data' be all (0, 0, .., 0) for k0=0
    # and all (1, 1, .., 1) 'const_1_minus_Cos_0_phi_data' for k0=constant
    if (
        const_1_minus_Cos_0_phi_data is not None
        and sum(const_1_minus_Cos_0_phi_data) == 0
        # or const_1_minus_Cos_0_phi_data is None
    ):
        k0 = 0

    if cos_powers_modified in ["1", 1, 1.0]:
        # NOTE: THE ALL BUT THE 'kx' VALUES ARE REPLACES WITH A CONSTANT,

        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1 / 2 * (k0 + k1 * (1 + np.cos(1 * scanned_phi * np.pi / 180)))
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = 1 / 2 * (k0 + k1 * const_1_plus_Cos_1_phi_data)

    elif cos_powers_modified in ["2", 2, 2.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1 / 2 * (k0 + k2 * (1 - np.cos(2 * scanned_phi * np.pi / 180)))
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = 1 / 2 * (k0 + k2 * const_1_minus_Cos_2_phi_data)

    elif cos_powers_modified in ["3", 3, 3.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1 / 2 * (k0 + k3 * (1 + np.cos(3 * scanned_phi * np.pi / 180)))
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:

            dihedral_energy = 1 / 2 * (k0 + k3 * const_1_plus_Cos_3_phi_data)

    elif cos_powers_modified in ["4", 4, 4.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1 / 2 * (k0 + k4 * (1 - np.cos(4 * scanned_phi * np.pi / 180)))
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = 1 / 2 * (k0 + k4 * const_1_minus_Cos_4_phi_data)

    elif cos_powers_modified in ["1_3", 13, 13.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + k1 * (1 + np.cos(1 * scanned_phi * np.pi / 180))
                    + k3 * (1 + np.cos(3 * scanned_phi * np.pi / 180))
                )
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k1 * const_1_plus_Cos_1_phi_data
                    + k3 * const_1_plus_Cos_3_phi_data
                )
            )

    elif cos_powers_modified in ["2_4", 24, 24.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k2 * (1 - np.cos(2 * scanned_phi * np.pi / 180))
                    + k4 * (1 - np.cos(4 * scanned_phi * np.pi / 180))
                )
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k2 * const_1_minus_Cos_2_phi_data
                    + k4 * const_1_minus_Cos_4_phi_data
                )
            )

    elif cos_powers_modified in ["3_4", 34, 34.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k3 * (1 + np.cos(3 * scanned_phi * np.pi / 180))
                    + k4 * (1 - np.cos(4 * scanned_phi * np.pi / 180))
                )
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k3 * const_1_plus_Cos_3_phi_data
                    + k4 * const_1_minus_Cos_4_phi_data
                )
            )

    elif cos_powers_modified in ["1_2", 12, 12.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + k1 * (1 + np.cos(1 * scanned_phi * np.pi / 180))
                    + k2 * (1 - np.cos(2 * scanned_phi * np.pi / 180))
                )
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = (
                1
                / 2
                * (
                    k0 * const_1_minus_Cos_0_phi_data
                    + k1 * const_1_plus_Cos_1_phi_data
                    + k2 * const_1_minus_Cos_2_phi_data
                )
            )

    elif cos_powers_modified in ["1_2_3", 123, 123.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k1 * (1 + np.cos(1 * scanned_phi * np.pi / 180))
                    + k2 * (1 - np.cos(2 * scanned_phi * np.pi / 180))
                    + k3 * (1 + np.cos(3 * scanned_phi * np.pi / 180))
                )
            )
        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k1 * const_1_plus_Cos_1_phi_data
                    + k2 * const_1_minus_Cos_2_phi_data
                    + k3 * const_1_plus_Cos_3_phi_data
                )
            )

    elif cos_powers_modified in ["1_2_3_4", 1234, 1234.0]:
        if use_const_1_plus_minus_Cos_x_values_bool is False:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k1 * (1 + np.cos(1 * scanned_phi * np.pi / 180))
                    + k2 * (1 - np.cos(2 * scanned_phi * np.pi / 180))
                    + k3 * (1 + np.cos(3 * scanned_phi * np.pi / 180))
                    + k4 * (1 - np.cos(4 * scanned_phi * np.pi / 180))
                )
            )

        elif use_const_1_plus_minus_Cos_x_values_bool is True:
            dihedral_energy = (
                1
                / 2
                * (
                    k0
                    + +k1 * const_1_plus_Cos_1_phi_data
                    + k2 * const_1_minus_Cos_2_phi_data
                    + k3 * const_1_plus_Cos_3_phi_data
                    + k4 * const_1_minus_Cos_4_phi_data
                )
            )

    return dihedral_energy


# periodic dihedral function using
def periodic_dihedral_n_1_2_3_4_5(
    phi_data,
    K_0,
    K_1,
    K_2,
    K_3,
    K_4,
    K_5,
    n_0,
    n_1,
    n_2,
    n_3,
    n_4,
    n_5,
    d_0,
    d_1,
    d_2,
    d_3,
    d_4,
    d_5,
):
    """Periodic or CHARMM style dihedral energy calculation from n=1 to n=5.

    This is the Periodic or CHARMM style dihedral energy calculation.
    The K_x-values are in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...),
    and the output ' dihedral_energy' energy are output in the same energy units.
    The 'n_x' values is the n-value cosine function scalar value or
    the cosine power representation, where n_x and d_x match the K_x x-values.
    The 'd_x' values is the d-value cosine phase angle or phase shift.

    NOTE: ALL THE K_x-VALUE ENERGY UNITS MUST BE THE SAME.

    NOTE: All the K_x, n_x, and d_x x-values are a pair (i.e., K_1, n_1, and d_1)

    .. math::
        periodic_dihedral &= K_0 * (1 + cos(n_0*t - d_0)) + \\
                          &= K_1 * (1 + cos(n_1*t - d_1)) + \\
                          &= K_2 * (1 + cos(n_2*t - d_2)) + \\
                          &= K_3 * (1 + cos(n_3*t - d_3)) + \\
                          &= K_4 * (1 + cos(n_4*t - d_4)) + \\
                          &= K_5 * (1 + cos(n_5*t - d_5))

    Parameters
    ----------
    phi_data: float, int, or list of floats/int, in degrees
        The 'phi_data' angle of the dihedral is in degrees
    K_0, K_1, K_2, K_3, K_4, and K_5: int or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'K_x' values is the k-value constant scalar for the
        Periodic or CHARMM style dihedral, where n_x and d_x match
        the K_x x-values.
    n_0, n_1, n_2, n_3, n_4, and n_5: int, or float
        The 'n_x' values is the n-value cosine function scalar value or
        the cosine power representation, where n_x and d_x match
        the K_x x-values.
    d_0, d_1, d_2, d_3, d_4, and d_5: int, or float in degrees
        The 'd_x' values is the d-value cosine phase angle or phase shift,
        where n_x and d_x match the K_x x-values.

    Returns
    -------
    dihedral_energy; float or list of floats, same as the k-value energy units
        The Periodic or CHARMM style dihedral energy from the entered K_x, n_x, and d_x values.
    """
    # phi_data in degrees
    dihedral_energy = (
        K_0 * (1 + np.cos((n_0 * phi_data - d_0) * np.pi / 180))
        + K_1 * (1 + np.cos((n_1 * phi_data - d_1) * np.pi / 180))
        + K_2 * (1 + np.cos((n_2 * phi_data - d_2) * np.pi / 180))
        + K_3 * (1 + np.cos((n_3 * phi_data - d_3) * np.pi / 180))
        + K_4 * (1 + np.cos((n_4 * phi_data - d_4) * np.pi / 180))
        + K_5 * (1 + np.cos((n_5 * phi_data - d_5) * np.pi / 180))
    )

    return dihedral_energy


# RB torsion function using
def RB_torsion_n_1_2_3_4_5(phi_data, k_0, k_1, k_2, k_3, k_4, k_5):
    """Ryckaert-Bellemans (RB) torsion energy calculation from n=1 to n=5.

    This is the Ryckaert-Bellemans (RB) torsion style energy calculation.
    The K_x-values are in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...),
    and the output ' dihedral_energy' energy are output in the same energy units.

    NOTE: ALL THE k_x-VALUE ENERGY UNITS MUST BE THE SAME.

    .. math::
        RB_torsions &= k_0 + k_1*cos(psi) + k_2*cos(psi)^2 + k_3*cos(psi)^3 + \\
                    &= k_4*cos(psi)^4 + k_5*cos(psi)^5

    Parameters
    ----------
    phi_data: float, int, or list of floats/int, in degrees
        The 'phi_data' angle of the dihedral is in degrees
    k_0, k_1, k_2, k_3, k_4, and k_5: int or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k_x' values is the k-value constant scalar for the Ryckaert-Bellemans (RB) torsion style.

    Returns
    -------
    torsion_energy; float or list of floats, same as the k-value energy units
        The Ryckaert-Bellemans (RB) torsion style energy from the entered k_x values.
    """
    # phi_data in degrees
    # psi_data = phi_data - Pi
    psi_data = phi_data * np.pi / 180 - np.pi

    torsion_energy = (
        k_0
        + k_1 * np.cos(psi_data)
        + k_2 * np.cos(psi_data) ** 2
        + k_3 * np.cos(psi_data) ** 3
        + k_4 * np.cos(psi_data) ** 4
        + k_5 * np.cos(psi_data) ** 5
    )

    return torsion_energy


# OPLS dihedral energy only function using
def opls_dihedral_n_1_2_3_4(phi_data, k_0, k_1, k_2, k_3, k_4):
    """OPLS dihedral energy calculation with only the selected k-values.

    This is the OPLS dihedral energy calculation done with only the selected k-values.
    However, all k-values are to be input as this makes entering the
    k-values in a standard method for all the different OPLS dihderal
    configurations, and more importantly allows it to be fit to the data
    properly, which is not possible if all the k-values are included.
    The k-values are in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...),
    and the output ' dihedral_energy' energy are output in the same energy units.

    NOTE: ALL THE K-VALUE ENERGY UNITS MUST BE THE SAME.

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
    phi_data: float, int, or list of floats/int, in degrees
        The 'phi_data' angle of the dihedral is in degrees
    k_0: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k0' value is the k-value for the opls dihedral where n=0,
        or the constant without a cosine multiple.
        NOTE: In this case, it is set to zero (0) regardless of the
        user entered value, because it is not in the equation form.
    k_1: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k1' value is the k-value for the opls dihedral where n=1
        in the cosine multiple.
    k_2: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k2' value is the k-value for the opls dihedral where n=2
        in the cosine multiple.
        NOTE: In this case, it is set to zero (0) regardless of the
        user entered value, because it is not in the equation form.
    k_3: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k3' value is the k-value for the opls dihedral where n=3
        in the cosine multiple.
        NOTE: In this case, it is set to zero (0) regardless of the
        user entered value, because it is not in the equation form.
    k_4: int, or float, in energy units (i.e., kcal/mol, kJ/mol, Kelvin, ...)
        The 'k4' value is the k-value for the opls dihedral where n=4
        in the cosine multiple.
        NOTE: In this case, it is set to zero (0) regardless of the
        user entered value, because it is not in the equation form.

    Returns
    -------
    dihedral_energy; float or list of floats, same as the k-value energy units
        The OPLS dihedral energy from the valid phi_data and k-values
        (i.e., k1 in this case).
    """
    dihedral_energy = (
        1
        / 2
        * (
            k_0
            + k_1 * (1 + np.cos(1 * (phi_data * np.pi / 180)))
            + k_2 * (1 - np.cos(2 * (phi_data * np.pi / 180)))
            + k_3 * (1 + np.cos(3 * (phi_data * np.pi / 180)))
            + k_4 * (1 - np.cos(4 * (phi_data * np.pi / 180)))
        )
    )

    return dihedral_energy


def get_r_squared(data_points, fitted_values):
    """Get the R**2 (R squared) values from the fitted equation.

    Get the R**2 (R squared) values from the fitted equation, which is
    utilized to determine how good a fit is.


    Parameters
    ----------
    data_points: list, tuple, or numpy.array of floats or integers (same length as 'fitted_values')
        The data points, simulated, or experimental values that is utilized
        to fit a given equation.
    fitted_values: list, tuple, or numpy.array of floats or integers (same length as 'data_points')
        The values that are calculated from the equation that was derived
        from the 'data_points' list.

    Returns
    -------
    r_squared: float
        The R**2 (R-squared) value determines how good the equation fit is
        compared to the data that was used in the fitting process.
        These values are typically 0-1, but can be negative if the intercept
        constant or other parameters are manually not used in the fitting process.
    """
    if not isinstance(
        data_points,
        (list, type((1, 2)), type(np.array([1])), type(np.ndarray([1]))),
    ) or not isinstance(
        fitted_values,
        (list, type((1, 2)), type(np.array([1])), type(np.array([1]))),
    ):
        raise TypeError(
            f"ERROR: Both the 'data_points' and 'fitted_values' must be "
            f"list, tuple, or numpy.array of the same length. The are the following: \n"
            f"- type(data_points) = {type(data_points)} \n"
            f"- len(data_points) = {len(data_points)} \n"
            f"- type(fitted_value) = {type(fitted_values)} \n"
            f"- len(fitted_value) = {len(fitted_values)}"
        )

    for r_i in range(0, len(data_points)):
        if not isinstance(data_points[r_i], (float, int)) or not isinstance(
            fitted_values[r_i], (float, int)
        ):
            raise TypeError(
                f"ERROR: The 'data_points' or 'fitted_values' lists do not contain all floats or integers."
            )

    rss = 0
    tss = 0
    for iter_i, value_i in enumerate(data_points):
        rss += (data_points[iter_i] - fitted_values[iter_i]) ** 2
        tss += (data_points[iter_i] - np.average(data_points)) ** 2

    r_squared_0_to_1 = 1 - rss / tss

    return r_squared_0_to_1
