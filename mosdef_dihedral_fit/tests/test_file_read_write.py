import os

import unyt as u

import mosdef_dihedral_fit.utils.math_operations as mdf_math
from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import (
    fit_dihedral_with_gomc,
)
from mosdef_dihedral_fit.tests.base_test import BaseTest
from mosdef_dihedral_fit.utils.file_read_and_write import (
    change_gomc_ff_file_dihedral_values,
    check_guassian_angle_energy_file_correct,
    check_guassian_optimized_coordinate_file_correct,
    get_atom_names_and_elements_from_mol2,
    get_atom_names_and_elements_from_pdb,
    get_final_gaussian_output_file_data,
    get_gaussian_log_file_data,
    get_matching_dihedral_info_and_opls_fitting_data,
    write_qm_data_files,
    write_restart_coor_from_xyz_file,
    write_xyz_file_from_gaussian_coordinates,
)
from mosdef_dihedral_fit.utils.io import get_mosdef_dihedral_fit_fn

# user changable variable, as it needs to be run locally
# gomc_binary_directory = "/Users/brad/Programs/GOMC/GOMC_2_75/bin"
gomc_binary_directory = "/Users/calcraven/Documents/Vanderbilt/Research/MoSDeF/Dihedral_Fitter/GOMC/bin"


class TestFileReading(BaseTest):
    def test_get_from_mol2(self):
        fn = "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
        full_fn = get_mosdef_dihedral_fit_fn(fn)
        (
            atom_namesList,
            element_namesList,
        ) = get_atom_names_and_elements_from_mol2(full_fn)

        assert atom_namesList == [
            "C1",
            "C2",
            "C3",
            "O1",
            "O2",
            "H1",
            "H2",
            "H3",
            "H4",
            "H5",
            "H6",
        ]
        assert element_namesList == [
            "C",
            "C",
            "C",
            "O",
            "O",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
        ]

    def test_get_from_pdb(self):
        fn = "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.pdb"
        full_fn = get_mosdef_dihedral_fit_fn(fn)
        (
            atom_namesList,
            element_namesList,
        ) = get_atom_names_and_elements_from_pdb(full_fn)

        assert atom_namesList == [
            "C1",
            "C2",
            "C3",
            "O1",
            "O2",
            "H1",
            "H2",
            "H3",
            "H4",
            "H5",
            "H6",
        ]
        assert element_namesList == [
            "C",
            "C",
            "C",
            "O",
            "O",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
        ]

    def test_write_xyz_file_from_gaussian(self):
        atom_namesList = [
            "C1",
            "C2",
            "C3",
            "O1",
            "O2",
            "H1",
            "H2",
            "H3",
            "H4",
            "H5",
            "H6",
        ]
        fn = "gaussian_style_output_files/CT_CT_C_OH/output/"
        full_fn = get_mosdef_dihedral_fit_fn(fn) + "dihedral_coords_position_"
        extension = ".txt"

        write_xyz_file_from_gaussian_coordinates(
            atom_namesList, full_fn, extension, "./", 2
        )
        assert os.listdir() == [
            "dihedral_coords_position_1.xyz",
            "dihedral_coords_position_2.xyz",
        ]
        write_restart_coor_from_xyz_file("./", 2)
        assert "dihedral_coords_position_1.coor" in os.listdir()
        assert "dihedral_coords_position_2.coor" in os.listdir()

    def test_check_guassian_angle_energy_file_correct(self):
        # load from tests/files
        pass

    def test_check_guassian_optimized_coordinate_file_correct(self):
        # load from tests/files
        pass

    def test_get_final_gaussian_output_file_data(self):
        # load from tests/files
        pass

    def test_get_gaussian_log_file_data(self):
        # load from tests/files
        pass

    def test_write_qm_data_files(self):
        # load from tests/files
        pass

    def test_get_matching_dihedral_info_and_opls_fitting_data(self):
        # load from tests/files
        pass

    def test_change_gomc_ff_file_dihedral_values(self):
        # load from tests/files
        pass
