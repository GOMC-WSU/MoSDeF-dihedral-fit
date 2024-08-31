import os

import numpy as np
import pytest

import mosdef_dihedral_fit.utils.math_operations as mdf_math
from mosdef_dihedral_fit.tests.base_test import BaseTest
from mosdef_dihedral_fit.utils.file_read_and_write import (
    change_gomc_ff_file_dihedral_values,
    check_gaussian_angle_energy_file_correct,
    check_gaussian_optimized_coordinate_file_correct,
    get_atom_names_and_elements_from_mol2,
    get_atom_names_and_elements_from_pdb,
    get_final_gaussian_output_file_data,
    get_gaussian_log_file_data,
    get_matching_dihedral_info_and_opls_fitting_data,
    write_qm_data_files,
    write_restart_coor_from_xyz_file,
    write_xyz_file_from_gaussian_coordinates,
)

# user changable variable, as it needs to be run locally
# gomc_binary_directory = "/Users/brad/Programs/GOMC/GOMC_2_75/bin"
gomc_binary_directory = "/Users/calcraven/Documents/Vanderbilt/Research/MoSDeF/Dihedral_Fitter/GOMC/bin"


class TestFileReading(BaseTest):
    def test_get_from_mol2(self):
        fn = "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
        full_fn = self.get_fn(fn)
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

        error_fn = "with_errors/CT_CT_C_3_OH_with_errors.mol2"
        full_error_fn = self.get_fn(error_fn)
        error_msg = (
            "ERROR: The provided mol2 format is not the required VMD TRIPOS format "
            "or the mol2 file has zero \(0\) atoms or beads in it."
        )
        with pytest.raises(TypeError, match=error_msg):
            get_atom_names_and_elements_from_mol2(full_error_fn)

    def test_get_from_pdb(self):
        fn = "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.pdb"
        full_fn = self.get_fn(fn)
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
        error_fn = "with_errors/CT_CT_C_3_OH_with_errors.pdb"
        full_error_fn = self.get_fn(error_fn)
        error_msg = (
            "ERROR: The provided pdb format is not the required PDB format "
            "or the pdb file has zero \(0\) atoms or beads in it."
        )
        with pytest.raises(TypeError, match=error_msg):
            get_atom_names_and_elements_from_pdb(full_error_fn)

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
        full_fn = os.path.join(self.get_fn(fn), "dihedral_coords_position_")
        extension = ".txt"

        write_xyz_file_from_gaussian_coordinates(
            atom_namesList, full_fn, extension, "./", 2
        )
        assert "dihedral_coords_position_1.xyz" in os.listdir()
        assert "dihedral_coords_position_2.xyz" in os.listdir()

        error_msg = "ERROR: The 'qm_coordinate_file_extension' variable extension is not listed or a an empty string."
        with pytest.raises(ValueError, match=error_msg):
            write_xyz_file_from_gaussian_coordinates(
                atom_namesList, full_fn, "", "./", 2
            )
        write_restart_coor_from_xyz_file("./", 2)
        assert "dihedral_coords_position_1.coor" in os.listdir()
        assert "dihedral_coords_position_2.coor" in os.listdir()

        error_fn = "with_errors/dihedral_coords_position_with_errors_"
        full_error_fn = self.get_fn(error_fn)
        with pytest.raises(TypeError):
            write_xyz_file_from_gaussian_coordinates(
                atom_namesList, full_error_fn, extension, "./", 1
            )

    def test_check_gaussian_angle_energy_file_correct(self):
        full_path = self.get_fn(
            "gaussian_style_output_files/CT_CT_C_OH/output/dihedral.txt"
        )
        assert check_gaussian_angle_energy_file_correct(full_path)

    def test_check_gaussian_optimized_coordinate_file_correct(self):
        full_path = self.get_fn(
            "gaussian_style_output_files/CT_CT_C_OH/output/dihedral_coords_position_36.txt"
        )
        assert check_gaussian_optimized_coordinate_file_correct(full_path)

    def test_get_final_gaussian_output_file_data(self):
        full_path = self.get_fn("qm_files/")
        in_indices = [3, 1, 2, 8]
        out = get_final_gaussian_output_file_data(
            {full_path: list(map(int, np.arange(32, dtype=int)))}, in_indices
        )
        (
            angles,
            energies,
            coords,
            elements,
            n_atoms,
            out_indices,
        ) = get_final_gaussian_output_file_data(
            {full_path: list(np.arange(32, dtype=int))}, [3, 1, 2, 8]
        )
        assert out
        assert in_indices == out_indices
        assert angles == ["150.0", "160.0", "170.0", "-180.0"]
        assert np.allclose(
            np.array(energies).astype(float),
            np.array(
                [-79.2264709151, -79.2276351404, -79.2284591083, -79.2287549865]
            ),
        )
        assert np.shape(coords) == (
            4,
            n_atoms,
            3,
        )  # 4 to grab from, 8 total atoms,
        assert elements == ["C", "C", "H", "H", "H", "H", "H", "H"]
        assert n_atoms == 8

        # Wrong input types
        with pytest.raises(TypeError):
            get_final_gaussian_output_file_data(None, [3, 1, 2, 8])

        # Dict value must be list of integers
        with pytest.raises(TypeError):
            get_final_gaussian_output_file_data(
                {full_path: [0.1, 0.2, 0.3, -0.4]}, [3, 1, 2, 8]
            )

        # Dict key must be a string
        with pytest.raises(TypeError):
            get_final_gaussian_output_file_data(
                {1: [0.1, 0.2, 0.3, 0.4]}, [3, 1, 2, 8]
            )

        # Dict value must be a list
        with pytest.raises(TypeError):
            get_final_gaussian_output_file_data(
                {full_path: "ERROR"}, [3, 1, 2, 8]
            )

        # manual_dihedral_atom_numbers_list must be a list of len 4
        with pytest.raises(TypeError):
            get_final_gaussian_output_file_data(
                {full_path: list(np.arange(32, dtype=int))}, None
            )

        # manual_dihedral_atom_numbers_list must be a list of integers
        with pytest.raises(TypeError):
            get_final_gaussian_output_file_data(
                {full_path: list(np.arange(32, dtype=int))},
                [0.1, 0.2, 0.3, 0.4],
            )

    def test_get_gaussian_log_file_data_1(self):
        full_path = self.get_fn(
            "gaussian/CT_CT_C_OH/output/CT_CT_C_OH_multiplicity_1.log"
        )
        out_indices = []
        out = get_gaussian_log_file_data({full_path: out_indices})
        anglesList = out[0]
        anglesList = list(map(float, anglesList))
        expected_angles = list(np.arange(0, 180, 10))
        expected_angles_rev = list(-1 * np.arange(0, 181, 10))
        expected_angles_rev.reverse()
        assert np.allclose(anglesList, expected_angles + expected_angles_rev)

        energyList = out[1]
        expected_energies = [
            -266.838410009,
            -266.843802829,
            -266.843993151,
            -266.844262131,
            -266.844544987,
            -266.844773399,
            -266.844890793,
            -266.844872383,
            -266.844739667,
            -266.844561254,
            -266.844428251,
            -266.844414788,
            -266.844555389,
            -266.844842982,
            -266.845243308,
            -266.845697317,
            -266.846125223,
            -266.846434752,
            -266.846547900,
            -266.846434830,
            -266.846125228,
            -266.845697147,
            -266.845242984,
            -266.844843027,
            -266.844555543,
            -266.844414778,
            -266.844428006,
            -266.844561584,
            -266.844739882,
            -266.844872319,
            -266.844891002,
            -266.844773677,
            -266.844545113,
            -266.844262294,
            -266.843993903,
            -266.843802982,
            -266.843733993,
        ]
        energyList = list(map(float, energyList))
        assert np.allclose(energyList, expected_energies)

        coordsList = out[2]
        assert np.shape(coordsList) == (37, 11, 3)

        elementsList = out[3]
        assert elementsList == [
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

        n_atoms = out[4]
        assert n_atoms == len(elementsList)

        dihedral_atoms = out[5]
        assert dihedral_atoms == [5, 1, 2, 3]

    def test_get_gaussian_log_file_data_2(self):
        full_path = self.get_fn(
            "gaussian/perfluorodimethylether/scan_perfluorodimethylether.log"
        )
        out_indices = []
        out = get_gaussian_log_file_data({full_path: out_indices})
        anglesList = out[0]
        anglesList = list(map(float, anglesList))
        expected_angles = np.arange(round(-77.802, 3), 2.198, 5)
        expected_angles = np.append(
            expected_angles, np.arange(2.198, 182.198, 5)
        )
        expected_angles = np.append(
            expected_angles, np.arange(-177.802, -72.802, 5)
        )
        expected_angles = list(expected_angles.round(decimals=3))
        assert np.allclose(anglesList, expected_angles)

        energyList = out[1]
        expected_energies = [
            -747.265516312,
            -747.265443737,
            -747.265268941,
            -747.265123260,
            -747.265140067,
            -747.265300729,
            -747.265463069,
            -747.265515062,
            -747.265426543,
            -747.265221767,
            -747.264947686,
            -747.264653189,
            -747.264377973,
            -747.264149999,
            -747.263987395,
            -747.263901001,
            -747.263896025,
            -747.263972758,
            -747.264126609,
            -747.264347481,
            -747.264618251,
            -747.264912293,
            -747.265191319,
            -747.265407157,
            -747.265511705,
            -747.265476532,
            -747.265323201,
            -747.265154212,
            -747.265115120,
            -747.265246555,
            -747.265426192,
            -747.265515212,
            -747.265459576,
            -747.265275962,
            -747.265010909,
            -747.264715229,
            -747.264431277,
            -747.264189813,
            -747.264011609,
            -747.263909826,
            -747.263891551,
            -747.263958437,
            -747.264106693,
            -747.264326464,
            -747.264600305,
            -747.264900629,
            -747.265187236,
            -747.265408739,
            -747.265512974,
            -747.265469070,
            -747.265304125,
            -747.265138596,
            -747.265126682,
            -747.265280857,
            -747.265454103,
            -747.265515723,
            -747.265428483,
            -747.265218549,
            -747.264936922,
            -747.264635857,
            -747.264357118,
            -747.264129715,
            -747.263972162,
            -747.263895133,
            -747.263903100,
            -747.263995060,
            -747.264164599,
            -747.264399353,
            -747.264679582,
            -747.264975866,
            -747.265247223,
            -747.265443415,
            -747.265516312,
        ]
        energyList = list(map(float, energyList))
        assert np.allclose(energyList, expected_energies)

        coordsList = out[2]
        assert np.shape(coordsList) == (73, 9, 3)

        elementsList = out[3]
        assert elementsList == ["C", "F", "F", "F", "C", "F", "F", "F", "O"]

        n_atoms = out[4]
        assert n_atoms == len(elementsList)

        dihedral_atoms = out[5]
        assert dihedral_atoms == [3, 1, 9, 5]

    def test_write_qm_data_files(self):
        full_path = self.get_fn(
            "gaussian/CT_CT_C_OH/output/CT_CT_C_OH_multiplicity_1.log"
        )
        out_indices = []
        write_qm_data_files({full_path: out_indices})
        assert os.path.exists("extracted_gaussian_data/dihedral.txt")

        # qm_engine not in list
        with pytest.raises(ValueError):
            write_qm_data_files(None, qm_engine="None")

        # when qm_engine="gaussian_style_final_files"
        # manual_dihedral_atom_numbers_list must be a list of len 4
        with pytest.raises(TypeError):
            write_qm_data_files(
                {full_path: out_indices},
                None,
                qm_engine="gaussian_style_final_files",
            )

    def test_get_matching_dihedral_info_and_opls_fitting_data(self):
        out = get_matching_dihedral_info_and_opls_fitting_data(
            fit_dihedral_atom_types=["HC", "CT", "CT", "HC"],
            psf_path_and_filename=self.get_fn(
                "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files.psf"
            ),
            qm_log_file_dict={
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                ): []
            },
            mol2_file=self.get_fn(
                "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
            ),
            qm_engine="gaussian",
            manual_dihedral_atom_numbers_list=None,
        )
        dihedral_numsList = out[0]
        expected_dihedralnumsList = [
            [3, 1, 2, 6],
            [3, 1, 2, 7],
            [3, 1, 2, 8],
            [4, 1, 2, 6],
            [4, 1, 2, 7],
            [4, 1, 2, 8],
            [5, 1, 2, 6],
            [5, 1, 2, 7],
            [5, 1, 2, 8],
        ]
        assert dihedral_numsList == expected_dihedralnumsList

        dihedral_typesList = out[1]
        expected_dihedraltypesList = [
            ["HC", "CT", "CT", "HC"],
            ["HC", "CT", "CT", "HC"],
            ["HC", "CT", "CT", "HC"],
            ["HC", "CT", "CT", "HC"],
            ["HC", "CT", "CT", "HC"],
            ["HC", "CT", "CT", "HC"],
            ["HC", "CT", "CT", "HC"],
            ["HC", "CT", "CT", "HC"],
            ["HC", "CT", "CT", "HC"],
        ]
        assert dihedral_typesList == expected_dihedraltypesList

        coordsList = out[2]
        expected_coords = [
            [-1.012159, 1.156346, 0.0],
            [0.0, 0.76376, 0.0],
            [-0.0, -0.76376, 0.0],
            [1.012159, -1.156346, 0.0],
        ]
        assert np.allclose(coordsList[0][0], expected_coords)
        assert np.shape(coordsList) == (37, 9, 4, 3)
        degreesList = out[3]
        expected_degrees = [
            -180.0,
            -60.00001300151252,
            60.00001300151252,
            -60.00001300151252,
            59.99997399697496,
            -180.0,
            60.00001300151252,
            -180.0,
            -59.99997399697496,
        ]
        assert np.allclose(degreesList[0], expected_degrees)
        assert np.shape(degreesList) == (37, 9)
        opls_paramsList = out[4]
        expected_opls_params = [
            1,
            8.999999999999845,
            8.99999999999938,
            2.780442542871242e-12,
            8.99999999999753,
        ]
        assert np.allclose(opls_paramsList[0], expected_opls_params)

        # fit_dihedral_atom_types must be a list of length 4
        with pytest.raises(TypeError):
            out = get_matching_dihedral_info_and_opls_fitting_data(
                fit_dihedral_atom_types=["HC", "CT", "CT"],
                psf_path_and_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files.psf"
                ),
                qm_log_file_dict={
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): []
                },
                mol2_file=self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                qm_engine="gaussian",
                manual_dihedral_atom_numbers_list=None,
            )

        # Each element of fit_dihedral_atom_types must be a string
        with pytest.raises(TypeError):
            out = get_matching_dihedral_info_and_opls_fitting_data(
                fit_dihedral_atom_types=[1, 2, 3, 4],
                psf_path_and_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files.psf"
                ),
                qm_log_file_dict={
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): []
                },
                mol2_file=self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                qm_engine="gaussian",
                manual_dihedral_atom_numbers_list=None,
            )

        # fit_dihedral_atom_types must be a list
        with pytest.raises(TypeError):
            out = get_matching_dihedral_info_and_opls_fitting_data(
                fit_dihedral_atom_types=None,
                psf_path_and_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files.psf"
                ),
                qm_log_file_dict={
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): []
                },
                mol2_file=self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                qm_engine="gaussian",
                manual_dihedral_atom_numbers_list=None,
            )

        # When qm_engine="gaussian_style_final_files"
        # manual_dihedral_atom_numbers_list must be a list of len 4
        with pytest.raises(TypeError):
            out = get_matching_dihedral_info_and_opls_fitting_data(
                fit_dihedral_atom_types=["HC", "CT", "CT", "HC"],
                psf_path_and_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files.psf"
                ),
                qm_log_file_dict={
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): []
                },
                mol2_file=self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                qm_engine="gaussian_style_final_files",
                manual_dihedral_atom_numbers_list=None,
            )

        # qm_engine must be "gaussian" or "gaussian_style_final_files"
        with pytest.raises(ValueError):
            out = get_matching_dihedral_info_and_opls_fitting_data(
                fit_dihedral_atom_types=["HC", "CT", "CT", "HC"],
                psf_path_and_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files.psf"
                ),
                qm_log_file_dict={
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): []
                },
                mol2_file=self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                qm_engine="None",
                manual_dihedral_atom_numbers_list=None,
            )

    def test_change_gomc_ff_file_dihedral_values(self):
        new_file = self.get_fn(
            "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files_dihedrals_zeroed.inp"
        )
        change_gomc_ff_file_dihedral_values(
            read_gomc_ff_filename=self.get_fn(
                "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp"
            ),
            new_gomc_ff_filename=new_file,
            fit_dihedral_atom_types=["HC", "CT", "CT", "HC"],
            fit_dihedral_opls_k_0_1_2_3_4_values=[1, 0, 0, 0, 0],
        )
        assert os.path.exists(new_file)
        checkFile = False
        with open(new_file, "r") as f:
            for line in f:
                if (
                    line.replace(" ", "")
                    == "HCCTCTHC0.5090.0!TMP_opls_140TMP_opls_135TMP_opls_135TMP_opls_140\n"
                ):
                    checkFile = True
                    break

        assert checkFile

        with pytest.raises(TypeError):
            change_gomc_ff_file_dihedral_values(
                read_gomc_ff_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp"
                ),
                new_gomc_ff_filename=new_file,
                fit_dihedral_atom_types=[],
                fit_dihedral_opls_k_0_1_2_3_4_values=[1, 0, 0, 0],
            )

        # zero_dihedral_atom_types must be a list or None
        with pytest.raises(TypeError):
            change_gomc_ff_file_dihedral_values(
                read_gomc_ff_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp"
                ),
                new_gomc_ff_filename=new_file,
                fit_dihedral_atom_types=["HC", "CT", "CT", "HC"],
                fit_dihedral_opls_k_0_1_2_3_4_values=[1, 0, 0, 0, 0],
                zero_dihedral_atom_types="INVALID",
            )

        with pytest.raises(TypeError):
            change_gomc_ff_file_dihedral_values(
                read_gomc_ff_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp"
                ),
                new_gomc_ff_filename=new_file,
                fit_dihedral_atom_types=["HC", "CT", "CT", "HC"],
                fit_dihedral_opls_k_0_1_2_3_4_values=[1, 0, 0, 0, 0],
                zero_dihedral_atom_types=["CT"],
            )

        with pytest.raises(TypeError):
            change_gomc_ff_file_dihedral_values(
                read_gomc_ff_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp"
                ),
                new_gomc_ff_filename=new_file,
                fit_dihedral_atom_types=["HC", "CT", "CT", "HC"],
                fit_dihedral_opls_k_0_1_2_3_4_values=[1, 0, 0, 0, 0],
                zero_dihedral_atom_types=[1, 2, 3, 4],
            )

        with pytest.raises(ValueError):
            change_gomc_ff_file_dihedral_values(
                read_gomc_ff_filename=self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/GOMC_pdb_psf_ff_files_dihedrals_per_xml.inp"
                ),
                new_gomc_ff_filename=new_file,
                fit_dihedral_atom_types=["HC", "CT", "CT", "HC"],
                fit_dihedral_opls_k_0_1_2_3_4_values=[1, 0, 0, 0, 0],
                zero_dihedral_atom_types=["HC", "CT", "CT", "HC"],
            )
