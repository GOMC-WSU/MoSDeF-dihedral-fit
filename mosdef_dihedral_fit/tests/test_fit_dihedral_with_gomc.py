import os

import numpy as np
import pytest
import unyt as u

import mosdef_dihedral_fit.utils.math_operations as mdf_math
from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import (
    fit_dihedral_with_gomc,
)
from mosdef_dihedral_fit.tests.base_test import BaseTest

# The below "gomc_binary_directory = "/opt/gomc_build/GOMC/bin"" is the
# automated GitHub test binary path. It has to be fully specified, even if it is in the bashrc file.
# The automated GitHub test filename and location = "MoSDeF-dihedral-fit/.github/workflows/CI.yml"
gomc_binary_directory = "/opt/gomc_build/GOMC/bin"

# User changable variable, as it needs to be run locally
# Examples
gomc_binary_directory = "/home/brad/Programs/GOMC/GOMC-2.75a/bin"
# gomc_binary_directory = "/Users/calcraven/Documents/Vanderbilt/Research/MoSDeF/Dihedral_Fitter/GOMC/bin"


# NOTE: When comparing fitted values with reference value,
# we are using numpy.isclose() with absolute tolerance of
# 0.02 and relative tolerance of 0.08 (8%) to account for
# difference that incur across operating system.
class TestFitDihedralWithGomc(BaseTest):
    def test_check_if_gomc_binary_exists(self):
        assert os.path.isfile(f"{gomc_binary_directory}/GOMC_CPU_NVT") is True

    def test_gaussian_log_file_fit_oplsaa_fit_ethane_HC_CT_CT_HC(self):
        fit_dihedral_with_gomc(
            ["HC", "CT", "CT", "HC"],
            self.get_fn(
                "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
            ),
            self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                ): [0],
            },
            zero_dihedral_atom_types=None,
            qm_engine="gaussian",
            combining_rule="lorentz",
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=1e-03,
            opls_force_k0_zero=True
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is False
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str(3), 0, 0, 0, 0.31400374842, 0, 0.998791145574],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_3"),
                        0.31400374842,
                        0,
                        0,
                        -0.15700187421,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998791145574,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )

        # check the periodic dihedral file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_3"),
                        0.15700187421,
                        0.47100562263,
                        0,
                        -0.62800749684,
                        0,
                        0,
                        0.998791145574,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )
    
    def test_gaussian_log_file_fit_oplsaa_fit_ethane_HC_CT_CT_HC_with_2_log_files(
        self,
    ):
        fit_dihedral_with_gomc(
            ["HC", "CT", "CT", "HC"],
            self.get_fn(
                "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
            ),
            self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                ): [0],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log"
                ): [0],
            },
            zero_dihedral_atom_types=None,
            qm_engine="gaussian",
            combining_rule="geometric",
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=1e-03,
            opls_force_k0_zero=True
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is False
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is False
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt")
            is False
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str(3), 0, 0, 0, 0.31400374842, 0, 0.998791145574],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_3"),
                        0.31400374842,
                        0,
                        0,
                        -0.15700187421,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998791145574,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )

        # check the periodic dihedral file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_3"),
                        0.15700187421,
                        0.47100562263,
                        0,
                        -0.62800749684,
                        0,
                        0,
                        0.998791145574,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )
    
    def test_gaussian_log_file_fit_oplsaa_protonated_fragment_CT_CT_C_OH_in_COOH(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CT", "CT", "C", "OH"],
            self.get_fn(
                "gaussian/CT_CT_C_OH/input/starting_coords/protonated_fragment_CT_CT_C_OH_in_COOH.mol2"
            ),
            self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian/CT_CT_C_OH/output/CT_CT_C_OH_multiplicity_1.log"
                ): [0],
            },
            zero_dihedral_atom_types=[["CT", "CT", "C", "O_3"]],
            qm_engine="gaussian",
            combining_rule="geometric",
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=1e-03,
            opls_force_k0_zero=True
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 2.00408568934, 0, 0, 0, -0.371866278201],
                    [str("2"), 0, 0, 2.13847194875, 0, 0, 0.264787158812],
                    [str("3"), 0, 0, 0, 1.79252468736, 0, -1.29046011346],
                    [str("4"), 0, 0, 0, 0, 1.71219173182, -1.61245472957],
                    [
                        str("1_3"),
                        0,
                        1.45632365052,
                        0,
                        0.821642101161,
                        0,
                        0.0570637738537,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.79461740051,
                        0,
                        0.515780222673,
                        0.433811455157,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        1.17191599323,
                        0.930913047097,
                        -0.739857285933,
                    ],
                    [
                        str("1_2"),
                        0,
                        1.04119433442,
                        1.444342119,
                        0,
                        0,
                        0.953575526751,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        0.925047594422,
                        1.32819621359,
                        0.290365264694, 
                        0,
                        0.998573206243,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        0.914079354784,
                        1.31722771457,
                        0.279396891914,
                        0.0383892606964,
                        0.999295535579,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        2.00408568934,
                        -1.00204284467,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.371866278201 ,
                    ],
                    [
                        str("0_2"),
                        2.13847194875,
                        0,
                        -1.06923597438,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.264787158812,
                    ],
                    [
                        str("0_3"),
                        1.79252468736,
                        0,
                        0,
                        -0.89626234368,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.29046011346,
                    ],
                    [
                        str("0_4"),
                        1.71219173182,
                        0,
                        0,
                        0,
                        -0.85609586591,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.61245472957,
                    ],
                    [
                        str("0_1_3"),
                        2.27796575168,
                        -0.72816182526,
                        0,
                        -0.41082105058,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.0570637738537,
                    ],
                    [
                        str("0_2_4"),
                        2.31039762318,
                        0,
                        -0.897308700255,
                        0,
                        -0.257890111337,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.433811455157,
                    ],
                    [
                        str("0_3_4"),
                        2.10282904033,
                        0,
                        0,
                        -0.585957996615,
                        -0.46545652354,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.739857285933,
                    ],
                    [
                        str("0_1_2"),
                        2.48553645342,
                        -0.52059716721,
                        -0.7221710595,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.953575526751,
                    ],
                    [
                        str("0_1_2_3"),
                        2.54360907271,
                        -0.462523797211,
                        -0.664098106795,
                        -0.145182632347,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998573206243,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.54909322196,
                        -0.457039677392,
                        -0.658613857285,
                        -0.139698445957,
                        -0.0191946303482,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.999295535579 ,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        1.0020428444,
                        -1.0020428444,
                        0,
                        0,
                        0,
                        0,
                        -0.371866278201,
                    ],
                    [
                        str("0_2"),
                        2.13847194795,
                        0,
                        -2.13847194795,
                        0,
                        0,
                        0,
                        0.264787158812,
                    ],
                    [
                        str("0_3"),
                        0.896262343155,
                        2.68878702946,
                        0,
                        -3.58504937262,
                        0,
                        0,
                        -1.29046011346,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.84876692728,
                        0,
                        -6.84876692728,
                        0,
                        -1.61245472957,
                    ],
                    [
                        str("0_1_3"),
                        1.13898287652,
                        0.504301323446,
                        0,
                        -1.64328419996,
                        0,
                        0,
                        0.0570637738537,
                    ],
                    [
                        str("0_2_4"),
                        1.79461739534,
                        0,
                        0.26850352834,
                        0,
                        -2.06312092368,
                        0,
                        0.433811455157,
                    ],
                    [
                        str("0_3_4"),
                        0.585957999835,
                        1.7578739995,
                        3.72365217288,
                        -2.34383199934,
                        -3.72365217288,
                        0,
                        -0.739857285933,
                    ],
                    [
                        str("0_1_2"),
                        1.96493928638,
                        -0.520597166825,
                        -1.44434211955,
                        0,
                        0,
                        0,
                        0.953575526751,
                    ],
                    [
                        str("0_1_2_3"),
                        1.93590264307,
                        -0.026975900098,
                        -1.32819621337,
                        -0.580730529606,
                        0,
                        0,
                        0.998573206243,
                    ],
                    [
                        str("0_1_2_3_4"),
                        1.91396583919,
                        -0.037944339811,
                        -1.1636706807,
                        -0.558793783876,
                        -0.153557034803,
                        0,
                        0.999295535579,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(
                                float(split_line_i[j]), number_sig_i
                            ) == mdf_math.round_to_sig_figs(
                                correct_line_values[i][j], number_sig_i
                            )
    
    def test_gaussian_style_files_fit_oplsaa_fit_CT_CT_C_OH_in_COOH_missing_1st_point(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CT", "CT", "C", "OH"],
            self.get_fn(
                "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
            ),
            self.get_fn("oplsaa_CT_CT_C_OH_in_COOH_zeroed.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn("gaussian_style_output_files/CT_CT_C_OH/output"): [
                    0
                ],
            },
            manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
            zero_dihedral_atom_types=None,
            qm_engine="gaussian_style_final_files",
            combining_rule="geometric",
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=0.02,
            opls_force_k0_zero=True
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 1.80634795241, 0, 0, 0, -0.268742410571],
                    [str("2"), 0, 0, 1.94879786206, 0, 0, 0.448421126946],
                    [str("3"), 0, 0, 0, 1.58083737128, 0, -1.50863674602],
                    [str("4"), 0, 0, 0, 0, 1.5121963942, -1.54542980028],
                    [
                        str("1_3"),
                        0,
                        1.41759854611,
                        0,
                        0.618102768602,
                        0,
                        -0.00379371805915,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.69725316221,
                        0,
                        0.37780223364,
                        0.551946157778,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        0.966769776613,
                        0.903643214082,
                        -0.897967473282,
                    ],
                    [
                        str("1_2"),
                        0,
                        0.882099545744,
                        1.34468116365,
                        0,
                        0,
                        0.986743110285,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        0.83099823616,
                        1.29385707544,
                        0.136793004885,
                        0,
                        0.997926382517,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        0.836836322925,
                        1.29869945716,
                        0.141702758859,
                        -0.0179044297589,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        1.80634795241,
                        -0.903173976205,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.268742410571,
                    ],
                    [
                        str("0_2"),
                        1.94879786206,
                        0,
                        -0.97439893103,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.448421126946,
                    ],
                    [
                        str("0_3"),
                        1.58083737128,
                        0,
                        0,
                        -0.79041868564,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.50863674602,
                    ],
                    [
                        str("0_4"),
                        1.5121963942,
                        0,
                        0,
                        0,
                        -0.7560981971,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.54542980028,
                    ],
                    [
                        str("0_1_3"),
                        2.03570131471,
                        -0.708799273055,
                        0,
                        -0.309051384301,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.00379371805915,
                    ],
                    [
                        str("0_2_4"),
                        2.07505539585,
                        0,
                        -0.848626581105,
                        0,
                        -0.18890111682,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.551946157778,
                    ],
                    [
                        str("0_3_4"),
                        1.8704129907,
                        0,
                        0,
                        -0.483384888306,
                        -0.451821607041,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.897967473282,
                    ],
                    [
                        str("0_1_2"),
                        2.22678070939,
                        -0.441049772872,
                        -0.672340581825,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.986743110285,
                    ],
                    [
                        str("0_1_2_3"),
                        2.26164831649,
                        -0.415499118084,
                        -0.64692853772,
                        -0.0683965024425,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.997926382517,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.25933410919,
                        -0.418418161462,
                        -0.64934972858,
                        -0.0708513794295,
                        0.00895221487945,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        0.903173976205,
                        -0.903173976205,
                        0,
                        0,
                        0,
                        0,
                        -0.268742410571,
                    ],
                    [
                        str("0_2"),
                        1.94879786206,
                        0,
                        -1.94879786206,
                        0,
                        0,
                        0,
                        0.448421126946,
                    ],
                    [
                        str("0_3"),
                        0.79041868564,
                        2.37125605692,
                        0,
                        -3.16167474256,
                        0,
                        0,
                        -1.50863674602,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.0487855768,
                        0,
                        -6.0487855768,
                        0,
                        -1.54542980028,
                    ],
                    [
                        str("0_1_3"),
                        1.01785065736,
                        0.218354879848,
                        0,
                        -1.2362055372,
                        0,
                        0,
                        -0.00379371805915,
                    ],
                    [
                        str("0_2_4"),
                        1.69725316221,
                        0,
                        -0.18604422765,
                        0,
                        -1.51120893456,
                        0,
                        0.551946157778,
                    ],
                    [
                        str("0_3_4"),
                        0.483384888306,
                        1.45015466492,
                        3.61457285633,
                        -1.93353955323,
                        -3.61457285633,
                        0,
                        -0.897967473282,
                    ],
                    [
                        str("0_1_2"),
                        1.78573093652,
                        -0.441049772872,
                        -1.34468116365,
                        0,
                        0,
                        0,
                        0.986743110285,
                    ],
                    [
                        str("0_1_2_3"),
                        1.77775269597,
                        -0.210309610756,
                        -1.29385707544,
                        -0.27358600977,
                        0,
                        0,
                        0.997926382517,
                    ],
                    [
                        str("0_1_2_3_4"),
                        1.78796899805,
                        -0.205864023174,
                        -1.3703171762,
                        -0.283405517718,
                        0.0716177190356,
                        0,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )
    
    def test_gaussian_style_files_fit_oplsaa_fit_CT_CT_C_OH_in_COOH(self):
        fit_dihedral_with_gomc(
            ["CT", "CT", "C", "OH"],
            self.get_fn(
                "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
            ),
            self.get_fn("oplsaa_CT_CT_C_OH_in_COOH_zeroed.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/output"
                ): [],
            },
            manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
            zero_dihedral_atom_types=None,
            qm_engine="gaussian_style_final_files",
            combining_rule="geometric",
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=5e-03,
            opls_force_k0_zero=True
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 1.87103578693, 0, 0, 0, -0.184214001408],
                    [str("2"), 0, 0, 2.01338731788, 0, 0, 0.468487865197],
                    [str("3"), 0, 0, 0, 1.64674772188, 0, -1.32197179229],
                    [str("4"), 0, 0, 0, 0, 1.5768734149, -1.35119291697],
                    [
                        str("1_3"),
                        0,
                        1.45815624076,
                        0,
                        0.656469146108,
                        0,
                        0.0786991845116,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.73604078977,
                        0,
                        0.416554745122,
                        0.579202191408,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        1.00514371097,
                        0.944164973554,
                        -0.735494529166,
                    ],
                    [
                        str("1_2"),
                        0,
                        0.920441371308,
                        1.38301174904,
                        0,
                        0,
                        0.984123002941,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        0.858889958603,
                        1.32179425273,
                        0.164766862242,
                        0,
                        0.998396296125,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        0.857564707345,
                        1.32069502687,
                        0.163652345714,
                        0.00406432148459,
                        0.998404336839,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        1.87103578693,
                        -0.935517893465,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.184214001408,
                    ],
                    [
                        str("0_2"),
                        2.01338731788,
                        0,
                        -1.00669365894,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.468487865197,
                    ],
                    [
                        str("0_3"),
                        1.64674772188,
                        0,
                        0,
                        -0.82337386094,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.32197179229,
                    ],
                    [
                        str("0_4"),
                        1.5768734149,
                        0,
                        0,
                        0,
                        -0.78843670745,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.35119291697,
                    ],
                    [
                        str("0_1_3"),
                        2.11462538687,
                        -0.72907812038,
                        0,
                        -0.328234573054,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.0786991845116,
                    ],
                    [
                        str("0_2_4"),
                        2.15259553489,
                        0,
                        -0.868020394885,
                        0,
                        -0.208277372561,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.579202191408,
                    ],
                    [
                        str("0_3_4"),
                        1.94930868452,
                        0,
                        0,
                        -0.502571855485,
                        -0.472082486777,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.735494529166,
                    ],
                    [
                        str("0_1_2"),
                        2.30345312035,
                        -0.460220685654,
                        -0.69150587452,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.984123002941,
                    ],
                    [
                        str("0_1_2_3"),
                        2.34545107358,
                        -0.429444979302,
                        -0.660897126365,
                        -0.082383431121,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998396296125,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.34597640141,
                        -0.428782353672,
                        -0.660347513435,
                        -0.081826172857,
                        -0.0020321607423,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998404336839,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        0.935517893465,
                        -0.935517893465,
                        0,
                        0,
                        0,
                        0,
                        -0.184214001408,
                    ],
                    [
                        str("0_2"),
                        2.01338731788,
                        0,
                        -2.01338731788,
                        0,
                        0,
                        0,
                        0.468487865197,
                    ],
                    [
                        str("0_3"),
                        0.82337386094,
                        2.47012158282,
                        0,
                        -3.29349544376,
                        0,
                        0,
                        -1.32197179229,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.3074936596,
                        0,
                        -6.3074936596,
                        0,
                        -1.35119291697,
                    ],
                    [
                        str("0_1_3"),
                        1.05731269343,
                        0.255625598782,
                        0,
                        -1.31293829222,
                        0,
                        0,
                        0.0786991845116,
                    ],
                    [
                        str("0_2_4"),
                        1.73604078977,
                        0,
                        -0.069821809282,
                        0,
                        -1.66621898049,
                        0,
                        0.579202191408,
                    ],
                    [
                        str("0_3_4"),
                        0.502571855485,
                        1.50771556646,
                        3.77665989422,
                        -2.01028742194,
                        -3.77665989422,
                        0,
                        -0.735494529166,
                    ],
                    [
                        str("0_1_2"),
                        1.84323243469,
                        -0.460220685654,
                        -1.38301174904,
                        0,
                        0,
                        0,
                        0.984123002941,
                    ],
                    [
                        str("0_1_2_3"),
                        1.83362266315,
                        -0.182294685938,
                        -1.32179425273,
                        -0.329533724484,
                        0,
                        0,
                        0.998396296125,
                    ],
                    [
                        str("0_1_2_3_4"),
                        1.8313035534,
                        -0.183303835102,
                        -1.30443774093,
                        -0.327304691428,
                        -0.0162572859384,
                        0,
                        0.998404336839,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

    def test_gaussian_style_files_fit_oplsaa_fit_CT_CT_C_OH_in_COOH_2_files_missing_1_first_point(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CT", "CT", "C", "OH"],
            self.get_fn(
                "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
            ),
            self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                ): [],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                ): [0],
            },
            manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
            zero_dihedral_atom_types=[["CT", "CT", "C", "O_3"]],
            qm_engine="gaussian_style_final_files",
            combining_rule="geometric",
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=0.02,
            opls_force_k0_zero=True
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 1.80634795241, 0, 0, 0, -0.268742410571],
                    [str("2"), 0, 0, 1.94879786206, 0, 0, 0.448421126946],
                    [str("3"), 0, 0, 0, 1.58083737128, 0, -1.50863674602],
                    [str("4"), 0, 0, 0, 0, 1.5121963942, -1.54542980028],
                    [
                        str("1_3"),
                        0,
                        1.41759854611,
                        0,
                        0.618102768602,
                        0,
                        -0.00379371805915,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.69725316221,
                        0,
                        0.37780223364,
                        0.551946157778,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        0.966769776613,
                        0.903643214082,
                        -0.897967473282,
                    ],
                    [
                        str("1_2"),
                        0,
                        0.882099545744,
                        1.34468116365,
                        0,
                        0,
                        0.986743110285,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        0.83099823616,
                        1.29385707544,
                        0.136793004885,
                        0,
                        0.997926382517,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        0.836836322925,
                        1.29869945716,
                        0.141702758859,
                        -0.0179044297589,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        1.80634795241,
                        -0.903173976205,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.268742410571,
                    ],
                    [
                        str("0_2"),
                        1.94879786206,
                        0,
                        -0.97439893103,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.448421126946,
                    ],
                    [
                        str("0_3"),
                        1.58083737128,
                        0,
                        0,
                        -0.79041868564,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.50863674602,
                    ],
                    [
                        str("0_4"),
                        1.5121963942,
                        0,
                        0,
                        0,
                        -0.7560981971,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.54542980028,
                    ],
                    [
                        str("0_1_3"),
                        2.03570131471,
                        -0.708799273055,
                        0,
                        -0.309051384301,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.00379371805915,
                    ],
                    [
                        str("0_2_4"),
                        2.07505539585,
                        0,
                        -0.848626581105,
                        0,
                        -0.18890111682,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.551946157778,
                    ],
                    [
                        str("0_3_4"),
                        1.8704129907,
                        0,
                        0,
                        -0.483384888306,
                        -0.451821607041,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.897967473282,
                    ],
                    [
                        str("0_1_2"),
                        2.22678070939,
                        -0.441049772872,
                        -0.672340581825,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.986743110285,
                    ],
                    [
                        str("0_1_2_3"),
                        2.26164831649,
                        -0.415499118084,
                        -0.64692853772,
                        -0.0683965024425,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.997926382517,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.25933410919,
                        -0.418418161462,
                        -0.64934972858,
                        -0.0708513794295,
                        0.00895221487945,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        0.903173976205,
                        -0.903173976205,
                        0,
                        0,
                        0,
                        0,
                        -0.268742410571,
                    ],
                    [
                        str("0_2"),
                        1.94879786206,
                        0,
                        -1.94879786206,
                        0,
                        0,
                        0,
                        0.448421126946,
                    ],
                    [
                        str("0_3"),
                        0.79041868564,
                        2.37125605692,
                        0,
                        -3.16167474256,
                        0,
                        0,
                        -1.50863674602,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.0487855768,
                        0,
                        -6.0487855768,
                        0,
                        -1.54542980028,
                    ],
                    [
                        str("0_1_3"),
                        1.01785065736,
                        0.218354879848,
                        0,
                        -1.2362055372,
                        0,
                        0,
                        -0.00379371805915,
                    ],
                    [
                        str("0_2_4"),
                        1.69725316221,
                        0,
                        -0.18604422765,
                        0,
                        -1.51120893456,
                        0,
                        0.551946157778,
                    ],
                    [
                        str("0_3_4"),
                        0.483384888306,
                        1.45015466492,
                        3.61457285633,
                        -1.93353955323,
                        -3.61457285633,
                        0,
                        -0.897967473282,
                    ],
                    [
                        str("0_1_2"),
                        1.78573093652,
                        -0.441049772872,
                        -1.34468116365,
                        0,
                        0,
                        0,
                        0.986743110285,
                    ],
                    [
                        str("0_1_2_3"),
                        1.77775269597,
                        -0.210309610756,
                        -1.29385707544,
                        -0.27358600977,
                        0,
                        0,
                        0.997926382517,
                    ],
                    [
                        str("0_1_2_3_4"),
                        1.78796899805,
                        -0.205864023174,
                        -1.3703171762,
                        -0.283405517718,
                        0.0716177190356,
                        0,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

    def test_gaussian_log_file_fit_oplsaa_protonated_fragment_CT_CT_C_OH_in_COOH_bad_element_order_mol2(
        self,
    ):
        mol2_elements_values = (
            "\['C', 'C', 'O', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H'\]"
        )
        qm_elements_values = (
            "\['C', 'C', 'C', 'O', 'O', 'H', 'H', 'H', 'H', 'H', 'H'\]"
        )
        with pytest.raises(
            ValueError,
            match=f"ERROR: The QM elements do not match the mol2 file elements, in order. \n"
            f"This does not guarantee that the element postions are correct. \n"
            f"mol2 file element names = {mol2_elements_values} \n"
            f"QM file element names = {qm_elements_values}",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian/CT_CT_C_OH/input/starting_coords/"
                    "protonated_fragment_CT_CT_C_OH_in_COOH_bad_element_order.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/CT_CT_C_OH/output/CT_CT_C_OH_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=[["CT", "CT", "C", "O_3"]],
                qm_engine="gaussian",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_gaussian_style_files_fit_oplsaa_protonated_fragment_CT_CT_C_OH_in_COOH_bad_element_order_mol2(
        self,
    ):
        mol2_elements_values = (
            "\['C', 'C', 'O', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H'\]"
        )
        qm_elements_values = (
            "\['C', 'C', 'C', 'O', 'O', 'H', 'H', 'H', 'H', 'H', 'H'\]"
        )
        with pytest.raises(
            ValueError,
            match=f"ERROR: The QM elements do not match the mol2 file elements, in order. \n"
            f"This does not guarantee that the element postions are correct. \n"
            f"mol2 file element names = {mol2_elements_values} \n"
            f"QM file element names = {qm_elements_values}",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/"
                    "starting_coords/CT_CT_C_3_OH_bad_element_order.mol2"
                ),
                self.get_fn("oplsaa_CT_CT_C_OH_in_COOH_zeroed.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH/output"
                    ): [],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=None,
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=5e-03,
                opls_force_k0_zero=True
            )
    
    def test_gaussian_log_file_variable_VDWGeometricSigma_default(self):
        fit_dihedral_with_gomc(
            ["HC", "CT", "CT", "HC"],
            self.get_fn(
                "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
            ),
            self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                ): [0],
            },
            zero_dihedral_atom_types=None,
            qm_engine="gaussian",
            combining_rule=None,
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=1e-03,
            opls_force_k0_zero=True
        )

        with open(
            "GOMC_simulations/GOMC_OPLS_fit_3_dihedral_coords_1.conf", "r"
        ) as fp:
            variables_read_dict = {
                "VDWGeometricSigma": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                    pass

            assert variables_read_dict == {
                "VDWGeometricSigma": True,
            }

    def test_gaussian_log_file_variable_VDWGeometricSigma_True(self):
        fit_dihedral_with_gomc(
            ["HC", "CT", "CT", "HC"],
            self.get_fn(
                "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
            ),
            self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                ): [0],
            },
            zero_dihedral_atom_types=None,
            qm_engine="gaussian",
            combining_rule=None,
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=1e-03,
            opls_force_k0_zero=True
        )

        with open(
            "GOMC_simulations/GOMC_OPLS_fit_3_dihedral_coords_1.conf", "r"
        ) as fp:
            variables_read_dict = {
                "VDWGeometricSigma": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                    pass

            assert variables_read_dict == {
                "VDWGeometricSigma": True,
            }

    def test_gaussian_log_file_variable_VDWGeometricSigma_False(self):
        fit_dihedral_with_gomc(
            ["HC", "CT", "CT", "HC"],
            self.get_fn(
                "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
            ),
            self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                ): [0],
            },
            zero_dihedral_atom_types=None,
            qm_engine="gaussian",
            combining_rule="lorentz",
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=1e-03,
            opls_force_k0_zero=True
        )

        with open(
            "GOMC_simulations/GOMC_OPLS_fit_3_dihedral_coords_1.conf", "r"
        ) as fp:
            variables_read_dict = {
                "VDWGeometricSigma": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                    pass

            assert variables_read_dict == {
                "VDWGeometricSigma": True,
            }

    def test_bad_fit_dihedral_atom_types_input_list_of_3(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The input 'fit_dihedral_atom_types' variable = \['HC', 'CT', 'CT'\], "
            r"but it needs to be a list of 4 strings, "
            r"where the strings are the atom types/classes. Example: \['HC', 'CT', 'CT', 'HC'\].",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_bad_fit_opls_force_k0_zero_not_bool(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the 'opls_force_k0_zero' file as a bool.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero='x'
            )

    def test_bad_fit_dihedral_atom_types_input_list_of_4_with_int_at_0(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The input 'fit_dihedral_atom_types' variable = \[0, 'CT', 'CT', 'HC'\], "
            r"but it needs to be a list of 4 strings, "
            r"where the strings are the atom types/classes. Example: \['HC', 'CT', 'CT', 'HC'\].",
        ):
            fit_dihedral_with_gomc(
                [0, "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_bad_fit_dihedral_atom_types_input_list_of_4_with_int_at_1(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The input 'fit_dihedral_atom_types' variable = \['HC', 1, 'CT', 'HC'\], "
            r"but it needs to be a list of 4 strings, "
            r"where the strings are the atom types/classes. Example: \['HC', 'CT', 'CT', 'HC'\].",
        ):
            fit_dihedral_with_gomc(
                ["HC", 1, "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_bad_fit_dihedral_atom_types_input_list_of_4_with_int_at_2(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The input 'fit_dihedral_atom_types' variable = \['HC', 'CT', 2, 'HC'\], "
            r"but it needs to be a list of 4 strings, "
            r"where the strings are the atom types/classes. Example: \['HC', 'CT', 'CT', 'HC'\].",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", 2, "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_bad_fit_dihedral_atom_types_input_list_of_4_with_int_at_3(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The input 'fit_dihedral_atom_types' variable = \['HC', 'CT', 'CT', 3\], "
            r"but it needs to be a list of 4 strings, "
            r"where the strings are the atom types/classes. Example: \['HC', 'CT', 'CT', 'HC'\].",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", 3],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_mol2_file_file_does_not_exist(self):
        value_path_mol2 = "bad_mol2_path.mol2"
        with pytest.raises(
            ValueError,
            match=f"ERROR: The {value_path_mol2} file "
            r"\('mol2_file'\) does not exists.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                value_path_mol2,
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_mol2_file_file_no_mol2_extention(self):
        value_path_mol2 = "bad_mol2_path"
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please enter enter mol2 file \('mol2_file'\) name with the .mol2 extension.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                value_path_mol2,
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )
    
    def test_mol2_file_file_not_a_string(self):
        value_path_mol2 = 1
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter mol2 file \('mol2_file'\) as a string.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                value_path_mol2,
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_xml_selection_file_does_not_exist(self):
        value_path_xml = "bad_xml_path.xml"
        with pytest.raises(
            ValueError,
            match=f"ERROR: The {value_path_xml} file "
            r"\('forcefield_file'\) does not exists.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                value_path_xml,
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_xml_selection_file_no_xml_extention(self):
        value_path_xml = "bad_xml_path"
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please enter enter xml file "
            r"\('forcefield_file'\) name with the .xml extension.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                value_path_xml,
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_xml_selection_file_not_a_string(self):
        value_path_xml = 1
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter xml file \('forcefield_file'\) as a string.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                value_path_xml,
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_temperature_unyt_units_not_a_temperture_but_pressure(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'temperature_unyt_units' is not temperature of type {type(u.unyt_quantity)}.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.bar,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_temperature_unyt_units_not_in_unyt_units(self):
        with pytest.raises(
            TypeError,
            match=f"ERROR: The 'temperature_unyt_units' is not temperature of type {type(u.unyt_quantity)}.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_qm_log_files_and_entries_not_a_dict(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                ["x"],
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_qm_log_files_and_entries_key_1_not_a_string(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    1: [],
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log"
                    ): [0, 1],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_qm_log_files_and_entries_key_2_not_a_string(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [],
                    2: [0, 1],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_qm_log_files_and_entries_value_1_not_a_list(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): "s",
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log"
                    ): [0, 1],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )
    
    def test_qm_log_files_and_entries_value_2_not_a_list(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [],
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log"
                    ): "x",
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_qm_log_files_and_entries_list_1_not_all_int(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [0, "s"],
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log"
                    ): "x",
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_qm_log_files_and_entries_list_2_not_all_int(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [],
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log"
                    ): [0, 5, "s"],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_qm_log_files_and_entries_list_1_int_less_than_0(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [-1],
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log"
                    ): [0, 5],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_qm_log_files_and_entries_list_2_int_less_than_0(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'qm_log_file_dict' is not a dict "
            r"with a string keys and list of int>=0 as the values. Example: "
            r"\{'path/HC_CT_CT_HC_part_1.log'\): \[\], 'path/HC_CT_CT_HC_part_2.log'\): \[0, 5\]\}",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [],
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log"
                    ): [0, -5],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_gomc_binary_path_not_a_string(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the 'gomc_binary_path' file as a string.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                99999,
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_gomc_binary_path_containing_the_GOMC_CPU_NVT_file_does_not_exist(
        self,
    ):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The 'gomc_binary_path' file does not exist or contain the GOMC 'GOMC_CPU_NVT' file.",
        ):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                f"gomc_binary_directory",
                {
                    self.get_fn(
                        "gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log"
                    ): [],
                },
                zero_dihedral_atom_types=None,
                qm_engine="gaussian",
                combining_rule="lorentz",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1e-03,
                opls_force_k0_zero=True
            )

    def test_zero_dihedral_atom_types_not_list(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'zero_dihedral_atom_types' is not None or a list containing "
            r"lists with 4 strings each. Example: "
            r"\[\['CT', 'CT, 'CT, 'HC'\], \['NT', 'CT, 'CT, 'HC'\]\].",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types="str",
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_zero_dihedral_atom_types_list_1_str(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'zero_dihedral_atom_types' is not None or a list containing "
            r"lists with 4 strings each. Example: "
            r"\[\['CT', 'CT, 'CT, 'HC'\], \['NT', 'CT, 'CT, 'HC'\]\].",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=["str", ["HC", "CT", "CT", "C"]],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_zero_dihedral_atom_types_list_2_str(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'zero_dihedral_atom_types' is not None or a list containing "
            r"lists with 4 strings each. Example: "
            r"\[\['CT', 'CT, 'CT, 'HC'\], \['NT', 'CT, 'CT, 'HC'\]\].",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[["CT", "CT", "C", "O_3"], "str"],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_zero_dihedral_atom_types_list_1_not_4_strings(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'zero_dihedral_atom_types' is not None or a list containing "
            r"lists with 4 strings each. Example: "
            r"\[\['CT', 'CT, 'CT, 'HC'\], \['NT', 'CT, 'CT, 'HC'\]\].",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", 1, "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_zero_dihedral_atom_types_list_2_not_4_strings(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'zero_dihedral_atom_types' is not None or a list containing "
            r"lists with 4 strings each. Example: "
            r"\[\['CT', 'CT, 'CT, 'HC'\], \['NT', 'CT, 'CT, 'HC'\]\].",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", 2, "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_zero_dihedral_atom_types_list_1_not_lenght_4(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'zero_dihedral_atom_types' is not None or a list containing "
            r"lists with 4 strings each. Example: "
            r"\[\['CT', 'CT, 'CT, 'HC'\], \['NT', 'CT, 'CT, 'HC'\]\].",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_zero_dihedral_atom_types_list_2_not_lenght_4(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The 'zero_dihedral_atom_types' is not None or a list containing "
            r"lists with 4 strings each. Example: "
            r"\[\['CT', 'CT, 'CT, 'HC'\], \['NT', 'CT, 'CT, 'HC'\]\].",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_qm_engine_not_correct_value(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'qm_engine' = {'x'}, which is not "
            f"any of the available options. "
            f"The options are 'gaussian' or 'gaussian_style_final_files'.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="x",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_qm_engine_not_a_string(self):
        with pytest.raises(
            TypeError,
            match=f"ERROR: The 'qm_engine' is a {type(['x'])}, but it needs to be a str.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine=["x"],
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_combining_rule_not_correct_value_is_list(self):
        with pytest.raises(
            ValueError,
            match="ERROR: Please enter the 'combining_rule' file as a string of 'geometric' or 'lorentz' or None.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule=["geometric"],
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_combining_rule_not_correct_value_wrong_string(self):
        with pytest.raises(
            ValueError,
            match="ERROR: Please enter the 'combining_rule' file as a string of 'geometric' or 'lorentz' or None.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="x",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_atom_type_naming_style_not_correct_value(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'atom_type_naming_style' = {'x'}, which is not "
            f"any of the available options. "
            f"The options are 'general' or 'all_unique'.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="x",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_atom_type_naming_style_not_a_string(self):
        with pytest.raises(
            TypeError,
            match=f"ERROR: The 'atom_type_naming_style' is a {type(['x'])}, but it needs to be a str.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style=["x"],
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_gomc_cpu_cores_not_correct_value(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'gomc_cpu_cores' = {0}, and it must be an int > 0.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=0,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_gomc_cpu_cores_not_a_int(self):
        with pytest.raises(
            TypeError,
            match=f"ERROR: The 'gomc_cpu_cores' is a {type(1.000)}, but it needs to be a int.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1.000,
                r_squared_min=0.99,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )
    
    def test_r_squared_min_not_correct_value_is_0(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'r_squared_min'= {0.00}, "
            f"but it must be a 0<float<1.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.00,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_r_squared_min_not_correct_value_is_1(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'r_squared_min'= {1.00}, "
            f"but it must be a 0<float<1.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=1.00,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_r_squared_min_not_a_float(self):
        with pytest.raises(
            TypeError,
            match=f"ERROR: The 'r_squared_min' is a {type(2)}, "
            f"but it must be a 0<float<1.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=2,
                r_squared_rtol=0.02,
                opls_force_k0_zero=True
            )

    def test_r_squared_rtol_not_correct_value_is_0(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'r_squared_rtol' = {0.00}, "
            f"but it must be a 0<float<1.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.00,
                opls_force_k0_zero=True
            )

    def test_r_squared_rtol_not_correct_value_is_1(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'r_squared_rtol' = {1.00}, "
            f"but it must be a 0<float<1.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=1.00,
                opls_force_k0_zero=True
            )
    
    def test_r_squared_rtol_not_a_float(self):
        with pytest.raises(
            TypeError,
            match=f"ERROR: The 'r_squared_rtol' is a {type(2)}, "
            f"but it must be a 0<float<1.",
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[
                    ["CT", "CT", "C", "O_3"],
                    ["HC", "CT", "CT", "C"],
                ],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=2,
                opls_force_k0_zero=True
            )

    def test_error_r_squared_min_and_r_squared_rtol_need_adjusted(
        self,
    ):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The calculated R-squared energy values from the fit types "
            f"do not match the validated case for 'r_squared_min' >= "
            f"{'0.99'}, "
            f"within the relative tolerance or 'r_squared_rtol' = "
            f"{'2e-07'}. \n"
            f"Constants_used = The calculated R-squared energy values from the fit type constants_used."
            f"R-squared_fitted = Gaussian minus GOMC with the selected dihedral set to zero \n"
            f"R-squared_new_dihedral = Gaussian minus GOMC with the selected individual dihedral added in GOMC \n"
            f"Abs\(delta\) = absolute_value\(R-squared_fitted - R-squared_new_dihedral\) \n"
            f"- The fits for all are shown here for all the dihedral combinations \n"
            f"\[opls_constants_used, R-squared_fitted, R-squared_new_dihedral, Abs\(delta\)\] \n"
            f"\[\['1', -0.27129793, -1.2824201, 1.0111222\], "
            f"\['2', 0.45284933, -0.069821698, 0.52267103\], "
            f"\['3', -1.5045949, -3.102543, 1.5979481\], "
            f"\['4', -1.5448789, -3.1964352, 1.6515563\], "
            f"\['1_3', -0.0060455503, -0.85802999, 0.85198444\], "
            f"\['2_4', 0.55469005, 0.092830036, 0.46186001\], "
            f"\['3_4', -0.8973394, -2.3260485, 1.4287091\], "
            f"\['1_2', 0.98670161, 0.96338608, 0.02331553\], "
            f"\['1_2_3', 0.99785712, 0.98674454, 0.01111258\], "
            f"\['1_2_3_4', 0.99807697, 0.98725523, 0.01082174\]\], \n"
            f"The 'r_squared_min' and 'r_squared_rtol' "
            f"variables may need to be adjusted, \n"
            f"there is likely something wrong with the fitting procedure, the "
            f"software parameters need tuned, or there is a bug in the software. \n\n "
            f"NOTE: Since the R-squared values are calculated via different parameters, \n"
            f"the compared R-squared values could be very different if they are not nearly \n"
            f"a perfect fit \(R-squared --> ~0.98 to 0.99999999\)."
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[["CT", "CT", "C", "O_3"]],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.99,
                r_squared_rtol=0.0000002,
            )

    def test_warning_r_squared_min_and_r_squared_rtol_need_adjusted(
        self,
    ):
        with pytest.warns(
            UserWarning,
            match=f"WARNING: The calculated R-squared energy values from the fit types "
            f"do match the validated case for 'r_squared_min' >= "
            f"{'0.98'}, "
            f"but do not fit within the relative tolerance of 'r_squared_rtol' = "
            f"{'0.015'}. \n"
            f"Constants_used = The calculated R-squared energy values from the fit type constants_used. \n"
            f"R-squared_fitted = Gaussian minus GOMC with the selected dihedral set to zero. \n"
            f"R-squared_new_dihedral = Gaussian minus GOMC with the selected individual dihedral added in GOMC. \n"
            f"Abs\(delta\) = Abs\(R-squared_fitted - R-squared_new_dihedral\). \n"
            f"- The fits for all are shown here for all the dihedral combinations \n"
            f"\[opls_constants_used, R-squared_fitted, R-squared_new_dihedral, Abs\(delta\)\] \n"
            f"\[\['1_2', 0.98670161, 0.96338608, 0.02331553\]\]. \n"
            f"The 'r_squared_min' and 'r_squared_rtol' "
            f"variables may need to be adjusted, \n"
            f"there is likely something wrong with the fitting procedure, the "
            f"software parameters need tuned, or there is a bug in the software. \n\n "
            f"NOTE: Since the R-squared values are calculated via different parameters, \n"
            f"the compared R-squared values could be very different if they are not nearly \n"
            f"a perfect fit \(R-squared --> ~0.98 to 0.99999999\)."
        ):
            fit_dihedral_with_gomc(
                ["CT", "CT", "C", "OH"],
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
                ),
                self.get_fn("gmso_oplsaa_CT_CT_C_OH_in_COOH.xml"),
                298.15 * u.Kelvin,
                gomc_binary_directory,
                {
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_2/output"
                    ): [],
                    self.get_fn(
                        "gaussian_style_output_files/CT_CT_C_OH_split_part_1/output"
                    ): [0],
                },
                manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
                zero_dihedral_atom_types=[["CT", "CT", "C", "O_3"]],
                qm_engine="gaussian_style_final_files",
                combining_rule="geometric",
                atom_type_naming_style="general",
                gomc_cpu_cores=1,
                r_squared_min=0.98,
                r_squared_rtol=0.015,
            )

    def test_gaussian_log_file_fit_oplsaa_protonated_fragment_CT_CT_C_OH_in_COOH_in_mie_form(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CT", "CT", "C", "OH"],
            self.get_fn(
                "gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2"
            ),
            self.get_fn("gmso_oplsaa_Mie_style_CT_CT_C_OH_in_COOH.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn("gaussian_style_output_files/CT_CT_C_OH/output"): [
                    0
                ],
            },
            manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
            zero_dihedral_atom_types=[["CT", "CT", "C", "O_3"]],
            qm_engine="gaussian_style_final_files",
            combining_rule="geometric",
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.99,
            r_squared_rtol=0.02,
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 1.80634795241, 0, 0, 0, -0.268742410571],
                    [str("2"), 0, 0, 1.94879786206, 0, 0, 0.448421126946],
                    [str("3"), 0, 0, 0, 1.58083737128, 0, -1.50863674602],
                    [str("4"), 0, 0, 0, 0, 1.5121963942, -1.54542980028],
                    [
                        str("1_3"),
                        0,
                        1.41759854611,
                        0,
                        0.618102768602,
                        0,
                        -0.00379371805915,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.69725316221,
                        0,
                        0.37780223364,
                        0.551946157778,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        0.966769776613,
                        0.903643214082,
                        -0.897967473282,
                    ],
                    [
                        str("1_2"),
                        0,
                        0.882099545744,
                        1.34468116365,
                        0,
                        0,
                        0.986743110285,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        0.83099823616,
                        1.29385707544,
                        0.136793004885,
                        0,
                        0.997926382517,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        0.836836322925,
                        1.29869945716,
                        0.141702758859,
                        -0.0179044297589,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        1.80634795241,
                        -0.903173976205,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.268742410571,
                    ],
                    [
                        str("0_2"),
                        1.94879786206,
                        0,
                        -0.97439893103,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.448421126946,
                    ],
                    [
                        str("0_3"),
                        1.58083737128,
                        0,
                        0,
                        -0.79041868564,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.50863674602,
                    ],
                    [
                        str("0_4"),
                        1.5121963942,
                        0,
                        0,
                        0,
                        -0.7560981971,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.54542980028,
                    ],
                    [
                        str("0_1_3"),
                        2.03570131471,
                        -0.708799273055,
                        0,
                        -0.309051384301,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.00379371805915,
                    ],
                    [
                        str("0_2_4"),
                        2.07505539585,
                        0,
                        -0.848626581105,
                        0,
                        -0.18890111682,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.551946157778,
                    ],
                    [
                        str("0_3_4"),
                        1.8704129907,
                        0,
                        0,
                        -0.483384888306,
                        -0.451821607041,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.897967473282,
                    ],
                    [
                        str("0_1_2"),
                        2.22678070939,
                        -0.441049772872,
                        -0.672340581825,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.986743110285,
                    ],
                    [
                        str("0_1_2_3"),
                        2.26164831649,
                        -0.415499118084,
                        -0.64692853772,
                        -0.0683965024425,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.997926382517,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.25933410919,
                        -0.418418161462,
                        -0.64934972858,
                        -0.0708513794295,
                        0.00895221487945,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        0.903173976205,
                        -0.903173976205,
                        0,
                        0,
                        0,
                        0,
                        -0.268742410571,
                    ],
                    [
                        str("0_2"),
                        1.94879786206,
                        0,
                        -1.94879786206,
                        0,
                        0,
                        0,
                        0.448421126946,
                    ],
                    [
                        str("0_3"),
                        0.79041868564,
                        2.37125605692,
                        0,
                        -3.16167474256,
                        0,
                        0,
                        -1.50863674602,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.0487855768,
                        0,
                        -6.0487855768,
                        0,
                        -1.54542980028,
                    ],
                    [
                        str("0_1_3"),
                        1.01785065736,
                        0.218354879848,
                        0,
                        -1.2362055372,
                        0,
                        0,
                        -0.00379371805915,
                    ],
                    [
                        str("0_2_4"),
                        1.69725316221,
                        0,
                        -0.18604422765,
                        0,
                        -1.51120893456,
                        0,
                        0.551946157778,
                    ],
                    [
                        str("0_3_4"),
                        0.483384888306,
                        1.45015466492,
                        3.61457285633,
                        -1.93353955323,
                        -3.61457285633,
                        0,
                        -0.897967473282,
                    ],
                    [
                        str("0_1_2"),
                        1.78573093652,
                        -0.441049772872,
                        -1.34468116365,
                        0,
                        0,
                        0,
                        0.986743110285,
                    ],
                    [
                        str("0_1_2_3"),
                        1.77775269597,
                        -0.210309610756,
                        -1.29385707544,
                        -0.27358600977,
                        0,
                        0,
                        0.997926382517,
                    ],
                    [
                        str("0_1_2_3_4"),
                        1.78796899805,
                        -0.205864023174,
                        -1.3703171762,
                        -0.283405517718,
                        0.0716177190356,
                        0,
                        0.998103759798,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )
    '''
    '''
    def test_gaussian_style_files_fit_amber_aa_fit_CT_CT_CT_CT_in_butane_files(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CT", "CT", "CT", "CT"],
            self.get_fn(
                "gaussian_style_output_files/CT_CT_CT_CT/input/starting_coords/butane_aa.mol2"
            ),
            self.get_fn("amber_aa_butane_CT_CT_CT_CT_gmso.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_CT_CT/output"
                ): [],
            },
            manual_dihedral_atom_numbers_list=[1, 2, 3, 4],
            zero_dihedral_atom_types=None,
            qm_engine="gaussian_style_final_files",
            combining_rule=None,
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.98,
            r_squared_rtol=0.02,
            opls_force_k0_zero=False
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 2.13492829593, 0, 0, 0, -0.851303269629],
                    [str("2"), 0, 0, 2.25256162343, 0, 0, -0.690454186065],
                    [str("3"), 6.12246736014, 0, 0, -3.06123368007, 0, 0.960991418867],
                    [str("4"), 0, 0, 0, 0, 2.0884153941, -0.91251897505],
                    [
                        str("1_3"),
                        6.5396241564,
                        -0.208579929672,
                        0,
                        -3.06123214853,
                        0,
                        0.977746134654,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.54851017553,
                        0,
                        1.05607698155,
                        -0.497352492204,
                    ],
                    [
                        str("3_4"),
                        6.66816846056,
                        0,
                        0,
                        -3.06123287378,
                        -0.348137574352,
                        0.985623284101,
                    ],
                    [
                        str("1_2"),
                        0,
                        1.13979595936,
                        1.49269914954,
                        0,
                        0,
                        -0.465523847167,
                    ],
                    [
                        str("1_2_3"),
                        6.53963027886,
                        -0.208580993867,
                        0.144337057476,
                        -3.06123414556,
                        0,
                        0.977793543598,
                    ],
                    [
                        str("1_2_3_4"),
                        6.55860366135,
                        -0.208578022741,
                        0.14433854159,
                        -3.06123334065,
                        -0.348136431047,
                        0.982572116976,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        2.1349282959,
                        -1.06746414796,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.851303269629,
                    ],
                    [
                        str("0_2"),
                        2.25256162343,
                        0,
                         -1.12628081172,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.690454186065,
                    ],
                    [
                        str("0_3"),
                        0,
                        0,
                        0,
                        1.53061684004,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.960991418867,
                    ],
                    [
                        str("0_4"),
                        2.0884153941,
                        0,
                        0,
                        0,
                        -1.04420769705,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.91251897505,
                    ],
                    [
                        str("0_1_3"),
                        -1.99973371195e-12,
                        0.104289964836,
                        0,
                        1.53061607426,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.977746134654,
                    ],
                    [
                        str("0_2_4"),
                        2.60458715708 ,
                        0,
                        -0.774255087765,
                        0,
                        -0.528038490775,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.497352492204,
                    ],
                    [
                        str("0_3_4"),
                        -0.075286217852,
                        0,
                        0,
                        1.53061643689,
                        0.174068787176,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.985623284101,
                    ],
                    [
                        str("0_1_2"),
                        2.6324951089,
                        -0.56989797968,
                        -0.74634957477,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.465523847167,
                    ],
                    [
                        str("0_1_2_3"),
                        0.144337057479,
                        0.104290496934,
                        -0.072168528738,
                        1.53061707278,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.977793543598,
                    ],
                    [
                        str("0_1_2_3_4"),
                        -0.194307422173,
                        0.10428901137,
                        -0.072169270795,
                        1.53061667032,
                        0.174068215524,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.982572116976,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        1.06746414796,
                        -1.06746414796,
                        0,
                        0,
                        0,
                        0,
                        -0.851303269629,
                    ],
                    [
                        str("0_2"),
                        2.25256162343,
                        0,
                        -2.25256162343,
                        0,
                        0,
                        0,
                        -0.690454186065,
                    ],
                    [
                        str("0_3"),
                        1.53061684004,
                        -4.5918505201,
                        0,
                        6.12246736014,
                        0,
                        0,
                        0.960991418867,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        8.3536615764,
                        0,
                        -8.3536615764,
                        0,
                        -0.91251897505 ,
                    ],
                    [
                        str("0_1_3"),
                        1.6349060391,
                        -4.48755825796,
                        0,
                        6.12246429706,
                        0,
                        0,
                        0.977746134654,
                    ],
                    [
                        str("0_2_4"),
                        1.54851017553,
                        0,
                        2.67579775067,
                        0,
                        -4.2243079262,
                        0,
                        -0.497352492204,
                    ],
                    [
                        str("0_3_4"),
                        1.80346779339,
                        -4.59184931067,
                        -1.39255029741,
                        6.12246574756,
                        1.39255029741,
                        0,
                        0.985623284101,
                    ],
                    [
                        str("0_1_2"),
                        2.06259712922,
                        -0.56989797968,
                        -1.49269914954,
                        0,
                        0,
                        0,
                        -0.465523847167,
                    ],
                    [
                        str("0_1_2_3"),
                        1.77924462719,
                        -4.48756072141,
                        -0.144337057476,
                        6.12246829112,
                        0,
                        0,
                        0.977793543598,
                    ],
                    [
                        str("0_1_2_3_4"),
                        1.78873469057,
                        -4.4875609996,
                        -1.53688426578,
                        6.1224666813,
                        1.39254572419,
                        0,
                        0.982572116976,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )
    
    def test_gaussian_style_files_fit_opls_ua_fit_CT_CT_CT_CT_in_butane_files(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CH3", "CH2", "CH2", "CH3"],
            self.get_fn(
                "gaussian_style_output_files/CT_CT_CT_CT/input/starting_coords/butane_aa.mol2"
            ),
            self.get_fn("trappe_ua_butane_CT_CT_CT_CT_gmso.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_CT_CT/output"
                ): [],
            },
            manual_dihedral_atom_numbers_list=[1, 2, 3, 4],
            zero_dihedral_atom_types=None,
            qm_engine="gaussian_style_final_files",
            combining_rule=None,
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.97,
            r_squared_rtol=0.02,
            opls_force_k0_zero=True
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 3.89019135095, 0, 0, 0, -0.15385629089],
                    [str("2"), 0, 0, 3.03580389236, 0, 0, -1.08874124105],
                    [str("3"), 0, 0, 0, 4.57026969264, 0, 0.755185044925],
                    [str("4"), 0, 0, 0, 0, 3.24432472808, -0.881847591138],
                    [
                        str("1_3"),
                        0,
                        1.51800454722,
                        0,
                        3.55826890289,
                        0,
                        0.957438274954,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.57125480195,
                        0,
                        2.19682320759,
                        -0.665155000765,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        4.33329974254,
                        0.355455279474,
                        0.76627477714,
                    ],
                    [
                        str("1_2"),
                        0,
                        3.35937487818,
                        0.796225057519,
                        0,
                        0,
                        -0.098211712848,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        1.81662366124,
                        -0.74655694329,
                        3.85689590923,
                        0,
                        0.998529953567,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        1.84350327627,
                        -0.719677798851,
                        3.88377504842,
                        -0.0940775749305,
                        0.999129219573,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        3.89019135095,
                        -1.94509567548,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.15385629089,
                    ],
                    [
                        str("0_2"),
                        3.03580389236,
                        0,
                        -1.51790194618,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.08874124105,
                    ],
                    [
                        str("0_3"),
                        4.57026969264,
                        0,
                        0,
                        -2.28513484632,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.755185044925,
                    ],
                    [
                        str("0_4"),
                        3.24432472808,
                        0,
                        0,
                        0,
                        -1.62216236404,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.881847591138,
                    ],
                    [
                        str("0_1_3"),
                        5.07627345011,
                        -0.75900227361,
                        0,
                        -1.77913445144,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.957438274954 ,
                    ],
                    [
                        str("0_2_4"),
                        3.76807800954,
                        0,
                        -0.785627400975,
                        0,
                        -1.0984116038,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.665155000765,
                    ],
                    [
                        str("0_3_4"),
                        4.68875502201,
                        0,
                        0,
                        -2.16664987127,
                        -0.177727639737,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.76627477714,
                    ],
                    [
                        str("0_1_2"),
                        4.1555999357,
                        -1.67968743909,
                        -0.39811252876,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.098211712848,
                    ],
                    [
                        str("0_1_2_3"),
                        4.92696262718,
                        -0.90831183062,
                        0.373278471645,
                        -1.92844795462,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998529953567,
                    ],
                    [
                        str("0_1_2_3_4"),
                        4.91352295091,
                        -0.921751638135,
                        0.359838899426,
                        -1.94188752421,
                        0.0470387874652,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.999129219573,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        1.94509567548,
                        -1.94509567548,
                        0,
                        0,
                        0,
                        0,
                        -0.15385629089,
                    ],
                    [
                        str("0_2"),
                        3.03580389236,
                        0,
                        -3.03580389236,
                        0,
                        0,
                        0,
                        -1.08874124105,
                    ],
                    [
                        str("0_3"),
                        2.28513484632,
                        6.85540453896,
                        0,
                        -9.14053938528,
                        0,
                        0,
                        0.755185044925,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        12.9772989123,
                        0,
                        -12.9772989123,
                        0,
                        -0.881847591138,
                    ],
                    [
                        str("0_1_3"),
                        2.53813672505,
                        4.57840108073,
                        0,
                        -7.11653780578,
                        0,
                        0,
                        0.957438274954,
                    ],
                    [
                        str("0_2_4"),
                        1.57125480195,
                        0,
                        7.21603802841,
                        0,
                        -8.78729283036,
                        0,
                        -0.665155000765,
                    ],
                    [
                        str("0_3_4"),
                        2.16664987127,
                        6.49994961381,
                        1.4218211179,
                        -8.66659948508,
                        -1.4218211179,
                        0,
                        0.76627477714,
                    ],
                    [
                        str("0_1_2"),
                        2.47591249661,
                        -1.67968743909,
                        -0.796225057519,
                        0,
                        0,
                        0,
                        -0.098211712848,
                    ],
                    [
                        str("0_1_2_3"),
                        2.09020284194,
                        4.877032033237,
                        0.74655694329,
                        -7.71379181846,
                        0,
                        0,
                        0.998529953567,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.14396136349,
                        4.90391093449,
                        0.343367499129,
                        -7.76755009684,
                        0.376310299722,
                        0,
                        0.999129219573,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

    def test_gaussian_style_files_fit_mia_ua_fit_CT_CT_CT_CT_in_butane_files(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CH3", "CH2", "CH2", "CH3"],
            self.get_fn(
                "gaussian_style_output_files/CT_CT_CT_CT/input/starting_coords/butane_aa.mol2"
            ),
            self.get_fn("mie_ua_butane_CT_CT_CT_CT_gmso.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_CT_CT/output"
                ): [],
            },
            manual_dihedral_atom_numbers_list=[1, 2, 3, 4],
            zero_dihedral_atom_types=None,
            qm_engine="gaussian_style_final_files",
            combining_rule=None,
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.98,
            r_squared_rtol=0.02,
            opls_force_k0_zero=False
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 3.89019135095, 0, 0, 0, -0.15385629089],
                    [str("2"), 0, 0, 3.03580389236, 0, 0, -1.08874124105],
                    [str("3"), 0, 0, 0, 4.57026969264, 0, 0.755185044925],
                    [str("4"), 0, 0, 0, 0, 3.24432472808, -0.881847591138],
                    [
                        str("1_3"),
                        0,
                        1.51800454722,
                        0,
                        3.55826890289,
                        0,
                        0.957438274954,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.57125480195,
                        0,
                        2.19682320759,
                        -0.665155000765,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        4.33329974254,
                        0.355455279474,
                        0.76627477714,
                    ],
                    [
                        str("1_2"),
                        0,
                        3.35937487818,
                        0.796225057519,
                        0,
                        0,
                        -0.098211712848,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        1.81662366124,
                        -0.74655694329,
                        3.85689590923,
                        0,
                        0.998529953567,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        1.84350327627,
                        -0.719677798851,
                        3.88377504842,
                        -0.0940775749305,
                        0.999129219573,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        3.89019135095,
                        -1.94509567548,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.15385629089,
                    ],
                    [
                        str("0_2"),
                        3.03580389236,
                        0,
                        -1.51790194618,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.08874124105,
                    ],
                    [
                        str("0_3"),
                        4.57026969264,
                        0,
                        0,
                        -2.28513484632,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.755185044925,
                    ],
                    [
                        str("0_4"),
                        3.24432472808,
                        0,
                        0,
                        0,
                        -1.62216236404,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.881847591138,
                    ],
                    [
                        str("0_1_3"),
                        5.07627345011,
                        -0.75900227361,
                        0,
                        -1.77913445144,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.957438274954 ,
                    ],
                    [
                        str("0_2_4"),
                        3.76807800954,
                        0,
                        -0.785627400975,
                        0,
                        -1.0984116038,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.665155000765,
                    ],
                    [
                        str("0_3_4"),
                        4.68875502201,
                        0,
                        0,
                        -2.16664987127,
                        -0.177727639737,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.76627477714,
                    ],
                    [
                        str("0_1_2"),
                        4.1555999357,
                        -1.67968743909,
                        -0.39811252876,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.098211712848,
                    ],
                    [
                        str("0_1_2_3"),
                        4.92696262718,
                        -0.90831183062,
                        0.373278471645,
                        -1.92844795462,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998529953567,
                    ],
                    [
                        str("0_1_2_3_4"),
                        4.91352295091,
                        -0.921751638135,
                        0.359838899426,
                        -1.94188752421,
                        0.0470387874652,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.999129219573,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        1.94509567548,
                        -1.94509567548,
                        0,
                        0,
                        0,
                        0,
                        -0.15385629089,
                    ],
                    [
                        str("0_2"),
                        3.03580389236,
                        0,
                        -3.03580389236,
                        0,
                        0,
                        0,
                        -1.08874124105,
                    ],
                    [
                        str("0_3"),
                        2.28513484632,
                        6.85540453896,
                        0,
                        -9.14053938528,
                        0,
                        0,
                        0.755185044925,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        12.9772989123,
                        0,
                        -12.9772989123,
                        0,
                        -0.881847591138,
                    ],
                    [
                        str("0_1_3"),
                        2.53813672505,
                        4.57840108073,
                        0,
                        -7.11653780578,
                        0,
                        0,
                        0.957438274954,
                    ],
                    [
                        str("0_2_4"),
                        1.57125480195,
                        0,
                        7.21603802841,
                        0,
                        -8.78729283036,
                        0,
                        -0.665155000765,
                    ],
                    [
                        str("0_3_4"),
                        2.16664987127,
                        6.49994961381,
                        1.4218211179,
                        -8.66659948508,
                        -1.4218211179,
                        0,
                        0.76627477714,
                    ],
                    [
                        str("0_1_2"),
                        2.47591249661,
                        -1.67968743909,
                        -0.796225057519,
                        0,
                        0,
                        0,
                        -0.098211712848,
                    ],
                    [
                        str("0_1_2_3"),
                        2.09020284194,
                        4.877032033237,
                        0.74655694329,
                        -7.71379181846,
                        0,
                        0,
                        0.998529953567,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.14396136349,
                        4.90391093449,
                        0.343367499129,
                        -7.76755009684,
                        0.376310299722,
                        0,
                        0.999129219573,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

    def test_gaussian_style_files_fit_amber_fit_CT_CT_CT_CT_in_butane_files(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CT", "CT", "CT", "CT"],
            self.get_fn(
                "gaussian_style_output_files/CT_CT_CT_CT/input/starting_coords/butane_aa.mol2"
            ),
            self.get_fn("exp6_aa_butane_CT_CT_CT_CT_gmso.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian_style_output_files/CT_CT_CT_CT/output"
                ): [],
            },
            manual_dihedral_atom_numbers_list=[1, 2, 3, 4],
            zero_dihedral_atom_types=None,
            qm_engine="gaussian_style_final_files",
            combining_rule=None,
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.98,
            r_squared_rtol=0.02,
            opls_force_k0_zero=False
        )

        assert (
            os.path.isfile("all_normalized_energies_in_kcal_per_mol.txt")
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt"
            )
            is True
        )
        assert (
            os.path.isfile(
                "all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt"
            )
            is True
        )
        assert os.path.isfile("gomc_raw_energies_in_Kelvin.txt") is True
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt")
            is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile("gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt") is True
        )
        assert (
            os.path.isfile(
                "opls_all_single_fit_dihedral_k_constants_figure.pdf"
            )
            is True
        )
        assert (
            os.path.isfile("opls_all_summed_dihedrals_k_constants_figure.pdf")
            is True
        )
        assert (
            os.path.isfile("opls_dihedral_k_constants_fit_energy.txt") is True
        )
        assert (
            os.path.isfile("periodic_dihedral_k_constants_fit_energy.txt")
            is True
        )
        assert os.path.isfile("RB_torsion_k_constants_fit_energy.txt") is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "r_squared",
                    ],
                    [str("1"), 0, 4.64189989445, 0, 0, 0, 0.0554910823082],
                    [str("2"), 0, 0, 3.62275916028, 0, 0, -1.23591651342],
                    [str("3"), 0, 0, 0, 4.92948888627, 0, 0.477550076455],
                    [str("4"), 0, 0, 0, 0, 3.7172458491, -1.12957843414],
                    [
                        str("1_3"),
                        0,
                        2.4400162262,
                        0,
                        3.30281500888,
                        0,
                        0.9846783871,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        2.06026818365,
                        0,
                        2.34373594682,
                        -0.768017540175,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        4.41238953459,
                        0.775649812101,
                        0.52879674349,
                    ],
                    [
                        str("1_2"),
                        0,
                        4.00810304165,
                        0.950695705123,
                        0,
                        0,
                        0.132478073226,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        2.61642006674,
                        -0.441015014606,
                        3.47922351339,
                        0,
                        0.998594497652,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        2.65142528583,
                        -0.406010416441,
                        3.51422811225,
                        -0.122516855939,
                        0.999580827406,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "n0_kcal_per_mol",
                        "n1_kcal_per_mol",
                        "n2_kcal_per_mol",
                        "n3_kcal_per_mol",
                        "n4_kcal_per_mol",
                        "n5_kcal_per_mol",
                        "d0_kcal_per_mol",
                        "d1_kcal_per_mol",
                        "d2_kcal_per_mol",
                        "d3_kcal_per_mol",
                        "d4_kcal_per_mol",
                        "d5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        4.64189989445,
                        -2.32094994722,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.0554910823082,
                    ],
                    [
                        str("0_2"),
                        3.62275916028,
                        0,
                        -1.81137958014,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.23591651342,
                    ],
                    [
                        str("0_3"),
                        4.92948888627,
                        0,
                        0,
                        -2.46474444314,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.477550076455,
                    ],
                    [
                        str("0_4"),
                        3.7172458491,
                        0,
                        0,
                        0,
                        -1.85862292455,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -1.12957843414,
                    ],
                    [
                        str("0_1_3"),
                        5.74283123508,
                        -1.2200081131,
                        0,
                        -1.65140750444,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                         0.9846783871,
                    ],
                    [
                        str("0_2_4"),
                        4.40400413047,
                        0,
                        -1.03013409182,
                        0,
                        -1.17186797341,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        -0.768017540175,
                    ],
                    [
                        str("0_3_4"),
                        5.18803934669,
                        0,
                        0,
                        -2.2061947673,
                        -0.38782490605,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.52879674349,
                    ],
                    [
                        str("0_1_2"),
                        4.95879874677,
                        -2.00405152082,
                        -0.475347852562,
                        0,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.132478073226,
                    ],
                    [
                        str("0_1_2_3"),
                        5.65462856552,
                        -1.3082100333,
                        0.220507507303,
                        -1.7396117567,
                        0,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.998594497652,
                    ],
                    [
                        str("0_1_2_3_4"),
                        5.6371261257,
                        -1.32571264291,
                        0.20300520822,
                        -1.75711405612,
                        0.0612584279695,
                        0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0,
                        180.0,
                        0,
                        180.0,
                        0.999580827406,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )

        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        "non_zero_k_constants",
                        "k0_kcal_per_mol",
                        "k1_kcal_per_mol",
                        "k2_kcal_per_mol",
                        "k3_kcal_per_mol",
                        "k4_kcal_per_mol",
                        "k5_kcal_per_mol",
                        "r_squared",
                    ],
                    [
                        str("0_1"),
                        2.32094994722,
                        -2.32094994722,
                        0,
                        0,
                        0,
                        0,
                        0.0554910823082,
                    ],
                    [
                        str("0_2"),
                        3.62275916028,
                        0,
                        -3.62275916028,
                        0,
                        0,
                        0,
                        -1.23591651342,
                    ],
                    [
                        str("0_3"),
                        2.46474444314,
                        7.3942333294,
                        0,
                        -9.85897777254,
                        0,
                        0,
                        0.477550076455,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        14.8689833964,
                        0,
                        -14.8689833964,
                        0,
                        -1.12957843414,
                    ],
                    [
                        str("0_1_3"),
                        2.87141561754,
                        3.73421440022,
                        0,
                        -6.60563001776,
                        0,
                        0,
                        0.9846783871,
                    ],
                    [
                        str("0_2_4"),
                        2.06026818365,
                        0,
                        7.31467560363,
                        0,
                        -9.37494378728,
                        0,
                        -0.768017540175,
                    ],
                    [
                        str("0_3_4"),
                        2.2061947673,
                        6.61858430188,
                        3.1025992484,
                        -8.82477906918,
                        -3.1025992484,
                        0,
                        0.52879674349,
                    ],
                    [
                        str("0_1_2"),
                        2.95474722595,
                        -2.00405152082,
                        -0.950695705123,
                        0,
                        0,
                        0,
                        0.132478073226,
                    ],
                    [
                        str("0_1_2_3"),
                        2.60680677546,
                        3.91062523672,
                        0.441015014606,
                        -6.95844702678,
                        0,
                        0,
                        0.998594497652,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.6768162826,
                        3.94562952546,
                        -0.084057007315,
                        -7.0284562245,
                        0.490067423756,
                        0,
                        0.999580827406,
                    ],
                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(
                                correct_line_values[i][j]
                            )

                        # check the k-values and the r-squared fit
                        else:
                            assert np.isclose(
                                mdf_math.round_to_sig_figs(
                                    float(split_line_i[j]), number_sig_i
                                ),
                                mdf_math.round_to_sig_figs(
                                    correct_line_values[i][j], number_sig_i
                                ),
                                atol=0.02,
                                rtol=0.08,
                            )