import os
import re
import subprocess

import numpy as np
import pytest
import unyt as u

import mosdef_dihedral_fit.utils.math_operations as mdf_math
from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import (
    fit_dihedral_with_gomc,
)
from mosdef_dihedral_fit.tests.base_test import BaseTest

# The below "gomc_binary_directory = "/opt/gomc_build/GOMC/bin"" is the
# automated GitHub test binary path.
# The automated GitHub test filename and location = "MoSDeF-dihedral-fit/.github/workflows/CI.yml"
# set GOMC binary path here if you want to manually specify it
gomc_binary_directory = ""
# User changable variable, if the binary is not in the environment variable "PATH"
# Examples
# gomc_binary_directory = "/home/brad/Programs/GOMC/GOMC-2.75a/bin"
# gomc_binary_directory = "/Users/calcraven/Documents/Vanderbilt/Research/MoSDeF/Dihedral_Fitter/GOMC/bin"
if gomc_binary_directory == "":
    gomc_binary_directory = subprocess.check_output(
        "echo ${CONDA_PREFIX}/bin", shell=True, text=True
    ).strip()
if not gomc_binary_directory or "GOMC_CPU_NVT" not in os.listdir(
    gomc_binary_directory
):
    raise OSError(
        "Missing GOMC installation. Please install from https://github.com/GOMC-WSU/GOMC and add to conda environment `$CONDA_PREFIX`, or manually set gomc_binary_directory variable in tests/test_fit_dihedral_with_gomc.py"
    )


# NOTE: When comparing fitted values with reference value,
# we are using numpy.isclose() with absolute tolerance of
# 0.02 and absolute tolerance of 0.08 (8%) to account for
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
            r_squared_atol=1e-03,
            opls_force_k0_zero=True,
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
            r_squared_atol=1e-03,
            opls_force_k0_zero=True,
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
            r_squared_atol=1e-03,
            opls_force_k0_zero=True,
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
                    [str("1"), 0, 2.0043843779, 0, 0, 0, -1],
                    [str("2"), 0, 0, 2.13887995673, 0, 0, 0.265058829483],
                    [str("3"), 0, 0, 0, 1.79266434371, 0, -1],
                    [str("4"), 0, 0, 0, 0, 1.71251524752, -1],
                    [
                        str("1_3"),
                        0,
                        1.45669370251,
                        0,
                        0.821535056147,
                        0,
                        0.0566682712428,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.79496359399,
                        0,
                        0.515872948092,
                        0.434086151613,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        1.17177915373,
                        0.931327790881,
                        -1,
                    ],
                    [
                        str("1_2"),
                        0,
                        1.04124236441,
                        1.44471810818,
                        0,
                        0,
                        0.953675473953,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        0.925218215852,
                        1.32869479285,
                        0.29005878853,
                        0,
                        0.998562878637,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        0.914195508542,
                        1.31767182551,
                        0.279035947504,
                        0.0385798985516,
                        0.99929215065,
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
                        2.0043843779,
                        -1.00219218895,
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
                        -1,
                    ],
                    [
                        str("0_2"),
                        2.13887995673,
                        0,
                        -1.06943997836,
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
                        0.265058829483,
                    ],
                    [
                        str("0_3"),
                        1.79266434371,
                        0,
                        0,
                        -0.896332171855,
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
                        -1,
                    ],
                    [
                        str("0_4"),
                        1.71251524752,
                        0,
                        0,
                        0,
                        -0.85625762376,
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
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        2.27822875866,
                        -0.728346851255,
                        0,
                        -0.410767528074,
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
                        0.0566682712428,
                    ],
                    [
                        str("0_2_4"),
                        2.31083654208,
                        0,
                        -0.897481796995,
                        0,
                        -0.257936474046,
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
                        0.434086151613,
                    ],
                    [
                        str("0_3_4"),
                        2.10310694461,
                        0,
                        0,
                        -0.585889576865,
                        -0.46566389544,
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
                        -1,
                    ],
                    [
                        str("0_1_2"),
                        2.48596047259,
                        -0.520621182205,
                        -0.72235905409,
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
                        0.953675473953,
                    ],
                    [
                        str("0_1_2_3"),
                        2.54397179723,
                        -0.462609107926,
                        -0.664347396425,
                        -0.145029394265,
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
                        0.998562878637,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.54948318011,
                        -0.457097754271,
                        -0.658835912755,
                        -0.139517973752,
                        -0.0192899492758,
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
                        0.99929215065,
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
                        1.00219218895,
                        -1.00219218895,
                        0,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_2"),
                        2.13887995673,
                        0,
                        -2.13887995673,
                        0,
                        0,
                        0,
                        0.265058829483,
                    ],
                    [
                        str("0_3"),
                        0.896332171855,
                        2.68899651556,
                        0,
                        -3.58532868742,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.85006099008,
                        0,
                        -6.85006099008,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        1.13911437933,
                        0.503955732966,
                        0,
                        -1.64307011229,
                        0,
                        0,
                        0.0566682712428,
                    ],
                    [
                        str("0_2_4"),
                        1.79496359399,
                        0,
                        0.268528198378,
                        0,
                        -2.06349179237,
                        0,
                        0.434086151613,
                    ],
                    [
                        str("0_3_4"),
                        0.585889576865,
                        1.75766873059,
                        3.72531116352,
                        -2.34355830746,
                        -3.72531116352,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_2"),
                        1.96533929038,
                        -0.520621182205,
                        -1.44471810818,
                        0,
                        0,
                        0,
                        0.953675473953,
                    ],
                    [
                        str("0_1_2_3"),
                        1.93633329504,
                        -0.027520925131,
                        -1.32869479285,
                        -0.58011757706,
                        0,
                        0,
                        0.998562878637,
                    ],
                    [
                        str("0_1_2_3_4"),
                        1.91428755353,
                        -0.038543833015,
                        -1.1633522313,
                        -0.558071895008,
                        -0.154319594206,
                        0,
                        0.99929215065,
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
            r_squared_atol=0.02,
            opls_force_k0_zero=True,
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
                    [str("1"), 0, 1.80634795241, 0, 0, 0, -1],
                    [str("2"), 0, 0, 1.94879786206, 0, 0, 0.448421126946],
                    [str("3"), 0, 0, 0, 1.58083737128, 0, -1],
                    [str("4"), 0, 0, 0, 0, 1.5121963942, -1],
                    [
                        str("1_3"),
                        0,
                        1.41759854611,
                        0,
                        0.618102768602,
                        0,
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.0487855768,
                        0,
                        -6.0487855768,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        1.01785065736,
                        0.218354879848,
                        0,
                        -1.2362055372,
                        0,
                        0,
                        -1,
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
                        -1,
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
            r_squared_atol=5e-03,
            opls_force_k0_zero=True,
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
                    [str("1"), 0, 1.87103578693, 0, 0, 0, -1],
                    [str("2"), 0, 0, 2.01338731788, 0, 0, 0.468487865197],
                    [str("3"), 0, 0, 0, 1.64674772188, 0, -1],
                    [str("4"), 0, 0, 0, 0, 1.5768734149, -1],
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.3074936596,
                        0,
                        -6.3074936596,
                        0,
                        -1,
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
                        -1,
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
            r_squared_atol=0.02,
            opls_force_k0_zero=True,
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
                    [str("1"), 0, 1.80634795241, 0, 0, 0, -1],
                    [str("2"), 0, 0, 1.94879786206, 0, 0, 0.448421126946],
                    [str("3"), 0, 0, 0, 1.58083737128, 0, -1],
                    [str("4"), 0, 0, 0, 0, 1.5121963942, -1],
                    [
                        str("1_3"),
                        0,
                        1.41759854611,
                        0,
                        0.618102768602,
                        0,
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        6.0487855768,
                        0,
                        -6.0487855768,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        1.01785065736,
                        0.218354879848,
                        0,
                        -1.2362055372,
                        0,
                        0,
                        -1,
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
                        -1,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=5e-03,
                opls_force_k0_zero=True,
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
            r_squared_atol=1e-03,
            opls_force_k0_zero=True,
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
            r_squared_atol=1e-03,
            opls_force_k0_zero=True,
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
            r_squared_atol=1e-03,
            opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
            )

    def test_bad_fit_opls_force_k0_zero_not_bool(self):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the 'opls_force_k0_zero' as a bool.",
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
                r_squared_atol=1e-03,
                opls_force_k0_zero="x",
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
            )

    def test_mol2_file_file_no_mol2_extention(self):
        value_path_mol2 = "bad_mol2_path"
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please enter mol2 file \('mol2_file'\) name with the .mol2 extension.",
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
            )

    def test_xml_selection_file_no_xml_extention(self):
        value_path_xml = "bad_xml_path"
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please enter xml file "
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
            )

    def test_gomc_binary_path_containing_the_GOMC_CPU_NVT_file_does_not_exist(
        self,
    ):
        errormsg = r"ERROR: The 'gomc_binary_path' of gomc_binary_path does not exist or contain the GOMC 'GOMC_CPU_NVT' file. If the path is incorrect, please pass the correct absolute path for gomc_binary_path. If GOMC is not installed, go to https://github.com/GOMC-WSU/MoSDeF-dihedral-fit/blob/main/docs/getting_started/installation/installation.rst#install-gomc README.md for this repository."
        with pytest.raises(OSError, match=re.escape(errormsg)):
            fit_dihedral_with_gomc(
                ["HC", "CT", "CT", "HC"],
                self.get_fn(
                    "gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2"
                ),
                self.get_fn("oplsaa_ethane_HC_CT_CT_HC.xml"),
                298.15 * u.Kelvin,
                f"gomc_binary_path",
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
                r_squared_atol=1e-03,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
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
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
            )

    def test_r_squared_min_not_a_float(self):
        with pytest.raises(
            TypeError,
            match=f"ERROR: The 'r_squared_min' is a {type(1)}, "
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
                r_squared_min=1,
                r_squared_atol=0.02,
                opls_force_k0_zero=True,
            )

    def test_r_squared_atol_not_correct_value_is_0(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'r_squared_atol' = {0.00}, "
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
                r_squared_atol=0.00,
                opls_force_k0_zero=True,
            )

    def test_r_squared_atol_not_correct_value_is_1(self):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The 'r_squared_atol' = {1.00}, "
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
                r_squared_atol=1.00,
                opls_force_k0_zero=True,
            )

    def test_r_squared_atol_not_a_float(self):
        with pytest.raises(
            TypeError,
            match=f"ERROR: The 'r_squared_atol' is a {type(1)}, "
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
                r_squared_atol=1,
                opls_force_k0_zero=True,
            )

    def test_error_r_squared_min_and_r_squared_atol_need_adjusted(
        self,
    ):
        with pytest.raises(
            ValueError,
            match=f"ERROR: The calculated R-squared energy values from the fit types "
            f"do not match the validated case for 'r_squared_min' >= "
            f"{'0.99'}, "
            f"within the absolute tolerance or 'r_squared_atol' = "
            f"{'1e-07'}. \n"
            f"Constants_used = The calculated R-squared energy values from the fit type constants_used."
            f"R-squared_fitted = Gaussian minus GOMC with the selected dihedral set to zero \n"
            f"R-squared_new_dihedral = Gaussian minus GOMC with the selected individual dihedral added in GOMC \n"
            f"Abs\(delta\) = absolute_value\(R-squared_fitted - R-squared_new_dihedral\) \n"
            f"- The fits for all are shown here for all the dihedral combinations \n"
            f"\[opls_constants_used, R-squared_fitted, R-squared_new_dihedral, Abs\(delta\)\] \n"
            f"\[\['1', -1.0, -1.0, 0\], "
            f"\['2', 0.44851424, -1.0, 1.4485142\], "
            f"\['3', -1.0, -1.0, 0\], "
            f"\['4', -1.0, -1.0, 0\], "
            f"\['1_3', -1.0, -1.0, 0\], "
            f"\['2_4', 0.55202043, -1.0, 1.5520204\], "
            f"\['3_4', -1.0, -1.0, 0\], "
            f"\['1_2', 0.98683182, 0.96185841, 0.02497341\], "
            f"\['1_2_3', 0.99796071, 0.9861541, 0.01180661\], "
            f"\['1_2_3_4', 0.99813699, 0.98663398, 0.01150301\]\], \n"
            f"The 'r_squared_min' and 'r_squared_atol' "
            f"variables may need to be adjusted, \n"
            f"there is likely something wrong with the fitting procedure, the "
            f"software parameters need tuned, or there is a bug in the software. \n\n "
            f"NOTE: Since the R-squared values are calculated via different parameters, \n"
            f"the compared R-squared values could be very different if they are not nearly \n"
            f"a perfect fit \(R-squared --> ~0.98 to 0.99999999\).",
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
                r_squared_atol=0.0000001,
            )

    def test_warning_r_squared_min_and_r_squared_atol_need_adjusted(
        self,
    ):
        with pytest.warns(
            UserWarning,
            match=f"WARNING: The calculated R-squared energy values from the fit types "
            f"do match the validated case for 'r_squared_min' >= "
            f"{'0.95'}, "
            f"but do not fit within the absolute tolerance of 'r_squared_atol' = "
            f"{'0.02'}. \n"
            f"Constants_used = The calculated R-squared energy values from the fit type constants_used. \n"
            f"R-squared_fitted = Gaussian minus GOMC with the selected dihedral set to zero. \n"
            f"R-squared_new_dihedral = Gaussian minus GOMC with the selected individual dihedral added in GOMC. \n"
            f"Abs\(delta\) = Abs\(R-squared_fitted - R-squared_new_dihedral\). \n"
            f"- The fits for all are shown here for all the dihedral combinations \n"
            f"\[opls_constants_used, R-squared_fitted, R-squared_new_dihedral, Abs\(delta\)\] \n"
            f"\[\['1_2', 0.98683182, 0.96185841, 0.02497341\]\]. \n"
            f"The 'r_squared_min' and 'r_squared_atol' "
            f"variables may need to be adjusted, \n"
            f"there is likely something wrong with the fitting procedure, the "
            f"software parameters need tuned, or there is a bug in the software. \n\n "
            f"NOTE: Since the R-squared values are calculated via different parameters, \n"
            f"the compared R-squared values could be very different if they are not nearly \n"
            f"a perfect fit \(R-squared --> ~0.98 to 0.99999999\).",
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
                r_squared_min=0.95,
                r_squared_atol=0.02,
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
            r_squared_atol=0.02,
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
                    [str("1"), 0, 0.840869276724, 0, 0, 0, -1],
                    [str("2"), 0, 0, 1.94841848488, 0, 0, 0.448514120427],
                    [str("3"), 0, 0, 0, 0.0278527540281, 0, -1],
                    [str("4"), 0.188139728024, 0, 0, 0, -0.0940698640118, -1],
                    [
                        str("1_3"),
                        0,
                        0.845246869908,
                        0,
                        0.0762104684416,
                        0,
                        -1,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        1.69688361533,
                        0,
                        0.377787471995,
                        0.552020352149,
                    ],
                    [
                        str("3_4"),
                        0.179654635283,
                        0,
                        0,
                        0.0227319855836,
                        -0.0928118945261,
                        -1,
                    ],
                    [
                        str("1_2"),
                        0,
                        0.8821416593,
                        1.34427294497,
                        0,
                        0,
                        0.986831815166,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        0.831162091238,
                        1.29356993782,
                        0.136467115148,
                        0,
                        0.997926382517,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        0.836982369298,
                        1.29839754787,
                        0.141361890758,
                        -0.0178498115248,
                        0.998136993063,
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
                        0.840869276724,
                        -0.420434638362,
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
                        -1,
                    ],
                    [
                        str("0_2"),
                        1.94841848488,
                        0,
                        -0.97420924244,
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
                        0.448514120427,
                    ],
                    [
                        str("0_3"),
                        0.0278527540281,
                        0,
                        0,
                        -0.013926377014,
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
                        -1,
                    ],
                    [
                        str("0_4"),
                        2.00006677886e-13,
                        0,
                        0,
                        0,
                        0.0470349320059,
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
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        0.92145733835,
                        -0.422623434954,
                        0,
                        -0.0381052342208,
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
                        -1,
                    ],
                    [
                        str("0_2_4"),
                        2.07467108732,
                        0,
                        -0.848441807665,
                        0,
                        -0.188893735998,
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
                        0.552020352149,
                    ],
                    [
                        str("0_3_4"),
                        0.019747408699,
                        0,
                        0,
                        -0.0113659927918,
                        0.046405947263,
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
                        -1,
                    ],
                    [
                        str("0_1_2"),
                        2.22641460427,
                        -0.44107082965,
                        -0.672136472485,
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
                        0.986831815166,
                    ],
                    [
                        str("0_1_2_3"),
                        2.26119914421,
                        -0.415581045619,
                        -0.64678496891,
                        -0.068233557574,
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
                        2.2588919964,
                        -0.418491184649,
                        -0.649198773935,
                        -0.070680945379,
                        0.0089249057624,
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
                        0.998136993063,
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
                        0.420434638362,
                        -0.420434638362,
                        0,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_2"),
                        1.94841848488,
                        0,
                        -1.94841848488,
                        0,
                        0,
                        0,
                        0.448514120427,
                    ],
                    [
                        str("0_3"),
                        0.013926377014,
                        0.0417791310421,
                        0,
                        -0.0557055080562,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_4"),
                        0.094069864012,
                        0,
                        -0.376279456047,
                        0,
                        0.376279456047,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        0.460728669175,
                        -0.308307732292,
                        0,
                        -0.152420936883,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_2_4"),
                        1.69688361533,
                        0,
                        -0.18573372735,
                        0,
                        -1.51114988798,
                        0,
                        0.552020352149,
                    ],
                    [
                        str("0_3_4"),
                        0.101193310433,
                        0.0340979783754,
                        -0.371247578104,
                        -0.0454639711672,
                        0.371247578104,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_2"),
                        1.78534377462,
                        -0.44107082965,
                        -1.34427294497,
                        0,
                        0,
                        0,
                        0.986831815166,
                    ],
                    [
                        str("0_1_2_3"),
                        1.77738454101,
                        -0.210880372897,
                        -1.29356993782,
                        -0.272934230296,
                        0,
                        0,
                        0.997960714354,
                    ],
                    [
                        str("0_1_2_3_4"),
                        1.7875696779,
                        -0.206448348512,
                        -1.36979679397,
                        -0.282723781516,
                        0.0713992460992,
                        0,
                        0.998136993063,
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
            r_squared_atol=0.02,
            opls_force_k0_zero=False,
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
                    [str("1"), 0.269321156456, -0.134660578228, 0, 0, 0, -1],
                    [str("2"), 0, 0, 0.173802523735, 0, 0, -1],
                    [str("3"), 0, 0, 0, 0.68578301424, 0, 0.538310820492],
                    [
                        str("4"),
                        0.549354203082,
                        0,
                        0,
                        0,
                        -0.274677101541,
                        -1,
                    ],
                    [
                        str("1_3"),
                        0,
                        0.0268690085126,
                        0,
                        0.667870382182,
                        0,
                        0.541423747235,
                    ],
                    [
                        str("2_4"),
                        0.389167581037,
                        0,
                        0.173803714089,
                        0,
                        -0.2746778490549,
                        -1,
                    ],
                    [
                        str("3_4"),
                        0.473988099803,
                        0,
                        0,
                        0.506335034115,
                        -0.274678271563,
                        0.833217913109,
                    ],
                    [
                        str("1_2"),
                        0.269323740886,
                        -0.134661870443,
                        0.173803524331,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        -0.0689389848116,
                        0.23952293664,
                        0.572059846546,
                        0,
                        0.74922032487,
                    ],
                    [
                        str("1_2_3_4"),
                        0.482200556331,
                        -0.134663245963,
                        0.173797700249,
                        0.506333618372,
                        -0.274677871599,
                        0.980910089059,
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
                        0,
                        0.067330289114,
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
                        -1,
                    ],
                    [
                        str("0_2"),
                        0.173802523735,
                        0,
                        -0.0869012618675,
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
                        -1,
                    ],
                    [
                        str("0_3"),
                        0.68578301424,
                        0,
                        0,
                        -0.34289150712,
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
                        0.538310820492,
                    ],
                    [
                        str("0_4"),
                        4.9998893914e-13,
                        0,
                        0,
                        0,
                        0.13733855077,
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
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        0.694739390695,
                        -0.0134345042563,
                        0,
                        -0.333935191091,
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
                        0.541423747235,
                    ],
                    [
                        str("0_2_4"),
                        0.0937096555535,
                        0,
                        -0.0869018570445,
                        0,
                        0.137338924527,
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
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        0.468650812454,
                        0,
                        0,
                        -0.253167517058,
                        0.137339135782,
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
                        0.833217913109,
                    ],
                    [
                        str("0_1_2"),
                        0.173803524331,
                        0.0673309352215,
                        -0.0869017621655,
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
                        -1,
                    ],
                    [
                        str("0_1_2_3"),
                        0.742643798374,
                        0.0344694924058,
                        -0.11976146832,
                        -0.286029923273,
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
                        0.74922032487,
                    ],
                    [
                        str("0_1_2_3_4"),
                        0.511890479224,
                        0.0673316229815,
                        -0.0868988501245,
                        -0.253166809186,
                        0.137338935799,
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
                        0.980910089059,
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
                        0.067330289114,
                        0.067330289114,
                        0,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_2"),
                        0.173802523735,
                        0,
                        -0.173802523735,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_3"),
                        0.34289150712,
                        1.02867452136,
                        0,
                        -1.37156602848,
                        0,
                        0,
                        0.538310820492,
                    ],
                    [
                        str("0_4"),
                        0.274677101541,
                        0,
                        -1.09870840616,
                        0,
                        1.09870840616,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        0.347369695347,
                        0.988371069017,
                        0,
                        -1.33574076436,
                        0,
                        0,
                        0.541423747235,
                    ],
                    [
                        str("0_2_4"),
                        0.368387504608,
                        0,
                        -1.27251511031,
                        0,
                        1.09871139622,
                        0,
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        0.490161566959,
                        0.759502551172,
                        -1.09871308625,
                        -1.01267006823,
                        1.09871308625,
                        0,
                        0.833217913109,
                    ],
                    [
                        str("0_1_2"),
                        0.241134459552,
                        0.0673309352215,
                        -0.173803524331,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_2_3"),
                        0.491083367507,
                        0.892559262225,
                        -0.23952293664,
                        -1.144119693092,
                        0,
                        0,
                        0.74922032487,
                    ],
                    [
                        str("0_1_2_3_4"),
                        0.600733164619,
                        0.826832050539,
                        -1.27250918664,
                        -1.01266723674,
                        1.0987114864,
                        0,
                        0.980910089059,
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

    def test_gaussian_style_files_fit_amber_aa_all_unique_fit_CT_CT_CT_CT_in_butane_files(
        self,
    ):
        fit_dihedral_with_gomc(
            ["CT1", "CT0", "CT0", "CT1"],
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
            atom_type_naming_style="all_unique",
            gomc_cpu_cores=1,
            r_squared_min=0.98,
            r_squared_atol=0.02,
            opls_force_k0_zero=False,
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
                    [str("1"), 0.269321156456, -0.134660578228, 0, 0, 0, -1],
                    [str("2"), 0, 0, 0.173802523735, 0, 0, -1],
                    [str("3"), 0, 0, 0, 0.68578301424, 0, 0.538310820492],
                    [
                        str("4"),
                        0.549354203082,
                        0,
                        0,
                        0,
                        -0.274677101541,
                        -1,
                    ],
                    [
                        str("1_3"),
                        0,
                        0.0268690085126,
                        0,
                        0.667870382182,
                        0,
                        0.541423747235,
                    ],
                    [
                        str("2_4"),
                        0.389167581037,
                        0,
                        0.173803714089,
                        0,
                        -0.2746778490549,
                        -1,
                    ],
                    [
                        str("3_4"),
                        0.473988099803,
                        0,
                        0,
                        0.506335034115,
                        -0.274678271563,
                        0.833217913109,
                    ],
                    [
                        str("1_2"),
                        0.269323740886,
                        -0.134661870443,
                        0.173803524331,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        -0.0689389848116,
                        0.23952293664,
                        0.572059846546,
                        0,
                        0.74922032487,
                    ],
                    [
                        str("1_2_3_4"),
                        0.482200556331,
                        -0.134663245963,
                        0.173797700249,
                        0.506333618372,
                        -0.274677871599,
                        0.980910089059,
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
                        0,
                        0.067330289114,
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
                        -1,
                    ],
                    [
                        str("0_2"),
                        0.173802523735,
                        0,
                        -0.0869012618675,
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
                        -1,
                    ],
                    [
                        str("0_3"),
                        0.68578301424,
                        0,
                        0,
                        -0.34289150712,
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
                        0.538310820492,
                    ],
                    [
                        str("0_4"),
                        4.9998893914e-13,
                        0,
                        0,
                        0,
                        0.13733855077,
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
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        0.694739390695,
                        -0.0134345042563,
                        0,
                        -0.333935191091,
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
                        0.541423747235,
                    ],
                    [
                        str("0_2_4"),
                        0.0937096555535,
                        0,
                        -0.0869018570445,
                        0,
                        0.137338924527,
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
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        0.468650812454,
                        0,
                        0,
                        -0.253167517058,
                        0.137339135782,
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
                        0.833217913109,
                    ],
                    [
                        str("0_1_2"),
                        0.173803524331,
                        0.0673309352215,
                        -0.0869017621655,
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
                        -1,
                    ],
                    [
                        str("0_1_2_3"),
                        0.742643798374,
                        0.0344694924058,
                        -0.11976146832,
                        -0.286029923273,
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
                        0.74922032487,
                    ],
                    [
                        str("0_1_2_3_4"),
                        0.511890479224,
                        0.0673316229815,
                        -0.0868988501245,
                        -0.253166809186,
                        0.137338935799,
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
                        0.980910089059,
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
                        0.067330289114,
                        0.067330289114,
                        0,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_2"),
                        0.173802523735,
                        0,
                        -0.173802523735,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_3"),
                        0.34289150712,
                        1.02867452136,
                        0,
                        -1.37156602848,
                        0,
                        0,
                        0.538310820492,
                    ],
                    [
                        str("0_4"),
                        0.274677101541,
                        0,
                        -1.09870840616,
                        0,
                        1.09870840616,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        0.347369695347,
                        0.988371069017,
                        0,
                        -1.33574076436,
                        0,
                        0,
                        0.541423747235,
                    ],
                    [
                        str("0_2_4"),
                        0.368387504608,
                        0,
                        -1.27251511031,
                        0,
                        1.09871139622,
                        0,
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        0.490161566959,
                        0.759502551172,
                        -1.09871308625,
                        -1.01267006823,
                        1.09871308625,
                        0,
                        0.833217913109,
                    ],
                    [
                        str("0_1_2"),
                        0.241134459552,
                        0.0673309352215,
                        -0.173803524331,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_2_3"),
                        0.491083367507,
                        0.892559262225,
                        -0.23952293664,
                        -1.144119693092,
                        0,
                        0,
                        0.74922032487,
                    ],
                    [
                        str("0_1_2_3_4"),
                        0.600733164619,
                        0.826832050539,
                        -1.27250918664,
                        -1.01266723674,
                        1.0987114864,
                        0,
                        0.980910089059,
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

    def test_gaussian_style_files_fit_amber_aa_fit_CT_CT_CT_CT_force_k0_to_True_in_butane_files(
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
            r_squared_min=0.8,
            r_squared_atol=0.2,
            opls_force_k0_zero=True,
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
                    [str("1"), 0, 0.4721804016, 0, 0, 0, -1],
                    [str("2"), 0, 0, 0.575010099536, 0, 0, -1],
                    [str("3"), 0, 0, 0, 0.685865110811, 0, 0.538213945525],
                    [str("4"), 0, 0, 0, 0, 0.425529576141, -1],
                    [
                        str("1_3"),
                        0,
                        0.0268839930837,
                        0,
                        0.667942488722,
                        0,
                        0.541330329235,
                    ],
                    [
                        str("2_4"),
                        0,
                        0,
                        0.524381976042,
                        0,
                        0.075942171368,
                        -1,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        0.723922225348,
                        -0.0570857290245,
                        0.552265364726,
                    ],
                    [
                        str("1_2"),
                        0,
                        0.15991259328,
                        0.468401914855,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        -0.0689348112445,
                        0.239549950189,
                        0.57212115096,
                        0,
                        0.749172728898,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        -0.0225155791286,
                        0.285968363141,
                        0.618539556698,
                        -0.162465438799,
                        0.836970467538,
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
                        0.4721804016,
                        -0.2360902008,
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
                        -1,
                    ],
                    [
                        str("0_2"),
                        0.575010099536,
                        0,
                        -0.287505049768,
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
                        -1,
                    ],
                    [
                        str("0_3"),
                        0.685865110811,
                        0,
                        0,
                        -0.342932555406,
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
                        0.5382139455257,
                    ],
                    [
                        str("0_4"),
                        0.425529576141,
                        0,
                        0,
                        0,
                        -0.21276478807,
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
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        0.694826481806,
                        -0.0134419965418,
                        0,
                        -0.333971244361,
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
                        0.541330329235,
                    ],
                    [
                        str("0_2_4"),
                        0.60032414741,
                        0,
                        -0.262190988021,
                        0,
                        -0.037971085684,
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
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        0.666836496324,
                        0,
                        0,
                        -0.361961112674,
                        0.0285428645122,
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
                        0.552265364726,
                    ],
                    [
                        str("0_1_2"),
                        0.628314508135,
                        -0.07995629664,
                        -0.234200957428,
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
                        -1,
                    ],
                    [
                        str("0_1_2_3"),
                        0.742736289904,
                        0.0344674056222,
                        -0.119774975094,
                        -0.28606057548,
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
                        0.749172728898,
                    ],
                    [
                        str("0_1_2_3_4"),
                        0.719526901911,
                        0.0112577895643,
                        -0.14298418157,
                        -0.309269778349,
                        0.0812327193995,
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
                        0.836970467538,
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
                        0.2360902008,
                        -0.2360902008,
                        0,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_2"),
                        0.575010099536,
                        0,
                        -0.575010099536,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_3"),
                        0.342932555406,
                        1.02879766622,
                        0,
                        -1.37173022162,
                        0,
                        0,
                        0.538213945525,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        1.70211830456,
                        0,
                        -1.70211830456,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        0.347413240903,
                        0.988471736541,
                        0,
                        -1.33588497744,
                        0,
                        0,
                        0.541330329235,
                    ],
                    [
                        str("0_2_4"),
                        0.524381976042,
                        0,
                        -0.22061329057,
                        0,
                        -0.303768685472,
                        0,
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        0.361961112674,
                        1.08588333802,
                        -0.228342916098,
                        -1.4478444507,
                        0.228342916098,
                        0,
                        0.552265364726,
                    ],
                    [
                        str("0_1_2"),
                        0.548358211495,
                        -0.07995629664,
                        -0.468401914855,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_2_3"),
                        0.491143120047,
                        0.892649132062,
                        -0.239549950189,
                        -1.14424230192,
                        0,
                        0,
                        0.749172728898,
                    ],
                    [
                        str("0_1_2_3_4"),
                        0.583980351926,
                        0.939067124611,
                        -0.935830118337,
                        -1.2370791134,
                        0.649861755196,
                        0,
                        0.836970467538,
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

    def test_gaussian_style_files_fit_trappe_ua_fit_CT_CT_CT_CT_in_butane_files(
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
            r_squared_atol=0.02,
            opls_force_k0_zero=True,
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
                    [str("1"), 0, 3.89019135095, 0, 0, 0, -1],
                    [str("2"), 0, 0, 3.03580389236, 0, 0, -1],
                    [str("3"), 0, 0, 0, 4.57026969264, 0, 0.755185044925],
                    [str("4"), 0, 0, 0, 0, 3.24432472808, -1],
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
                        0.957438274954,
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
                        -1,
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
                        -1,
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
                        -1,
                    ],
                    [
                        str("0_2"),
                        3.03580389236,
                        0,
                        -3.03580389236,
                        0,
                        0,
                        0,
                        -1,
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
                        -1,
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
                        -1,
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
                        -1,
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
            r_squared_atol=0.02,
            opls_force_k0_zero=False,
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
                    [str("1"), 0, 1.83985239565, 0, 0, 0, -1],
                    [str("2"), 1.44677123926, 0, -0.723385619628, 0, 0, -1],
                    [str("3"), 0, 0, 0, 4.5703138082, 0, 0.755170256901],
                    [str("4"), 0.19559788237, 0, 0, 0, -0.0977989411849, -1],
                    [
                        str("1_3"),
                        0,
                        1.51806767337,
                        0,
                        3.55827093343,
                        0,
                        0.957436115154,
                    ],
                    [
                        str("2_4"),
                        1.44677039049,
                        0,
                        -0.723385195245,
                        0,
                        -0.097795846002,
                        -1,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        4.33332965501,
                        0.35547658682,
                        0.766261090206,
                    ],
                    [
                        str("1_2"),
                        0.191804720892,
                        1.83985769901,
                        -0.723399284839,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        1.81669776403,
                        -0.746584400207,
                        3.85690892695,
                        0,
                        0.998529964618,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        1.84357638993,
                        -0.719706246233,
                        3.88378707851,
                        -0.0940741148625,
                        0.999129174099,
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
                        1.83985239565,
                        -0.919926197825,
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
                        -1,
                    ],
                    [
                        str("0_2"),
                        2.00006677886e-12,
                        0,
                        0.361692809814,
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
                        -1,
                    ],
                    [
                        str("0_3"),
                        4.5703138082,
                        0,
                        0,
                        -2.2851569041,
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
                        0.755170256901,
                    ],
                    [
                        str("0_4"),
                        1.00003338943e-13,
                        0,
                        0,
                        0,
                        0.0488994705924,
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
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        5.0763386068,
                        -0.759033836685,
                        0,
                        -1.77913546672,
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
                        0.957436115154,
                    ],
                    [
                        str("0_2_4"),
                        -0.097795846002,
                        0,
                        0.361692597622,
                        0,
                        0.048897923001,
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
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        4.68880624183,
                        0,
                        0,
                        -2.1666648275,
                        -0.177738293412,
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
                        0.766261090206,
                    ],
                    [
                        str("0_1_2"),
                        1.21236077462,
                        -0.919928849505,
                        0.36169964242,
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
                        -1,
                    ],
                    [
                        str("0_1_2_3"),
                        4.92702229077,
                        -0.908348882015,
                        0.373292200104,
                        -1.92845446348,
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
                        0.998529964618,
                    ],
                    [
                        str("0_1_2_3_4"),
                        4.91358310734,
                        -0.921788194965,
                        0.359853123116,
                        -1.94189353926,
                        0.0470370574312,
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
                        0.999129174099,
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
                        0.919926197825,
                        -0.919926197825,
                        0,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_2"),
                        2.00006677886e-12,
                        0,
                        0.723385619628,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_3"),
                        2.2851569041,
                        6.8554707123,
                        0,
                        -9.1406276164,
                        0,
                        0,
                        0.755170256901,
                    ],
                    [
                        str("0_4"),
                        0.097798941185,
                        0,
                        -0.39119576474,
                        0,
                        0.39119576474,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        2.5381693034,
                        4.57837256346,
                        0,
                        -7.11654186686,
                        0,
                        0,
                        0.957436115154,
                    ],
                    [
                        str("0_2_4"),
                        0,
                        0,
                        0.332201811237,
                        0,
                        0.391183384008,
                        0,
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        2.1666648275,
                        6.49999448252,
                        1.4219063473,
                        -8.66665931002,
                        -1.4219063473,
                        0,
                        0.766261090206,
                    ],
                    [
                        str("0_1_2"),
                        0.292431925112,
                        -0.919928849505,
                        0.723399284839,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_2_3"),
                        2.09021894528,
                        4.87701450841,
                        0.746584400207,
                        -7.7138178539,
                        0,
                        0,
                        0.998529964618,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.14397548799,
                        4.9038924228,
                        0.343409786783,
                        -7.76757415702,
                        0.37629645945,
                        0,
                        0.999129174099,
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

    def test_gaussian_style_files_fit_exp6_fit_CT_CT_CT_CT_in_butane_files(
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
            r_squared_atol=0.02,
            opls_force_k0_zero=False,
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
                    [str("1"), 0, 4.64199086498, 0, 0, 0, 0.0554477339676],
                    [str("2"), 0.812204320046, 0, -0.406102160023, 0, 0, -1],
                    [str("3"), 0, 0, 0, 4.92958479165, 0, 0.477527082073],
                    [str("4"), 0.245247516662, 0, 0, 0, -0.122623758331, -1],
                    [
                        str("1_3"),
                        0,
                        2.44006484197,
                        0,
                        3.30287852792,
                        0,
                        0.984681393995,
                    ],
                    [
                        str("2_4"),
                        0.819462098421,
                        0,
                        -0.40610160909,
                        0,
                        -0.122622010902,
                        -1,
                    ],
                    [
                        str("3_4"),
                        0,
                        0,
                        0,
                        4.41244213668,
                        0.775714796511,
                        0.528782922306,
                    ],
                    [
                        str("1_2"),
                        0,
                        4.00815613582,
                        0.950752493066,
                        0,
                        0,
                        0.132444797386,
                    ],
                    [
                        str("1_2_3"),
                        0,
                        2.61645375375,
                        -0.440977692642,
                        3.47927210449,
                        0,
                        0.998595308163,
                    ],
                    [
                        str("1_2_3_4"),
                        0,
                        2.65144541937,
                        -0.4059866415761,
                        3.51426315104,
                        -0.122469426688,
                        0.999580885743,
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
                        4.64199086498,
                        -2.32099543249,
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
                        0.0554477339676,
                    ],
                    [
                        str("0_2"),
                        0,
                        0,
                        0.203051080012,
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
                        -1,
                    ],
                    [
                        str("0_3"),
                        4.92958479165,
                        0,
                        0,
                        -2.46479239582,
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
                        0.477527082073,
                    ],
                    [
                        str("0_4"),
                        0,
                        0,
                        0,
                        0,
                        0.0613118791655,
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
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        5.74294336989,
                        -1.22003242098,
                        0,
                        -1.65143926396,
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
                        0.984681393995,
                    ],
                    [
                        str("0_2_4"),
                        -0.118992570782,
                        0,
                        0.203050804545,
                        0,
                        0.061311005451,
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
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        5.18815693319,
                        0,
                        0,
                        -2.20622106834,
                        -0.387857398256,
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
                        0.528782922306,
                    ],
                    [
                        str("0_1_2"),
                        4.95890862889,
                        -2.00407806791,
                        -0.475376246533,
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
                        0.132444797386,
                    ],
                    [
                        str("0_1_2_3"),
                        5.6547481656,
                        -1.30822687688,
                        0.220488846321,
                        -1.73963605224,
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
                        0.998595308163,
                    ],
                    [
                        str("0_1_2_3_4"),
                        5.63725250215,
                        -1.32572270968,
                        0.202993320788,
                        -1.75713157552,
                        0.061234713344,
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
                        0.999580885743,
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
                        2.32099543249,
                        -2.32099543249,
                        0,
                        0,
                        0,
                        0,
                        0.0554477339676,
                    ],
                    [
                        str("0_2"),
                        0,
                        0,
                        0.406102160023,
                        0,
                        0,
                        0,
                        -1,
                    ],
                    [
                        str("0_3"),
                        2.46479239582,
                        7.39437718748,
                        0,
                        -9.8591695833,
                        0,
                        0,
                        0.477527082073,
                    ],
                    [
                        str("0_4"),
                        0.122623758331,
                        0,
                        -0.490495033324,
                        0,
                        0.490495033324,
                        0,
                        -1,
                    ],
                    [
                        str("0_1_3"),
                        2.87147168494,
                        3.73428537089,
                        0,
                        -6.60575705584,
                        0,
                        0,
                        0.984681393995,
                    ],
                    [
                        str("0_2_4"),
                        0.0036294401205,
                        0,
                        -0.084386434518,
                        0,
                        0.490488043608,
                        0,
                        -1,
                    ],
                    [
                        str("0_3_4"),
                        2.20622106834,
                        6.61866320502,
                        3.10285918604,
                        -8.82488427336,
                        -3.10285918604,
                        0,
                        0.528782922306,
                    ],
                    [
                        str("0_1_2"),
                        2.95483056098,
                        -2.00407806791,
                        -0.950752493066,
                        0,
                        0,
                        0,
                        0.132444797386,
                    ],
                    [
                        str("0_1_2_3"),
                        2.60688523648,
                        3.91068127986,
                        0.440977692642,
                        -6.95854420898,
                        0,
                        0,
                        0.998595308163,
                    ],
                    [
                        str("0_1_2_3_4"),
                        2.67686764363,
                        3.94567201688,
                        -0.083891065176,
                        -7.02852630208,
                        0.489877706752,
                        0,
                        0.999580885743,
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

    def test_gaussian_log_file_fit_butane_HC_CT_CT_HC_check_all_summed_same_angles_working(
        self,
    ):
        fit_dihedral_with_gomc(
            ["HC", "CT", "CT", "HC"],
            self.get_fn(
                "gaussian_style_output_files/HC_CT_CT_HC_butane/input/starting_coords/HC_CT_CT_HC_butane.mol2"
            ),
            self.get_fn("amber_aa_butane_CT_CT_CT_CT_gmso.xml"),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                self.get_fn(
                    "gaussian_style_output_files/HC_CT_CT_HC_butane/output"
                ): [],
            },
            manual_dihedral_atom_numbers_list=[6, 1, 2, 9],
            qm_engine="gaussian_style_final_files",
            combining_rule=None,
            atom_type_naming_style="general",
            gomc_cpu_cores=1,
            r_squared_min=0.98,
            r_squared_atol=0.02,
            opls_force_k0_zero=False,
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
                    [str(3), 0, 0, 0, 0.272625771312, 0, 0.995091320332],
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
                        0.272625771312,
                        0,
                        0,
                        -0.136312885656,
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
                        0.995091320332,
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
                        0.136312885656,
                        0.408938656968,
                        0,
                        -0.545251542624,
                        0,
                        0,
                        0.995091320332,
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
