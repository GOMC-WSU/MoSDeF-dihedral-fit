import os
import unyt as u
from mosdef_dihedral_fit.tests.base_test import BaseTest
from mosdef_dihedral_fit.utils.io import get_mosdef_dihedral_fit_fn
from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import fit_dihedral_with_gomc
import mosdef_dihedral_fit.utils.math_operations as mdf_math

# user changable variable, as it needs to be run locally
gomc_binary_directory = "/Users/brad/Programs/GOMC/GOMC_2_75/bin"

class TestFitDihedralWithGomc(BaseTest):
    def test_gaussian_log_file_fit_oplsaa_fit_ethane_HC_CT_CT_HC(self):
        fit_dihedral_with_gomc(
            ['HC', 'CT', 'CT', 'HC'],
            get_mosdef_dihedral_fit_fn('gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2'),
            get_mosdef_dihedral_fit_fn('oplsaa_ethane_HC_CT_CT_HC.xml'),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                get_mosdef_dihedral_fit_fn('gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log'): [0],
            },
            zeroed_dihedral_atom_types=None,
            qm_engine="gaussian",
            override_VDWGeometricSigma=True,
            atom_type_naming_style='general',
            gomc_cpu_cores=1,
            fit_min_validated_r_squared=0.99,
            fit_validation_r_squared_rtol=1e-03
        )

        assert os.path.isfile('all_normalized_energies_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt') is False
        assert os.path.isfile('gomc_raw_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt') is False
        assert os.path.isfile('opls_all_single_fit_dihedral_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_all_summed_dihedrals_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('periodic_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('RB_torsion_k_constants_fit_energy.txt') is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str(3),
                        0,
                        0,
                        0,
                        0.314003748427,
                        0,
                        0.998791145574
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'n0_kcal_per_mol',
                        'n1_kcal_per_mol',
                        'n2_kcal_per_mol',
                        'n3_kcal_per_mol',
                        'n4_kcal_per_mol',
                        'n5_kcal_per_mol',
                        'd0_kcal_per_mol',
                        'd1_kcal_per_mol',
                        'd2_kcal_per_mol',
                        'd3_kcal_per_mol',
                        'd4_kcal_per_mol',
                        'd5_kcal_per_mol',
                        'r_squared',
                    ],
                    [
                        str('0_3'),
                        0.314003748427,
                        0,
                        0,
                        -0.157001874214,
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
                        0.998791145574
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

        # check the periodic dihedral file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'r_squared',
                    ],
                    [
                        str('0_3'),
                        0.157001874212,
                        0.471005622638,
                        0,
                        -0.62800749685,
                        0,
                        0,
                        0.998791145574
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

    def test_gaussian_log_file_fit_oplsaa_fit_ethane_HC_CT_CT_HC_with_2_log_files(self):
        fit_dihedral_with_gomc(
            ['HC', 'CT', 'CT', 'HC'],
            get_mosdef_dihedral_fit_fn('gaussian/HC_CT_CT_HC/input/starting_coords/ethane_aa.mol2'),
            get_mosdef_dihedral_fit_fn('oplsaa_ethane_HC_CT_CT_HC.xml'),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                get_mosdef_dihedral_fit_fn(
                    'gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1.log'): [0],
                get_mosdef_dihedral_fit_fn(
                    'gaussian/HC_CT_CT_HC/output/HC_CT_CT_HC_multiplicity_1_copy_for_test.log'): [0],
            },
            zeroed_dihedral_atom_types=None,
            qm_engine="gaussian",
            override_VDWGeometricSigma=True,
            atom_type_naming_style='general',
            gomc_cpu_cores=1,
            fit_min_validated_r_squared=0.99,
            fit_validation_r_squared_rtol=1e-03
        )

        assert os.path.isfile('all_normalized_energies_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt') is False
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt') is False
        assert os.path.isfile('gomc_raw_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt') is False
        assert os.path.isfile('gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt') is False
        assert os.path.isfile('opls_all_single_fit_dihedral_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_all_summed_dihedrals_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('periodic_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('RB_torsion_k_constants_fit_energy.txt') is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str(3),
                        0,
                        0,
                        0,
                        0.314003748425,
                        0,
                        0.998791145574
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'n0_kcal_per_mol',
                        'n1_kcal_per_mol',
                        'n2_kcal_per_mol',
                        'n3_kcal_per_mol',
                        'n4_kcal_per_mol',
                        'n5_kcal_per_mol',
                        'd0_kcal_per_mol',
                        'd1_kcal_per_mol',
                        'd2_kcal_per_mol',
                        'd3_kcal_per_mol',
                        'd4_kcal_per_mol',
                        'd5_kcal_per_mol',
                        'r_squared',
                    ],
                    [
                        str('0_3'),
                        0.314003748425,
                        0,
                        0,
                        -0.157001874212,
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
                        0.998791145574
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

        # check the periodic dihedral file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'r_squared',
                    ],
                    [
                        str('0_3'),
                        0.157001874212,
                        0.471005622638,
                        0,
                        -0.62800749685,
                        0,
                        0,
                        0.998791145574
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)


    def test_gaussian_log_file_fit_oplsaa_protonated_fragment_CT_CT_C_OH_in_COOH(self):
        fit_dihedral_with_gomc(
            ['CT', 'CT', 'C', 'OH'],
            get_mosdef_dihedral_fit_fn(
                'gaussian/CT_CT_C_OH/input/starting_coords/protonated_fragment_CT_CT_C_OH_in_COOH.mol2'),
            get_mosdef_dihedral_fit_fn('oplsaa_CT_CT_C_OH_in_COOH.xml'),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                get_mosdef_dihedral_fit_fn(
                    'gaussian/CT_CT_C_OH/output/CT_CT_C_OH_multiplicity_1.log'): [0],
            },
            zeroed_dihedral_atom_types=[['CT', 'CT', 'C', 'O_3']],
            qm_engine="gaussian",
            override_VDWGeometricSigma=True,
            atom_type_naming_style='general',
            gomc_cpu_cores=1,
            fit_min_validated_r_squared=0.99,
            fit_validation_r_squared_rtol=1e-03
        )

        assert os.path.isfile('all_normalized_energies_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('gomc_raw_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('opls_all_single_fit_dihedral_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_all_summed_dihedrals_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('periodic_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('RB_torsion_k_constants_fit_energy.txt') is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str('1'),
                        0,
                        2.0040856888,
                        0,
                        0,
                        0,
                        -0.371866278201
                    ],
                    [
                        str('2'),
                        0,
                        0,
                        2.13847194795,
                        0,
                        0,
                        0.264787158812
                    ],
                    [
                        str('3'),
                        0,
                        0,
                        0,
                        1.79252468631,
                        0,
                        -1.29046011346
                    ],
                    [
                        str('4'),
                        0,
                        0,
                        0,
                        0,
                        1.71219173182,
                        -1.61245472957
                    ],
                    [
                        str('1_3'),
                        0,
                        1.45632365305,
                        0,
                        0.821642099981,
                        0,
                        0.0570637738537
                    ],
                    [
                        str('2_4'),
                        0,
                        0,
                        1.79461739534,
                        0,
                        0.51578023092,
                        0.433811455157
                    ],
                    [
                        str('3_4'),
                        0,
                        0,
                        0,
                        1.17191599967,
                        0.93091304322,
                        -0.739857285933
                    ],
                    [
                        str('1_2'),
                        0,
                        1.04119433365,
                        1.44434211955,
                        0,
                        0,
                        0.953575526751
                    ],
                    [
                        str('1_2_3'),
                        0,
                        0.925047594605,
                        1.32819621337,
                        0.290365264803,
                        0,
                        0.998573206243
                    ],
                    [
                        str('1_2_3_4'),
                        0,
                        0.914079355436,
                        1.3172277155,
                        0.279396891938,
                        0.0383892587007,
                        0.999295535579
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:

                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'n0_kcal_per_mol',
                        'n1_kcal_per_mol',
                        'n2_kcal_per_mol',
                        'n3_kcal_per_mol',
                        'n4_kcal_per_mol',
                        'n5_kcal_per_mol',
                        'd0_kcal_per_mol',
                        'd1_kcal_per_mol',
                        'd2_kcal_per_mol',
                        'd3_kcal_per_mol',
                        'd4_kcal_per_mol',
                        'd5_kcal_per_mol',
                        'r_squared',
                    ],
                    [
                        str('0_1'),
                        2.0040856888,
                        -1.0020428444,
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
                        -0.371866278201,
                    ],
                    [
                        str('0_2'),
                        2.13847194795,
                        0,
                        -1.06923597397,
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
                        str('0_3'),
                        1.79252468631,
                        0,
                        0,
                        -0.896262343155,
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
                        str('0_4'),
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
                        str('0_1_3'),
                        2.27796575303,
                        -0.728161826525,
                        0,
                        -0.41082104999,
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
                        str('0_2_4'),
                        2.31039762626,
                        0,
                        -0.8973086976,
                        0,
                        -0.25789011546,
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
                        str('0_3_4'),
                        2.10282904289,
                        0,
                        0,
                        -0.585957999835,
                        -0.46545652161,
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
                        str('0_1_2'),
                        2.4855364532,
                        -0.520597166825,
                        -0.722171059775,
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
                        str('0_1_2_3'),
                        2.54360907278,
                        -0.462523797302,
                        -0.664098106685,
                        -0.145182632402,
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
                        str('0_1_2_3_4'),
                        2.54909322157,
                        -0.457039677718,
                        -0.65861385775,
                        -0.139698445969,
                        -0.0191946293504,
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
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)


        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str('0_1'),
                        1.0020428444,
                        -1.0020428444,
                        0,
                        0,
                        0,
                        0,
                        -0.371866278201
                    ],
                    [
                        str('0_2'),
                        2.13847194795,
                        0,
                        -2.13847194795,
                        0,
                        0,
                        0,
                        0.264787158812
                    ],
                    [
                        str('0_3'),
                        0.896262343155,
                        2.68878702946,
                        0,
                        -3.58504937262,
                        0,
                        0,
                        -1.29046011346
                    ],
                    [
                        str('0_4'),
                        0,
                        0,
                        6.84876692728,
                        0,
                        -6.84876692728,
                        0,
                        -1.61245472957
                    ],
                    [
                        str('0_1_3'),
                        1.13898287652,
                        0.504301323446,
                        0,
                        -1.64328419996,
                        0,
                        0,
                        0.0570637738537
                    ],
                    [
                        str('0_2_4'),
                        1.79461739534,
                        0,
                        0.26850352834,
                        0,
                        -2.06312092368,
                        0,
                        0.433811455157
                    ],
                    [
                        str('0_3_4'),
                        0.585957999835,
                        1.7578739995,
                        3.72365217288,
                        -2.34383199934,
                        -3.72365217288,
                        0,
                        -0.739857285933
                    ],
                    [
                        str('0_1_2'),
                        1.96493928638,
                        -0.520597166825,
                        -1.44434211955,
                        0,
                        0,
                        0,
                        0.953575526751
                    ],
                    [
                        str('0_1_2_3'),
                        1.93590264307,
                        -0.026975900098,
                        -1.32819621337,
                        -0.580730529606,
                        0,
                        0,
                        0.998573206243
                    ],
                    [
                        str('0_1_2_3_4'),
                        1.91396583919,
                        -0.037944339811,
                        -1.1636706807,
                        -0.558793783876,
                        -0.153557034803,
                        0,
                        0.999295535579
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)


    def test_gaussian_style_files_fit_oplsaa_fit_CT_CT_C_OH_in_COOH_missing_1st_point(self):
        fit_dihedral_with_gomc(
            ['CT', 'CT', 'C', 'OH'],
            get_mosdef_dihedral_fit_fn(
                'gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2'),
            get_mosdef_dihedral_fit_fn('oplsaa_CT_CT_C_OH_in_COOH_zeroed.xml'),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                get_mosdef_dihedral_fit_fn('gaussian_style_output_files/CT_CT_C_OH/output'): [0],
            },
            manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
            zeroed_dihedral_atom_types=None,
            qm_engine="gaussian_style_final_files",
            override_VDWGeometricSigma=True,
            atom_type_naming_style='general',
            gomc_cpu_cores=1,
            fit_min_validated_r_squared=0.99,
            fit_validation_r_squared_rtol=0.02
        )

        assert os.path.isfile('all_normalized_energies_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('gomc_raw_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('opls_all_single_fit_dihedral_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_all_summed_dihedrals_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('periodic_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('RB_torsion_k_constants_fit_energy.txt') is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str('1'),
                        0,
                        1.80634795241,
                        0,
                        0,
                        0,
                        -0.268742410571
                    ],
                    [
                        str('2'),
                        0,
                        0,
                        1.94879786206,
                        0,
                        0,
                        0.448421126946
                    ],
                    [
                        str('3'),
                        0,
                        0,
                        0,
                        1.58083737128,
                        0,
                        -1.50863674602
                    ],
                    [
                        str('4'),
                        0,
                        0,
                        0,
                        0,
                        1.5121963942,
                        -1.54542980028
                    ],
                    [
                        str('1_3'),
                        0,
                        1.41759854611,
                        0,
                        0.618102768602,
                        0,
                        -0.00379371805915
                    ],
                    [
                        str('2_4'),
                        0,
                        0,
                        1.69725316221,
                        0,
                        0.37780223364,
                        0.551946157778
                    ],
                    [
                        str('3_4'),
                        0,
                        0,
                        0,
                        0.966769776613,
                        0.903643214082,
                        -0.897967473282
                    ],
                    [
                        str('1_2'),
                        0,
                        0.882099545744,
                        1.34468116365,
                        0,
                        0,
                        0.986743110285
                    ],
                    [
                        str('1_2_3'),
                        0,
                        0.83099823616,
                        1.29385707544,
                        0.136793004885,
                        0,
                        0.997926382517
                    ],
                    [
                        str('1_2_3_4'),
                        0,
                        0.836836322925,
                        1.29869945716,
                        0.141702758859,
                        -0.0179044297589,
                        0.998103759798
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:

                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'n0_kcal_per_mol',
                        'n1_kcal_per_mol',
                        'n2_kcal_per_mol',
                        'n3_kcal_per_mol',
                        'n4_kcal_per_mol',
                        'n5_kcal_per_mol',
                        'd0_kcal_per_mol',
                        'd1_kcal_per_mol',
                        'd2_kcal_per_mol',
                        'd3_kcal_per_mol',
                        'd4_kcal_per_mol',
                        'd5_kcal_per_mol',
                        'r_squared',
                    ],
                    [
                        str('0_1'),
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
                        str('0_2'),
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
                        str('0_3'),
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
                        str('0_4'),
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
                        str('0_1_3'),
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
                        str('0_2_4'),
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
                        str('0_3_4'),
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
                        str('0_1_2'),
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
                        str('0_1_2_3'),
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
                        str('0_1_2_3_4'),
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
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)


        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str('0_1'),
                        0.903173976205,
                        -0.903173976205,
                        0,
                        0,
                        0,
                        0,
                        -0.268742410571
                    ],
                    [
                        str('0_2'),
                        1.94879786206,
                        0,
                        -1.94879786206,
                        0,
                        0,
                        0,
                        0.448421126946
                    ],
                    [
                        str('0_3'),
                        0.79041868564,
                        2.37125605692,
                        0,
                        -3.16167474256,
                        0,
                        0,
                        -1.50863674602
                    ],
                    [
                        str('0_4'),
                        0,
                        0,
                        6.0487855768,
                        0,
                        -6.0487855768,
                        0,
                        -1.54542980028
                    ],
                    [
                        str('0_1_3'),
                        1.01785065736,
                        0.218354879848,
                        0,
                        -1.2362055372,
                        0,
                        0,
                        -0.00379371805915
                    ],
                    [
                        str('0_2_4'),
                        1.69725316221,
                        0,
                        -0.18604422765,
                        0,
                        -1.51120893456,
                        0,
                        0.551946157778
                    ],
                    [
                        str('0_3_4'),
                        0.483384888306,
                        1.45015466492,
                        3.61457285633,
                        -1.93353955323,
                        -3.61457285633,
                        0,
                        -0.897967473282
                    ],
                    [
                        str('0_1_2'),
                        1.78573093652,
                        -0.441049772872,
                        -1.34468116365,
                        0,
                        0,
                        0,
                        0.986743110285
                    ],
                    [
                        str('0_1_2_3'),
                        1.77775269597,
                        -0.210309610756,
                        -1.29385707544,
                        -0.27358600977,
                        0,
                        0,
                        0.997926382517
                    ],
                    [
                        str('0_1_2_3_4'),
                        1.78796899805,
                        -0.205864023174,
                        -1.3703171762,
                        -0.283405517718,
                        0.0716177190356,
                        0,
                        0.998103759798
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)


    def test_gaussian_style_files_fit_oplsaa_fit_CT_CT_C_OH_in_COOH(self):
        fit_dihedral_with_gomc(
            ['CT', 'CT', 'C', 'OH'],
            get_mosdef_dihedral_fit_fn(
                'gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2'),
            get_mosdef_dihedral_fit_fn('oplsaa_CT_CT_C_OH_in_COOH_zeroed.xml'),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                get_mosdef_dihedral_fit_fn('gaussian_style_output_files/CT_CT_C_OH/output'): [],
            },
            manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
            zeroed_dihedral_atom_types=None,
            qm_engine="gaussian_style_final_files",
            override_VDWGeometricSigma=True,
            atom_type_naming_style='general',
            gomc_cpu_cores=1,
            fit_min_validated_r_squared=0.99,
            fit_validation_r_squared_rtol=5e-03
        )

        assert os.path.isfile('all_normalized_energies_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('gomc_raw_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('opls_all_single_fit_dihedral_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_all_summed_dihedrals_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('periodic_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('RB_torsion_k_constants_fit_energy.txt') is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str('1'),
                        0,
                        1.87103578693,
                        0,
                        0,
                        0,
                        -0.184214001408
                    ],
                    [
                        str('2'),
                        0,
                        0,
                        2.01338731788,
                        0,
                        0,
                        0.468487865197
                    ],
                    [
                        str('3'),
                        0,
                        0,
                        0,
                        1.64674772188,
                        0,
                        -1.32197179229
                    ],
                    [
                        str('4'),
                        0,
                        0,
                        0,
                        0,
                        1.5768734149,
                        -1.35119291697
                    ],
                    [
                        str('1_3'),
                        0,
                        1.45815624076,
                        0,
                        0.656469146108,
                        0,
                        0.0786991845116
                    ],
                    [
                        str('2_4'),
                        0,
                        0,
                        1.73604078977,
                        0,
                        0.416554745122,
                        0.579202191408
                    ],
                    [
                        str('3_4'),
                        0,
                        0,
                        0,
                        1.00514371097,
                        0.944164973554,
                        -0.735494529166
                    ],
                    [
                        str('1_2'),
                        0,
                        0.920441371308,
                        1.38301174904,
                        0,
                        0,
                        0.984123002941
                    ],
                    [
                        str('1_2_3'),
                        0,
                        0.858889958603,
                        1.32179425273,
                        0.164766862242,
                        0,
                        0.998396296125
                    ],
                    [
                        str('1_2_3_4'),
                        0,
                        0.857564707345,
                        1.32069502687,
                        0.163652345714,
                        0.00406432148459,
                        0.998404336839
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:

                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'n0_kcal_per_mol',
                        'n1_kcal_per_mol',
                        'n2_kcal_per_mol',
                        'n3_kcal_per_mol',
                        'n4_kcal_per_mol',
                        'n5_kcal_per_mol',
                        'd0_kcal_per_mol',
                        'd1_kcal_per_mol',
                        'd2_kcal_per_mol',
                        'd3_kcal_per_mol',
                        'd4_kcal_per_mol',
                        'd5_kcal_per_mol',
                        'r_squared',
                    ],
                    [
                        str('0_1'),
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
                        str('0_2'),
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
                        str('0_3'),
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
                        str('0_4'),
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
                        str('0_1_3'),
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
                        str('0_2_4'),
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
                        str('0_3_4'),
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
                        str('0_1_2'),
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
                        str('0_1_2_3'),
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
                        str('0_1_2_3_4'),
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
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)


        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str('0_1'),
                        0.935517893465,
                        -0.935517893465,
                        0,
                        0,
                        0,
                        0,
                        -0.184214001408
                    ],
                    [
                        str('0_2'),
                        2.01338731788,
                        0,
                        -2.01338731788,
                        0,
                        0,
                        0,
                        0.468487865197
                    ],
                    [
                        str('0_3'),
                        0.82337386094,
                        2.47012158282,
                        0,
                        -3.29349544376,
                        0,
                        0,
                        -1.32197179229
                    ],
                    [
                        str('0_4'),
                        0,
                        0,
                        6.3074936596,
                        0,
                        -6.3074936596,
                        0,
                        -1.35119291697
                    ],
                    [
                        str('0_1_3'),
                        1.05731269343,
                        0.255625598782,
                        0,
                        -1.31293829222,
                        0,
                        0,
                        0.0786991845116
                    ],
                    [
                        str('0_2_4'),
                        1.73604078977,
                        0,
                        -0.069821809282,
                        0,
                        -1.66621898049,
                        0,
                        0.579202191408
                    ],
                    [
                        str('0_3_4'),
                        0.502571855485,
                        1.50771556646,
                        3.77665989422,
                        -2.01028742194,
                        -3.77665989422,
                        0,
                        -0.735494529166
                    ],
                    [
                        str('0_1_2'),
                        1.84323243469,
                        -0.460220685654,
                        -1.38301174904,
                        0,
                        0,
                        0,
                        0.984123002941
                    ],
                    [
                        str('0_1_2_3'),
                        1.83362266315,
                        -0.182294685938,
                        -1.32179425273,
                        -0.329533724484,
                        0,
                        0,
                        0.998396296125
                    ],
                    [
                        str('0_1_2_3_4'),
                        1.8313035534 ,
                        -0.183303835102,
                        -1.30443774093,
                        -0.327304691428,
                        -0.0162572859384,
                        0,
                        0.998404336839
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)


    def test_gaussian_style_files_fit_oplsaa_fit_CT_CT_C_OH_in_COOH_2_files_missing_1_first_point(self):
        fit_dihedral_with_gomc(
            ['CT', 'CT', 'C', 'OH'],
            get_mosdef_dihedral_fit_fn(
                'gaussian_style_output_files/CT_CT_C_OH/input/starting_coords/CT_CT_C_3_OH.mol2'),
            get_mosdef_dihedral_fit_fn('oplsaa_CT_CT_C_OH_in_COOH.xml'),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                get_mosdef_dihedral_fit_fn('gaussian_style_output_files/CT_CT_C_OH_split_part_2/output'): [],
                get_mosdef_dihedral_fit_fn('gaussian_style_output_files/CT_CT_C_OH_split_part_1/output'): [0],
            },
            manual_dihedral_atom_numbers_list=[3, 2, 1, 4],
            zeroed_dihedral_atom_types=[['CT', 'CT', 'C', 'O_3']],
            qm_engine="gaussian_style_final_files",
            override_VDWGeometricSigma=True,
            atom_type_naming_style='general',
            gomc_cpu_cores=1,
            fit_min_validated_r_squared=0.99,
            fit_validation_r_squared_rtol=0.02
        )

        assert os.path.isfile('all_normalized_energies_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_2_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_1_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_2_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('all_normalized_energies_OPLS_fit_4_in_kcal_per_mol.txt') is True
        assert os.path.isfile('gomc_raw_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_2_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_1_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_2_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_2_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_4_energies_in_Kelvin.txt') is True
        assert os.path.isfile('opls_all_single_fit_dihedral_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_all_summed_dihedrals_k_constants_figure.pdf') is True
        assert os.path.isfile('opls_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('periodic_dihedral_k_constants_fit_energy.txt') is True
        assert os.path.isfile('RB_torsion_k_constants_fit_energy.txt') is True

        # check the OPLS dihedral file
        with open("opls_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str('1'),
                        0,
                        1.80634795241,
                        0,
                        0,
                        0,
                        -0.268742410571
                    ],
                    [
                        str('2'),
                        0,
                        0,
                        1.94879786206,
                        0,
                        0,
                        0.448421126946
                    ],
                    [
                        str('3'),
                        0,
                        0,
                        0,
                        1.58083737128,
                        0,
                        -1.50863674602
                    ],
                    [
                        str('4'),
                        0,
                        0,
                        0,
                        0,
                        1.5121963942,
                        -1.54542980028
                    ],
                    [
                        str('1_3'),
                        0,
                        1.41759854611,
                        0,
                        0.618102768602,
                        0,
                        -0.00379371805915
                    ],
                    [
                        str('2_4'),
                        0,
                        0,
                        1.69725316221,
                        0,
                        0.37780223364,
                        0.551946157778
                    ],
                    [
                        str('3_4'),
                        0,
                        0,
                        0,
                        0.966769776613,
                        0.903643214082,
                        -0.897967473282
                    ],
                    [
                        str('1_2'),
                        0,
                        0.882099545744,
                        1.34468116365,
                        0,
                        0,
                        0.986743110285
                    ],
                    [
                        str('1_2_3'),
                        0,
                        0.83099823616,
                        1.29385707544,
                        0.136793004885,
                        0,
                        0.997926382517
                    ],
                    [
                        str('1_2_3_4'),
                        0,
                        0.836836322925,
                        1.29869945716,
                        0.141702758859,
                        -0.0179044297589,
                        0.998103759798
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:

                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)

        # check the periodic dihedral file
        with open("periodic_dihedral_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'n0_kcal_per_mol',
                        'n1_kcal_per_mol',
                        'n2_kcal_per_mol',
                        'n3_kcal_per_mol',
                        'n4_kcal_per_mol',
                        'n5_kcal_per_mol',
                        'd0_kcal_per_mol',
                        'd1_kcal_per_mol',
                        'd2_kcal_per_mol',
                        'd3_kcal_per_mol',
                        'd4_kcal_per_mol',
                        'd5_kcal_per_mol',
                        'r_squared',
                    ],
                    [
                        str('0_1'),
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
                        str('0_2'),
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
                        str('0_3'),
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
                        str('0_4'),
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
                        str('0_1_3'),
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
                        str('0_2_4'),
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
                        str('0_3_4'),
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
                        str('0_1_2'),
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
                        str('0_1_2_3'),
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
                        str('0_1_2_3_4'),
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
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)


        # check the RB torsion file
        with open("RB_torsion_k_constants_fit_energy.txt", "r") as fp:
            number_sig_i = 4
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                split_line_i = line.split()
                correct_line_values = [
                    [
                        'non_zero_k_constants',
                        'k0_kcal_per_mol',
                        'k1_kcal_per_mol',
                        'k2_kcal_per_mol',
                        'k3_kcal_per_mol',
                        'k4_kcal_per_mol',
                        'k5_kcal_per_mol',
                        'r_squared'
                    ],
                    [
                        str('0_1'),
                        0.903173976205,
                        -0.903173976205,
                        0,
                        0,
                        0,
                        0,
                        -0.268742410571
                    ],
                    [
                        str('0_2'),
                        1.94879786206,
                        0,
                        -1.94879786206,
                        0,
                        0,
                        0,
                        0.448421126946
                    ],
                    [
                        str('0_3'),
                        0.79041868564,
                        2.37125605692,
                        0,
                        -3.16167474256,
                        0,
                        0,
                        -1.50863674602
                    ],
                    [
                        str('0_4'),
                        0,
                        0,
                        6.0487855768,
                        0,
                        -6.0487855768,
                        0,
                        -1.54542980028
                    ],
                    [
                        str('0_1_3'),
                        1.01785065736,
                        0.218354879848,
                        0,
                        -1.2362055372,
                        0,
                        0,
                        -0.00379371805915
                    ],
                    [
                        str('0_2_4'),
                        1.69725316221,
                        0,
                        -0.18604422765,
                        0,
                        -1.51120893456,
                        0,
                        0.551946157778
                    ],
                    [
                        str('0_3_4'),
                        0.483384888306,
                        1.45015466492,
                        3.61457285633,
                        -1.93353955323,
                        -3.61457285633,
                        0,
                        -0.897967473282
                    ],
                    [
                        str('0_1_2'),
                        1.78573093652,
                        -0.441049772872,
                        -1.34468116365,
                        0,
                        0,
                        0,
                        0.986743110285
                    ],
                    [
                        str('0_1_2_3'),
                        1.77775269597,
                        -0.210309610756,
                        -1.29385707544,
                        -0.27358600977,
                        0,
                        0,
                        0.997926382517
                    ],
                    [
                        str('0_1_2_3_4'),
                        1.78796899805,
                        -0.205864023174,
                        -1.3703171762,
                        -0.283405517718,
                        0.0716177190356,
                        0,
                        0.998103759798
                    ],

                ]

                if i == 0:
                    assert split_line_i == correct_line_values[i]

                else:
                    assert len(split_line_i) == len(correct_line_values[i])
                    for j in range(0, len(correct_line_values[i])):
                        # check the string listing what cosin powers are used
                        if j == 0:
                            assert str(split_line_i[j]) == str(correct_line_values[i][j])

                        # check the k-values and the r-squared fit
                        else:
                            assert mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i) == \
                                   mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)