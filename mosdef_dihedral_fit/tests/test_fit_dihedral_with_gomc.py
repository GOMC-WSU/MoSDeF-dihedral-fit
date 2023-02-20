import os
import unyt as u
from mosdef_dihedral_fit.tests.base_test import BaseTest
from mosdef_dihedral_fit.utils.io import get_mosdef_dihedral_fit_fn
from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import fit_dihedral_with_gomc
import mosdef_dihedral_fit.utils.math_operations as mdf_math

# user changable variable, as it needs to be run locally
gomc_binary_directory = "/Users/brad/Programs/GOMC/GOMC_2_75/bin"

class TestFitDihedralWithGomc(BaseTest):
    def test_fit_oplsaa_fit_ethane_HC_CT_CT_HC(self):
        fit_dihedral_with_gomc(
            ['HC', 'CT', 'CT', 'HC'],
            get_mosdef_dihedral_fit_fn('ethane_aa.mol2'),
            get_mosdef_dihedral_fit_fn('oplsaa_ethane_HC_CT_CT_HC.xml'),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                get_mosdef_dihedral_fit_fn('HC_CT_CT_HC_multiplicity_1.log'): [0],
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
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('gomc_raw_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt') is True
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
                        0.0,
                        0.0,
                        0.0,
                        0.3140037484218711,
                        0.0,
                        0.9987911455740498
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
                        0.3140037484278298,
                        -0.0,
                        -0.0,
                        -0.1570018742139149,
                        -0.0,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        0.9987911455740498
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
                        str('0_3'),
                        0.1570018742126558,
                        0.4710056226379674,
                        0.0,
                        -0.6280074968506232,
                        -0.0,
                        0.0,
                        0.9987911455740498
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


    def test_fit_oplsaa_protonated_fragment_CT_CT_C_OH_in_COOH(self):
        fit_dihedral_with_gomc(
            ['CT', 'CT', 'C', 'OH'],
            get_mosdef_dihedral_fit_fn('protonated_fragment_CT_CT_C_OH_in_COOH.mol2'),
            get_mosdef_dihedral_fit_fn('oplsaa_CT_CT_C_OH_in_COOH.xml'),
            298.15 * u.Kelvin,
            gomc_binary_directory,
            {
                get_mosdef_dihedral_fit_fn('CT_CT_C_OH_multiplicity_1.log'): [0],
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
        assert os.path.isfile('all_normalized_energies_OPLS_fit_3_in_kcal_per_mol.txt') is True
        assert os.path.isfile('gomc_raw_energies_in_Kelvin.txt') is True
        assert os.path.isfile('gomc_raw_OPLS_fit_3_energies_in_Kelvin.txt') is True
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
                        0.0,
                        2.0040856887989777,
                        0.0,
                        0.0,
                        0.0,
                        -0.37186627820142193
                    ],
                    [
                        str('2'),
                        0.0,
                        0.0,
                        2.13847194795373,
                        0.0,
                        0.0,
                        0.26478715881215786
                    ],
                    [
                        str('3'),
                        0.0,
                        0.0,
                        0.0,
                        1.79252468631398,
                        0.0,
                        -1.29046011346467
                    ],
                    [
                        str('4'),
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        1.7121917318224398,
                        -1.612454729569004
                    ],
                    [
                        str('1_3'),
                        0.0,
                        1.4563236530479318,
                        0.0,
                        0.8216420999814841,
                        0.0,
                        0.05706377385372252
                    ],
                    [
                        str('2_4'),
                        0.0,
                        0.0,
                        1.7946173953391606,
                        0.0,
                        0.5157802309204645,
                        0.43381145515687736
                    ],
                    [
                        str('3_4'),
                        0.0,
                        0.0,
                        0.0,
                        1.1719159996737067,
                        0.9309130432196181,
                        -0.739857285932813
                    ],
                    [
                        str('1_2'),
                        0.0,
                        1.0411943336512077,
                        1.4443421195517043,
                        0.0,
                        0.0,
                        0.953575526750655
                    ],
                    [
                        str('1_2_3'),
                        0.0,
                        0.9250475946051446,
                        1.3281962133659146,
                        0.2903652648033817,
                        0.0,
                        0.9985732062427662
                    ],
                    [
                        str('1_2_3_4'),
                        0.0,
                        0.9140793554355323,
                        1.3172277155021408,
                        0.2793968919383311,
                        0.03838925870074258,
                        0.9992955355787176
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
                            print(f'split_line_i               = {split_line_i}')
                            print( f'correct_line_values[i] = {correct_line_values[i]}')

                            print(f'mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i)    = {mdf_math.round_to_sig_figs(float(split_line_i[j]), number_sig_i)}')
                            print(f'mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i) = {mdf_math.round_to_sig_figs(correct_line_values[i][j], number_sig_i)}')


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
                        2.0040856887989777,
                        -1.0020428443994889,
                        -0.0,
                        -0.0,
                        -0.0,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        -0.3718662782014219,
                    ],
                    [
                        str('0_2'),
                        2.13847194795373,
                        -0.0,
                        -1.069235973976865,
                        -0.0,
                        -0.0,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        0.2647871588121578,
                    ],
                    [
                        str('0_3'),
                        1.79252468631398,
                        -0.0,
                        -0.0,
                        -0.89626234315699,
                        -0.0,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        -1.29046011346467,
                    ],
                    [
                        str('0_4'),
                        1.712191731822439,
                        -0.0,
                        -0.0,
                        -0.0,
                        -0.8560958659112199,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        -1.612454729569004,
                    ],
                    [
                        str('0_1_3'),
                        2.2779657530294157,
                        -0.7281618265239659,
                        -0.0,
                        -0.41082104999074204,
                        -0.0,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        0.0570637738537225,
                    ],
                    [
                        str('0_2_4'),
                        2.310397626259625,
                        -0.0,
                        -0.8973086976695803,
                        -0.0,
                        -0.25789011546023227,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        0.4338114551568773,
                    ],
                    [
                        str('0_3_4'),
                        2.1028290428933247,
                        -0.0,
                        -0.0,
                        -0.5859579998368534,
                        -0.465456521609809,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        -0.739857285932813,
                    ],
                    [
                        str('0_1_2'),
                        2.485536453202912,
                        -0.5205971668256038,
                        -0.7221710597758522,
                        -0.0,
                        -0.0,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        0.953575526750655,
                    ],
                    [
                        str('0_1_2_3'),
                        2.543609072774441,
                        -0.4625237973025723,
                        -0.6640981066829573,
                        -0.14518263240169085,
                        -0.0,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        0.9985732062427662,
                    ],

                    [
                        str('0_1_2_3_4'),
                        2.549093221576747,
                        -0.4570396777177662,
                        -0.6586138577510704,
                        -0.13969844596916556,
                        -0.01919462935037125,
                        0.0,
                        0,
                        1,
                        2,
                        3,
                        4,
                        5,
                        90.0,
                        180.0,
                        0.0,
                        180.0,
                        0,
                        180.0,
                        0.9992955355787176,
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
                        1.0020428443994889,
                         -1.0020428443994889,
                        0.0,
                        -0.0,
                        -0.0,
                        0.0,
                        -0.3718662782014219
                    ],
                    [
                        str('0_2'),
                        2.13847194795373,
                        0.0,
                        -2.13847194795373,
                        -0.0,
                        -0.0,
                        0.0,
                        0.2647871588121578
                    ],
                    [
                        str('0_3'),
                        0.89626234315699,
                        2.68878702947097,
                        0.0,
                        -3.58504937262796,
                        -0.0,
                        0.0,
                        -1.29046011346467
                    ],
                    [
                        str('0_4'),
                        0.0,
                        0.0,
                        6.848766927289759,
                        -0.0,
                        -6.848766927289759,
                        0.0,
                        -1.612454729569004
                    ],
                    [
                        str('0_1_3'),
                        1.1389828765147079,
                        0.5043013234482602,
                        0.0,
                        -1.6432841999629682,
                        -0.0,
                        0.0,
                        0.0570637738537225
                    ],
                    [
                        str('0_2_4'),
                        1.7946173953391606,
                        0.0,
                        0.26850352834269753,
                        -0.0,
                        -2.063120923681858,
                        0.0,
                        0.4338114551568773
                    ],
                    [
                        str('0_3_4'),
                        0.5859579998368534,
                        1.7578739995105601,
                        3.723652172878472,
                        -2.3438319993474135,
                        -3.723652172878472,
                        0.0,
                        -0.739857285932813
                    ],
                    [
                        str('0_1_2'),
                        1.964939286377308,
                        -0.5205971668256038,
                        -1.4443421195517043,
                        -0.0,
                        -0.0,
                        0.0,
                        0.953575526750655
                    ],
                    [
                        str('0_1_2_3'),
                        1.9359026430701778,
                        -0.026975900097499728,
                        -1.3281962133659146,
                        -0.5807305296067634,
                        -0.0,
                        0.0,
                        0.9985732062427662
                    ],
                    [
                        str('0_1_2_3_4'),
                        1.9139658391890726,
                        -0.03794433981026957,
                        -1.1636706806991708,
                        -0.5587937838766622,
                        -0.15355703480297,
                        0.0,
                        0.9992955355787176
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
