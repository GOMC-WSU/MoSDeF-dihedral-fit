.. _secondary_data_structures:

=========================
Secondary Data Structures
=========================

.. warning::
	These secondary functions, mostly from the **utils** folder,= 
	**do not contain Type or Value Input Checks**.  Therefore, they are being
	listed with a warning that only advanced users should use these functions, as  
	they were not designed with input error checks.  

	Advanced users who wish to use these functions can find them 
	`here <https://github.com/GOMC-WSU/MoSDeF-dihedral-fit/tree/main/mosdef_dihedral_fit/utils>`_. 


Read and Write Functions (For Advanded Users)
---------------------------------------------

.. warning::
	These secondary functions, mostly from the **utils** folder,= 
	**do not contain Type or Value Input Checks**.  Therefore, they are being
	listed with a warning that only advanced users should use these functions, as  
	they were not designed with input error checks.  

	Advanced users who wish to use these functions can find them 
	`here <https://github.com/GOMC-WSU/MoSDeF-dihedral-fit/tree/main/mosdef_dihedral_fit/utils>`_. 


.. automodule:: mosdef_dihedral_fit.utils.file_read_and_write
	:members: get_atom_names_and_elements_from_mol2, get_atom_names_and_elements_from_pdb, write_xyz_file_from_gaussian_coordinates, write_restart_coor_from_xyz_file, check_gaussian_angle_energy_file_correct, check_gaussian_optimized_coordinate_file_correct, get_final_gaussian_output_file_data, get_gaussian_log_file_data, write_qm_data_files, get_matching_dihedral_info_and_opls_fitting_data
	

Math and Operation Functions (For Advanded Users)
-------------------------------------------------

.. warning::
	These secondary functions, mostly from the **utils** folder,= 
	**do not contain Type or Value Input Checks**.  Therefore, they are being
	listed with a warning that only advanced users should use these functions, as  
	they were not designed with input error checks.  

	Advanced users who wish to use these functions can find them 
	`here <https://github.com/GOMC-WSU/MoSDeF-dihedral-fit/tree/main/mosdef_dihedral_fit/utils>`_. 


.. automodule:: mosdef_dihedral_fit.utils.math_operations
	:noindex: dihedral_angle
	:members: round_to_sig_figs, normalize_vector, angle_between_2_vectors, check_previous_qm_values_match, sum_opls_const_1_plus_or_minus_cos_n_values, periodic_dihedral_n_1_2_3_4_5, RB_torsion_n_1_2_3_4_5, opls_dihedral_n_1_2_3_4
	
