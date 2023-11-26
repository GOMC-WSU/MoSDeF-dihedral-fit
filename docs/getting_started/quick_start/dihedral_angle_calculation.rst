Calculate a Single Dihedral Angle
=================================

Calculate a single dihedral angle from the coordinates of the four (4) atoms or beads.

Import the required packages.

.. code:: ipython3

    import unyt as u
    from mosdef_dihedral_fit.utils.math_operations import dihedral_angle

Select the desired coordinates of each of the 4 atoms and beads, in order, and calculate the dihedral angle.

.. code:: ipython3

    # The atom or bead coordinates in the x, y, and z coordinates
    atom_xyz_coord_1 = [0.3256   -0.9455    0.0000]
    atom_xyz_coord_2 = [0.0000    0.0000    0.0000]
    atom_xyz_coord_3 = [-1.5000    0.0000    0.0000]
    atom_xyz_coord_4 = [-1.8256   -0.4728   -0.8188]

    dihedral_angle_output = dihedral_angle(
        atom_xyz_coord_1, 
        atom_xyz_coord_2,       
        atom_xyz_coord_3, 
        atom_xyz_coord_4
    )

    print(f"dihedral_angle_output = {dihedral_angle_output} degrees")