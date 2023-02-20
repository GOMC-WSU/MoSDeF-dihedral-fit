Select the Inputs and Fit the Dihedral
======================================
.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT


Fit the dihedral for Ethane
---------------------------
The ethane HC-CT-CT-HC dihedral fit, where nine (9) of the HC-CT-CT-HC dihedral are fit simultaneously.


.. note::
    The GOMC software need to be installed manually, outside of this Python install,
    with it's directory/path specified in the dihedral fit function.

Import the required mbuild package.

.. code:: ipython3

    import unyt as u
    from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import fit_dihedral_with_gomc

Select the desired variables, file, and set the temperature

.. code:: ipython3

    # The MoSDeF force field (FF) XML file which will be used
    FF_XML_file = 'path_to_file/oplsaa_ethane_HC_CT_CT_HC.xml'

    # The mol2 file which is used.
    mol2_file = 'path_to_file/ethane_aa.mol2'

    # The temperature of the Molecular Mechanics (MM) simulation.
    temperature_in_unyt_units = 298.15 * u.Kelvin

    # Choose mixing rule (override_VDWGeometricSigma_bool) to override use the bool (True or False)
    # that was chosen in the force field (FF) XML file.  This variable is not required and will be
    # selected automatically; however, you should override to the correct setting if unsure.
    override_VDWGeometricSigma_bool = True

    # Atom type naming convention ( str, optional, default=’all_unique’, (‘general’ or ‘all_unique’) )
    # General is safe and recommended since we are using a single FF XML file.
    atom_type_naming_style = 'general'

    # Load 1 or more Gaussian file (key), and the value ([0]) is a list of Gaussian point to remove
    # from the fitting process.  More than 1 Gaussian file can be loaded, allowing the user
    # to run multiple dihedral angles in separate file, minimizing the time required to run
    # the simimulaiton (i.e., the user can split them up into many simulations to obtain
    # the full dihedral rotation.)
    log_files_and_removed_points = {
        'path_to_file/HC_CT_CT_HC_multiplicity_1.log': [0],
    }

    # The dihedral which is being fit
    fit_dihedral_atom_types = ['HC', 'CT', 'CT', 'HC']

    # All the other dihedrals which can be zeroed in the fitting process, in a nexted list
    zeroed_dihedrals = None

Perform the dihedral fit using the Gaussian log file as the Quantum Mechanics (QM) engine
and GOMC Molecular Mechanics (MM) engine, and write out all the files fit via the standard
OPLS equation.  It also outputs the periodic dihedral format and the RB torsions format,
which are analytically converted from the standard OPLS dihedral fit.

.. code:: ipython3

    # run the "fit_dihedral_with_gomc" command
    fit_dihedral_with_gomc(
        fit_dihedral_atom_types,
        mol2_file,
        FF_XML_file,
        temperature_in_unyt_units,
        gomc_binary_directory,
        log_files_and_removed_points,
        zeroed_dihedral_atom_types=zeroed_dihedrals,
        qm_engine="gaussian",
        override_VDWGeometricSigma=override_VDWGeometricSigma_bool,
        atom_type_naming_style='general',
        gomc_cpu_cores=1,
        fit_min_validated_r_squared=0.99,
        fit_validation_r_squared_rtol=1e-03
    )

