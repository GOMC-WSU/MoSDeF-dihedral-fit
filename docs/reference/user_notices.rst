============
User Notices
============

There are some critical items to complete and consider when using the **MoSDeF-dihedral-fit** package:

    #. The **MoSDeF** XML file, or force field file, must be of an existing force field and have the dihedral parameters entered for the selected dihedral.  The selected force field also must be compatible with **MoSDeF-GOMC**.  These dihedral parameters should all be set all zeros, but if not, this program will automatically zero  the selected dihedral for the fitting process.

    #. All the parameters in the **MoSDeF** XML file must be correct to ensure a correct dihedral fit, aside from the dihedral which is being fit.  This is entirely up to the user to determine the correctness of the **MoSDeF** XML file, as this software has no way to determine the correct values of a force field.

    #. The new **MoSDeF-GOMC** version using `GMSO <https://gmso.mosdef.org/en/stable/>`_ allows the users to enter the force field parameter in their original form and units, which minimizes user errors.  Therefore, this **GMSO** force field file form is always recommended when creating a new **MoSDeF** XML file.

    #. The mixing/combining rules and 1-4 non-bonded dihedral interactions are in the **MoSDeF** XML file and will automatically be input into **GOMC**, the Molecular Mechanics (MM) engine, so please be sure these are entered properly in the XML file or it may set the default parameters, which will produce a wrong result.

    #. The atomic order and molecule must exactly match for the user inputted `Gaussian <https://gaussian.com>`_ log files and the **mol2** file.  There are checks to ensure that the same elements and the number of atoms are used, but there is no simple way to determine if the atoms are out of order when these criteria match.  Therefore, it is up to the user to ensure these are the same.  It is recommended that the user create the **mol2** and use its elements and coordinates to set up the **Gaussian** simulation, as this will reduce the chance of errors.
