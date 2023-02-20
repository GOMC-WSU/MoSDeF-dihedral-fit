
General Information
===================
.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT


MoSDeF-dihedral-fit Basics
--------------------------
The **MoSDeF-dihedral-fit** program is an open-source, transparent, and lightweight Python software package capable
of fitting dihedrals for existing forces fields.  This software fits the
`Optimized Potentials for Liquid Simulations (OPLS) <https://pubs.acs.org/doi/10.1021/ja9621760>`_ style
dihedrals, then also analytically converts them to the periodic/ `CHARMM <https://www.charmm.org>`_ dihedral and
Ryckaert-Bellemans (RB) torsion.

The **MoSDeF-dihedral-fit** software produces dihedral fits for existing force fields, compatible with
`GPU Optimized Monte Carlo (GOMC) <http://gomc.eng.wayne.edu>`_ and
`Molecular Simulation Design Framework (MoSDeF) <https://mosdef.org>`_, with only tens of lines of python code,
the `Gaussian software <https://www.gaussin.com>`_ log files and a **mol2** file.

.. note::
    Currently, this means that only the fourth (4th) cosine multiple or power is utilized in the dihedral fit.

.. note::
    If using the **CHARMM** dihedral format, the K0 term (constant term)
    will have to be disregarded because the **CHARMM** recognizes this as a harmonic dihedral,
    not a periodic dihedral.


The `MoSDeF-GOMC <https://github.com/GOMC-WSU/MoSDeF-GOMC/tree/master/mosdef_gomc>`_ software package is used
for the Molecular Mechanics (MM) calculation, which utilizes the
`GPU Optimized Monte Carlo (GOMC) <http://gomc.eng.wayne.edu>`_, the
`Molecular Simulation Design Framework (MoSDeF) <https://mosdef.org>`_, and the
`vmd-python <https://github.com/Eigenstate/vmd-python>`_ core software packages. For the Quantum Mechanics calculations,
the `Gaussian software <https://gaussian.com>`_ software is used, reading the **Gaussian** log files.
From the **Gaussian** log file, a user created **mol2** file, and a few user inputs, this software automatically
fits the desired dihedral, accounting for multiple dihedrals simultaneously. The software output provides
a wide range of allowable dihedral fits with different cosine term combinations, including plots for visual reference;
from this information, the user can then select the best dihedral fit for the specific application.
Additionally, the 1-4 interactions for the force fields can be explicitly set in the force field XML file,
allowing the flexibility that some other dihedral fitters lack.


MoSDeF-dihedral-fit is a part of the MoSDeF ecosystem
-----------------------------------------------------
The **MoSDeF-dihedral-fit** software acts as manager between
`GPU Optimized Monte Carlo (GOMC) <http://gomc.eng.wayne.edu>`_, the
`Molecular Simulation Design Framework (MoSDeF) <https://mosdef.org>`_, and the
`vmd-python <https://github.com/Eigenstate/vmd-python>`_ core software packages,
providing a user-friendly package that automatically produces one or more reproducible
dihedral fits, with visualizations, which the user can then select the best fit for their system.


The **MoSDeF** software consists of the following core packages:
	* `mBuild <https://mbuild.mosdef.org/en/stable/>`_ -- A hierarchical, component based molecule builder

	* `foyer <https://foyer.mosdef.org/en/stable/>`_ -- A package for atom-typing as well as applying and disseminating forcefields

	* `GMSO <https://gmso.mosdef.org/en/stable/>`_ -- Flexible storage of chemical topology for molecular simulation

    * `forcefield-utilities <https://github.com/mosdef-hub/forcefield-utilities/>`_
