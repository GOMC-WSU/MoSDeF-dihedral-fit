MoSDeF-dihedral-fit: A simple software package to fit dihedrals via the MoSDeF software
=======================================================================================


The **MoSDeF-dihedral-fit** program is an open-source, transparent, and lightweight Python software package capable
of fitting dihedrals for existing forces fields.  This software fits the
`Optimized Potentials for Liquid Simulations (OPLS) <https://pubs.acs.org/doi/10.1021/ja9621760>`_ style
dihedrals, then also analytically converts them to the periodic/`CHARMM <https://www.charmm.org>`_ dihedral and
Ryckaert-Bellemans (RB) torsion.

The **MoSDeF-dihedral-fit** software produces dihedral fits for existing force fields, compatible with
`GPU Optimized Monte Carlo (GOMC) >= v2.75 <http://gomc.eng.wayne.edu>`_ and
`Molecular Simulation Design Framework (MoSDeF) <https://mosdef.org>`_, with only tens of lines of python code,
the `Gaussian 16 <https://gaussian.com>`_ log files, and a **mol2** file.

.. note::
    Currently, this means that only the fourth (4th) cosine multiple or power is utilized in the dihedral fit.

.. note::
    If using the **CHARMM** dihedral format, the K0 term (constant term)
    will have to be disregarded because the **CHARMM** recognizes this as a harmonic dihedral,
    not a periodic dihedral.


The `MoSDeF-GOMC <https://github.com/GOMC-WSU/MoSDeF-GOMC/tree/master/mosdef_gomc>`_ software package is used
for the Molecular Mechanics (MM) calculation, which utilizes
`GPU Optimized Monte Carlo (GOMC) >= v2.75 <http://gomc.eng.wayne.edu>`_, the
`Molecular Simulation Design Framework (MoSDeF) <https://mosdef.org>`_, and the
`vmd-python <https://github.com/Eigenstate/vmd-python>`_ core software packages. For the Quantum Mechanics calculations,
the **Gaussian 16** software is used, reading the **Gaussian 16** log files.
From the **Gaussian 16** log file, a user created **mol2** file, and a few user inputs, this software automatically
fits the desired dihedral, accounting for multiple dihedrals simultaneously. The software output provides
a wide range of allowable dihedral fits with different cosine term combinations, including plots for visual reference;
from this information, the user can then select the best dihedral fit for the specific application.
Additionally, the 1-4 interactions for the force fields can be explicitly set in the force field XML file,
allowing the flexibility that some other dihedral fitters lack. Lastly, the dihedral fits are compared by recalculating
the dihedral in **GOMC >= v2.75** and comparing it to the original **Gaussian 16** energies, ensuring a correct dihedral fit.

**Dihedral Equations**:

OPLS-dihedral:

.. math::
    U_{OPLS} = \frac{k_0}{2}
                    + \frac{k_1}{2}*(1+cos(\theta))
                    + \frac{k_2}{2}*(1-cos(2*\theta))

.. math::
                        + \frac{k_3}{2}*(1+cos(3*\theta))
                        + \frac{k_4}{2}*(1-cos(4*\theta))

Ryckaert-Bellemans (RB)-torsions:

.. math::
    U_{RB} = C_0 + C_1*cos(\psi)
                  + C_2*cos(\psi)^2
                  + C_3*cos(\psi)^3
                  + C_4*cos(\psi)^4

.. math::
   \psi = \theta - 180^o

Periodic-dihedral:

.. math::
    U_{Periodic} = K_0 * (1 + cos(n_0*\theta - 90^o))

.. math::
                            + K_1 * (1 + cos(n_1*\theta - 180^o))
                            + K_2 * (1 + cos(n_2*\theta))

.. math::
                            + K_3 * (1 + cos(n_3*\theta - 180^o))
                            + K_4 * (1 + cos(n_4*\theta))

**MoSDeF-dihedral-fit Highlights**:
   #. With a **Gaussian 16** log file and a few user inputs, the user can easily fit a dihedral.

   #. **MoSDeF-dihedral-fit** is designed to automate the dihedral fit

   #. **MoSDeF-dihedral-fit** permits a transparent dihedral fitting process, where the user can alter the input force field XML file, which is the standard input for all of the **MoSDeF** software packages.

**MoSDeF-dihedral-fit with United-Atom (UA) force fields**:
   #.  For the QM All-Atom (AA) and UA atoms parameterization to work together, the UA atoms must be entered as single atoms in the sudo UA force field XML and mol2 files.
   #. The individual atoms in the force field XML and mol2 files used to fit the dihedral shall have force field parameters for the central atoms (i.e., C or CH_3 in CH3) and other non-UA force field atoms alike (atoms with zero nonbonded interactions in a UA force field). The user should enter zeros for the force field terms wherever possible for the atoms not part of the UA forcefield (Example: atoms with zero nonbonded interactions in a UA force field), as setting the non-central atoms parameters to zero is especially true and critically important for the nonbonded interactions.

**MoSDeF-dihedral-fit License**:

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

**vmd-python License**:

The **vmd-python** software is licensed by the Theoretical and Computational Biophysics Group at the Beckman Institute for Advanced Science and Technology at the University of Illinois at Urbana-Champaign, which is a modified version of VMD.

The vmd-python website is https://github.com/Eigenstate/vmd-python

The official VMD web page is http://www.ks.uiuc.edu/Research/vmd



.. toctree::
	:caption: Overview
    	:maxdepth: 2

    	overview/general_info

.. toctree::
	:caption: Getting Started
    	:maxdepth: 2

	getting_started/installation/installation
        getting_started/quick_start/quick_start

.. toctree::
	:caption: Topic Guides
    	:maxdepth: 2

    	topic_guides/data_structures
	topic_guides/secondary_structures

.. toctree::
    	:caption: Reference
    	:maxdepth: 2

        reference/units
        reference/user_notices
        reference/contributing
        reference/credits
        reference/citing_mosdef_dihedral_fit_python
