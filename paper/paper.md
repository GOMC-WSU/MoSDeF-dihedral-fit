
---
title: 'MoSDeF-dihedral-fit: A simple software package to fit dihedrals via the MoSDeF software'

tags:
  - Python
  - Molecular simulations
  - Molecular mechanics
  - Molecular dynamics
  - Monte Carlo
  - Quantum mechanics
  - dihedral fitting
  - torsion fitting
  - force field
  - MoSDeF
  - GOMC
  - MoSDeF-GOMC

authors:
  - name: Brad Crawford
    orcid: 0000-0003-0638-7333
    equal-contrib: true
    affiliation: "1, 2"
  - name: Co D. Quach
    orcid: 0000-0002-1255-4161
    affiliation: "3, 4"
  - name: Nicholas C. Craven
    orcid: 0000-0002-4607-4377
    affiliation: "4, 5"
  - name: Christopher R. Iacovella
    orcid: 0000-0003-0557-0427
    affiliation: "3, 4, 5"
  - name: Clare McCabe
    orcid: 0000-0002-8552-9135
    affiliation: "3, 4"
  - name: Peter T. Cummings
    orcid: 0000-0002-9766-2216
    affiliation: "3, 4"
  - name: Jeffrey J. Potoff
    orcid: 0000-0002-4421-8787
    equal-contrib: true
    affiliation: 2

affiliations:
 - name: Atomfold LLC, PA, USA
   index: 1
 - name: Department of Chemical Engineering, Wayne State University, Detroit, MI 48202-4050, USA
   index: 2
 - name: Department of Chemical and Biomolecular Engineering, Vanderbilt University, Nashville, TN 37235-1604, USA
   index: 3
 - name: Multiscale Modeling and Simulation (MuMS) Center, Vanderbilt University, Nashville, TN 37212, USA
   index: 4
 - name: Interdisciplinary Material Science Program, Vanderbilt University, Nashville, TN 37235-0106, USA
   index: 5

date: 31 December 2023
bibliography: paper.bib

---

# Summary

Molecular Mechanics (MM) simulations (molecular dynamics and Monte Carlo) provide a third method of scientific discovery, simulation modeling, adding to the traditional theoretical and experimental scientific methods.  These molecular simulations provide visualizations and calculated properties that are difficult, too expensive, or unattainable by conventional methods.  Additionally, molecular simulations can be utilized to obtain insights and properties on chemicals or materials that do not currently exist, are not easily attainable, or require hard-to-achieve conditions (i.e., very high pressures and temperatures).  However, these MM models operate from force field parameters determined from Quantum Mechanics (QM) simulations, where the MM proper dihedrals (i.e., dihedrals) are the most difficult to obtain if they don't currently exist for the chosen force field. These MM dihedrals are also not easily transferable between different force fields.  Accurate dihedrals are required to produce the correct energy landscape for the molecules with the proper molecular orientations in the simulations. These dihedrals can be harder to manually fit without software if there are multiple instances of the same dihedral in a molecule, or some cosine power series value should removed from the dihedral fit due to molecular symmetry.

The `MoSDeF-Dihedral-Fit` [@Crawford:2023b] library lets users quickly calculate the MM proper dihedrals (dihedrals) directly from the QM simulations for several force fields (OPLS, CHARMM, TraPPE, AMBER, Mie, and Exp6) [@Jorgensen:1996; @Brooks:2009; @Lee:2016-CG; @Martin:1998; @Weiner:1984; @Weiner:1986; @Mie:1903; @Buckingham:1938].  The user simply has to generate or use an existing Molecular Simulation Design Framework (MoSDeF) force field XML file [@Cummings:2021; @Summers:2020; @GMSO:2019; @forcefield-utilities:2022], provide a Gaussian 16 or Gaussian-style Quantum Mechanics (QM) simulation files that cover the dihedral rotation (typically, 0-360 degrees) and provide the molecular structure information in a mol2 format [@Gaussian16:2016].  The `MoSDeF-Dihedral-Fit` software utilizes the QM and MM data to fit the dihedral to the specific force field, fitting the constants for the OPLS dihedral equation form with the correct combining rules and 1-4 scaling factors, as specified in the MoSDeF XML force field file. This software also accounts for multiple instances of the dihedral and the molecular symmetry in the dihedral, automatically removing all the unusable cosine power series values due to this symmetry. The user can also set other dihedral energies in the molecule to zero, allowing for a more flexible and accurate dihedral fit; this allows the multiple dihedral's conformational energies to be calculated from a single dihedral angle, a trick that was utilized in some of the original OPLS dihedral fits.  The `MoSDeF-Dihedral-Fit` software then analytically calculates the Ryckaert-Bellemans (RB)-torsions and the periodic dihedral from the OPLS dihedral.  If another dihedral equation form is needed, the software outputs the raw data points to fit any other dihedral form.  Therefore, the `MoSDeF-Dihedral-Fit` software allows the fitting of any dihedral form, provided the force fields are supported by MoSDeF and MoSDeF-GOMC (uses the GPU Optimized Monte Carlo - GOMC MM engine) software [Crawford:2023a; Crawford:2022; @Crawford:2023b; Nejahi:2019; Nejahi:2021], and the QM data is provided as a Gaussian output file or a generalized Gaussian-style output form [Gaussian16:2016].


# Statement of need

Many different types of Molecular Mechanics (MM) simulation models exist, also called force fields.  While many of these force field parameters can be transferred between force fields, such as bonds, angles, and improper dihedrals (impropers), the proper dihedrals (dihedrals) can not be easily transferred between force fields due to the different combining rules (arithmetic and geometric) and 1-4 scaling factors (i.e., scaling factors between the 1st and 4th atoms) differing between the force fields [@Berthelot:1898; @Good:1970; @Lorentz:1881]. The accuracy of these dihedral parameters for each force field is critical for the molecular simulation to obtain the correct molecular configuration and orientations, which is absolutely required for understanding and analyzing the system's orientation, behavior, and properties (Example: free energy calculations, viscosity, adsorption concentrations, diffusion constants, and many more).

While some dihedral fitting software currently exists, they are not generalized and only fit the CHARMM-style force fields [@Mayne:2013], or only fit the dihedral constants to the final MM and QM energies that need to be calculated by other means [@Guvench:2008].  Therefore, the molecular simulation community needs a generalized software package that imports the QM and MM files, automatically reads and organizes the QM data, calculates the MM energies, and automatically fits the dihedral.  Additionally, the molecular simulation community needs software that fits the dihedral to any force field style and auto-corrects the fit to account for multiple instances of the dihedral and molecular symmetry since fitting these dihedrals is a high barrier to simulating new chemistry and materials if these parameters do not exist for the chosen force field.  The `MoSDeF-dihedral-fit` software accomplishes all this and automatically accounts for and fits the dihedral for any common combining rules, allows the zeroing out any other dihedrals in the molecule (i.e., allowing a better fit and user control by setting the other dihedral energies to zero), and accounts for any 1-4 scaling factors used via the MoSDeF XML files [Cummings:2021; Summers:2020; GMSO:2019; forcefield-utilities:2022], which contain the force fields. The `MoSDeF-dihedral-fit` [@Crawford:2023b] API fills the missing gap by providing a generalized and easy solution to fitting dihedrals for any dihedral form that is allowable in the Molecular Simulation Design Framework (MoSDeF) and MoSDeF-GOMC (uses the GPU Optimized Monte Carlo - GOMC MM engine) software [Crawford:2023a; Crawford:2022; @Crawford:2023b; Nejahi:2019; Nejahi:2021].

The publicly available force field dihedrals (Example: OPLS) were originally fit when the QM simulations used to fit the dihedrals were computationally prohibitive.  Due to these limitations, scientists assumed that the dihedral fits were transferable with all the atom classes in the dihedral fit; however, this is not always an accurate assumption.  Some of the dihedrals were only fit to the first minimum and not the entire dihedral landscape, which can create issues when the simulation does not accurately calculate the whole energy landscape.  These prior assumptions in the dihedral fits can also lead to problems in reproducibility between force fields.  Today, with more advanced hardware and software, QM simulations can be conducted with more complex molecules, allowing for higher quality and customized dihedral fits.  The `MoSDeF-dihedral-fit` software will enable scientists to create a general or custom dihedral quickly and accurately for any force field, providing the new fit or the entire and correct energy landscape to an existing fit (i.e., producing the correct molecular orientations), and if applied correctly, this can make simulations using different forcefield more reproducible.


# Acknowledgements

This research was partially supported by the National Science Foundation (grants OAC-1835713, OAC-1835874, and CBET 2052438).  Atomfold LLC also donated research and development time and computational resources for this research and software.  Wayne State University Grid provided some of the computational resources used in this work.

# Mathematics

Proper dihedral (dihedral) forms.


<u>OPLS-dihedral</u>:

$$OPLS_{Energy} = \frac{f_0}{2}$$

$$+ \frac{f_1}{2} * (1 + cos(\theta)) + \frac{f_2}{2} * (1-cos(2 * \theta))$$

$$+ \frac{f_3}{2} * (1 + cos(3 * \theta)) + \frac{f_4}{2}  *(1-cos(4 * \theta))$$

<u>Ryckaert-Bellemans (RB)-torsions</u>:

$$RB_{Energy} = C_0$$

$$+ C_1 * cos(\psi) + C_2 * cos(\psi)^2$$

$$+ C_3 * cos(\psi)^3 + C_4 * cos(\psi)^4$$

$$\psi = \theta - 180^o$$

<u>Periodic-dihedral</u>:

$$Periodic_{Energy} = K_0 * (1 + cos(n_0*\theta - 90^o))$$

$$+ K_1 * (1 + cos(n_1*\theta - 180^o)) + K_2 * (1 + cos(n_2*\theta))$$

$$+  K_3 * (1 + cos(n_3*\theta - 180^o)) +  K_4 * (1 + cos(n_4*\theta))$$
