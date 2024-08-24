
---
title: 'MoSDeF-dihedral-fit: A lightweight software for fitting dihedrals within MoSDeF'

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
 - name: Atomfold, PA, USA
   index: 1
 - name: Department of Chemical Engineering and Materials Science, Wayne State University, Detroit, MI 48202-4050, USA
   index: 2
 - name: Department of Chemical and Biomolecular Engineering, Vanderbilt University, Nashville, TN 37235-1604, USA
   index: 3
 - name: Multiscale Modeling and Simulation (MuMS) Center, Vanderbilt University, Nashville, TN 37212, USA
   index: 4
 - name: Interdisciplinary Material Science Program, Vanderbilt University, Nashville, TN 37235-0106, USA
   index: 5

date: 24 August 2024
bibliography: paper.bib

---

# Summary

Molecular Mechanics (MM) simulations (molecular dynamics and Monte Carlo) provide a third method of scientific discovery, simulation modeling, adding to the traditional theoretical and experimental scientific methods.  These molecular simulations provide visualizations and calculated properties that are difficult, too expensive, or unattainable by conventional methods.  Additionally, molecular simulations can be used to obtain insights and properties on chemicals or materials that do not currently exist, are not easily attainable, or require hard-to-achieve conditions (i.e., very high pressures and temperatures).  However, these MM models require force field parameters to be determined ideally from Quantum Mechanics (QM) simulations or other methods, including the vibrational spectrum and machine learning methods [@Kania:2021; @Friederich:2018; @Vermeyen:2023; @Mayne:2013; @Schmid:2011], where the MM proper dihedrals (i.e., dihedrals) are challenging to obtain if they don't currently exist for the chosen force field or they do not properly scale-up to larger molecules or moiety combinations from their separately derived parameters using small molecules [@Kania:2021: @Mayne:2013]. While the same QM simulations can be used to fit the dihedrals in all of the force field types, these MM dihedrals are also not easily transferable between different force fields due to the differing parameters and formulas used in each force field (e.g., Combining rules and 1-4 scaling factors).


The `MoSDeF-Dihedral-Fit` [@Crawford:2023b] library lets users quickly calculate the MM proper dihedrals (dihedrals) directly from the QM simulations for several force fields (OPLS, TraPPE, AMBER, Mie, and Exp6) [@Jorgensen:1996; @Martin:1998; @Weiner:1984; @Weiner:1986; @Potoff:2009; @Hemmen:2015; @Errington:1999].  The user simply has to generate or use an existing Molecular Simulation Design Framework (MoSDeF) force field XML file [@Cummings:2021; @Summers:2020; @GMSO:2019; @forcefield-utilities:2022], provide Gaussian 16 or Gaussian-style Quantum Mechanics (QM) simulation files that cover the dihedral rotation (typically, 0-360 degrees), and provide the molecular structure information in a mol2 format [@Gaussian16:2016].  The `MoSDeF-Dihedral-Fit` software uses the QM and MM data to fit the dihedral for the specific force field, fitting the constants for the OPLS dihedral equation form with the correct combining rules and 1-4 scaling factors, as specified in the MoSDeF XML force field file. This software also accounts for multiple instances of the dihedral and the molecular symmetry in the molecule, and automatically removes all the unusable cosine power series values due to this symmetry.  The user can set other dihedral energies in the molecule to zero, allowing for a more flexible and accurate dihedral fit; this allows the multiple dihedral's conformational energies to be calculated from a single dihedral angle, a strategy that was used in some of the original OPLS dihedral fits.  The `MoSDeF-Dihedral-Fit` software analytically calculates the Ryckaert-Bellemans (RB)-torsions and the periodic dihedral from the OPLS dihedral.  If another form of the dihedral equation that is not currently supported is needed, the software outputs the raw data points to enable users to fit any other dihedral form.  Therefore, the `MoSDeF-Dihedral-Fit` software allows the fitting of any dihedral form, provided the force fields and software it utilizes are supported by MoSDeF and MoSDeF-GOMC (which uses GPU Optimized Monte Carlo - GOMC) [@Crawford:2023a; @Crawford:2022; @Crawford:2023b; @Nejahi:2019; @Nejahi:2021], vmd-python [@vmd-python:2016] (a derivative or modified version of the original VMD software [@Humphrey:1996; @Stone:2001]), and the QM data is provided as a Gaussian output file, or a generalized Gaussian-style output form [@Gaussian16:2016].


# Statement of need

Many different types of Molecular Mechanics (MM) simulation models exist, which are also known as "force fields".  While many of these force field parameters can be transferred between force fields, such as bonds, angles, and improper dihedrals (impropers), the proper dihedrals (dihedrals) can not be easily transferred due to the different combining rules (arithmetic and geometric) and 1-4 scaling factors (i.e., scaling factors between the 1st and 4th atoms) that were used in the development of the original parameters [@Berthelot:1898; @Good:1970; @Lorentz:1881]. The accuracy of these dihedral parameters for each force field is critical for molecular simulations to obtain the correct molecular conformations and configurations, which are absolutely required for understanding and analyzing the system's microstructure and physical properties (e.g., free energies, viscosities, adsorption loading, diffusion constants, and many more).

While some dihedral fitting software currently exists, they only fit the CHARMM-style force fields [@Mayne:2013], or fit the dihedral constants to the final MM and QM energies, which need to be calculated by other means [@Guvench:2008].  Therefore, the molecular simulation community needs a generalized software package that imports QM and MM files, automatically reads and organizes the QM data, calculates the MM energies, and automatically fits the dihedral.  Additionally, the molecular simulation community needs software that fits the dihedral to any force field style and auto-corrects the fit to account for multiple instances of the dihedral and molecular symmetry, since fitting these dihedrals is a high barrier to simulating new chemistry and materials if these parameters do not exist for the chosen force field.  The `MoSDeF-dihedral-fit` software accomplishes all this and automatically accounts for any of the common combining rules, allows the zeroing out any other dihedrals in the molecule (i.e., allowing a better fit and user control by setting the other dihedral energies to zero), and accounts for any 1-4 scaling factors specified via the MoSDeF XML files [@Cummings:2021; @Summers:2020; @GMSO:2019; @forcefield-utilities:2022], which contain the force fields. For example, two dihedrals in the same rotation cycle that both start at the carbon chain and end with a carboxylic acid demonstrate that setting one dihedral to zero is a typical practice, as the C-C-C-O: (O: = oxygen without hydrogen) dihedral is set to zero while the C-C-C-O (O = oxygen with hydrogen) is the only fit or non-zero dihedral [@Jorgensen:1996; @Chao:2021; @Ganesh:2004]. Setting one dihedral to zero avoids potentially incorrect forces and a problematic or bad cosine series dihedral fit; otherwise, it would require fitting the first dihedral with a keytone or alcohol and then fitting the remaining dihedral in the carboxylic acid.  The `MoSDeF-dihedral-fit` [@Crawford:2023b] API fills the missing gap by providing a generalized and easy solution to fitting dihedrals for any dihedral form that is allowable in the Molecular Simulation Design Framework (MoSDeF) and MoSDeF-GOMC (uses the GPU Optimized Monte Carlo - GOMC MM engine) software [@Crawford:2023a; @Crawford:2022; @Crawford:2023b; @Nejahi:2019; @Nejahi:2021].

Commonly used force field dihedrals, such as OPLS, were originally fit when the QM simulations used to fit the dihedrals were computationally prohibitive.  Due to these limitations, scientists assumed that the dihedral fits were transferable with all the atom classes in the dihedral fit; however, this is not always an accurate assumption.  Some of the dihedrals were only fit to the first minimum and not the entire dihedral landscape, which can lead to errors in the predicted molecular conformations.  These prior assumptions in the dihedral fits may also lead to problems in reproducibility in modified force fields.  Today, with more advanced hardware and software, QM simulations can be conducted with more complex molecules, allowing for higher quality and customized dihedral fits.  The `MoSDeF-dihedral-fit` software will enable scientists to create a generalized, or molecule-specific dihedral parameters, quickly, accurately and reproducibly for any force field.

# Acknowledgements

This research was partially supported by the National Science Foundation (grants OAC-1835713, OAC-1835874, and CBET 2052438).  Atomfold LLC also donated research and development time and computational resources for this research and software.  Wayne State University Grid provided some of the computational resources used in this work.

# Mathematics

Proper dihedral (dihedral) forms.

<u>OPLS-dihedral</u>:

$$U_{OPLS} = \frac{k_0}{2}$$

$$+ \frac{k_1}{2} * (1 + cos(\theta)) + \frac{k_2}{2} * (1-cos(2 * \theta))$$

$$+ \frac{k_3}{2} * (1 + cos(3 * \theta)) + \frac{k_4}{2}  *(1-cos(4 * \theta))$$

<u>Ryckaert-Bellemans (RB)-torsions</u>:

$$U_{RB} = C_0$$

$$+ C_1 * cos(\psi) + C_2 * cos(\psi)^2$$

$$+ C_3 * cos(\psi)^3 + C_4 * cos(\psi)^4$$

$$\psi = \theta - 180^o$$

<u>Periodic-dihedral</u>:

$$U_{Periodic} = K_0 * (1 + cos(n_0*\theta - 90^o))$$

$$+ K_1 * (1 + cos(n_1*\theta - 180^o)) + K_2 * (1 + cos(n_2*\theta))$$

$$+  K_3 * (1 + cos(n_3*\theta - 180^o)) +  K_4 * (1 + cos(n_4*\theta))$$
