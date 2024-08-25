
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

Molecular Mechanics (MM) simulations (e.g., molecular dynamics and Monte Carlo) provide a third method of scientific discovery, adding to the traditional theoretical and experimental scientific methods [@Mielke:2019; @Siegfried:2014].  Experimental methods measure the data under set conditions (e.g., temperature and pressure), whereas the traditional theoretical methods are based on analytical equations, and sometimes their constants are fitted to experimental data.  The MM simulations are deterministic and stochastic, and their models, commonly known as "force fields", can be optimized to match experimental data, similar to analytical theory-based methods [@Allen:2017; @Frekel:2002; @Jorgensen:1996; @Martin:1998; @Weiner:1984; @Weiner:1986; @Potoff:2009; @Hemmen:2015; @Errington:1999].  In larger, more complex systems, the stochastic simulation's molecules can jump large energy barriers that deterministic simulations may not be able to overcome in a reasonable timeframe, even with modern computing capabilities [@Allen:2017; @Frekel:2002].  However, deterministic and stochastic systems that provide adequate sampling for calculating a given property can provide critical insights into the system's phase space, which are not obtainable via traditional theoretical and experimental methods.  Additionally, molecular simulations provide critical insights from visualizations and by obtaining chemical or material properties that do not currently exist, are not easily attainable (e.g., too expensive or dangerous) by traditional theoretical and experimental methods [@Hollingsworth:2018; @Hirst:2014], or require hard-to-achieve conditions, such as very high pressures and temperatures [@Yu:2023; @Koneru:2022; @Swai:2020; @Kumar:2022; @Louie:2021].  However, the force field parameters are ideally determined from Quantum Mechanics (QM) simulations or other methods, including the vibrational spectrum and machine learning methods [@Kania:2021; @Friederich:2018; @Vermeyen:2023; @Mayne:2013; @Schmid:2011; @Vanommeslaeghe:2014].  The MM proper dihedrals (i.e., dihedrals) are challenging to obtain if they don't currently exist for the chosen force field, inaccurately scale-up in larger molecules, or misbehave with other moiety combinations, provided some were separately derived using small molecules [@Kania:2021; @Mayne:2013].  While the same QM simulations can fit the dihedrals in most force field types, these dihedrals aren't easily transferable between force fields due to the differing parameters and formulas, including the combining rules and 1-4 scaling factors. [@Huang:2013; @Vanommeslaeghe:2010; @Vanommeslaeghe:2014; @Chen:2015].

The `MoSDeF-Dihedral-Fit` [@Crawford:2023b] library lets users quickly calculate the MM dihedrals directly from the QM simulations for several force fields (OPLS, TraPPE, AMBER, Mie, and Exp6) [@Jorgensen:1996; @Martin:1998; @Weiner:1984; @Weiner:1986; @Potoff:2009; @Hemmen:2015; @Errington:1999].  The user simply has to generate or use an existing Molecular Simulation Design Framework (MoSDeF) force field XML file [@Cummings:2021; @Summers:2020; @GMSO:2019; @forcefield-utilities:2022], provide Gaussian 16 log or Gaussian-style QM simulation files that cover the dihedral rotation (typically, 0-360 degrees), and provide the molecular structure information in a mol2 format [@Gaussian16:2016].  The `MoSDeF-Dihedral-Fit` software uses the QM and MM data to produce the dihedral for the specific force field, fitting the constants for the OPLS dihedral form and then analytically converting them to the Ryckaert-Bellemans (RB)-torsions and the periodic dihedrals.  The software outputs the calculated MM dihedral points, enabling users to fit unsupported dihedral forms, provided the force fields are supported by the MoSDeF, GPU Optimized Monte Carlo (GOMC),and MoSDeF-GOMC [@Crawford:2023a; @Crawford:2022; @Crawford:2023b; @Nejahi:2019; @Nejahi:2021], and vmd-python [@vmd-python:2016] software (a derivative of the VMD software [@Humphrey:1996; @Stone:2001]).


# Statement of need

While many of these Molecular Mechanics (MM) force field parameters can be transferred between force fields, such as bonds, angles, and improper dihedrals (impropers), the proper dihedrals (dihedrals) can not be easily transferred due to the different combining rules (arithmetic and geometric) and 1-4 scaling factors (i.e., between the 1st and 4th bonded atoms) that were used in the development of the original parameters [@Berthelot:1898; @Good:1970; @Lorentz:1881]. The accuracy of these dihedral parameters is critical in obtaining the correct molecular conformations and configurations, which are absolutely required for understanding and analyzing the system's microstructure and physical properties (e.g., free energies, viscosities, adsorption loading, diffusion constants, and many more).

While some dihedral fitting software currently exists, they only fit the CHARMM-style force fields [@Mayne:2013], or fit the dihedral constants to the final MM and QM energies, which need to be calculated by other means [@Guvench:2008].  Therefore, the molecular simulation community needs a generalized software package that imports QM and MM files, automatically reads and organizes the QM data, calculates the MM energies, auto-corrects the dihedral fit to account for multiple instances of the dihedral, and automatically removes the unusable cosine power series combinations due to this symmetry.  The `MoSDeF-dihedral-fit` software accomplishes all this and automatically accounts for any of the common combining rules and the 1-4 scaling factors specified via the MoSDeF XML (i.e., force field) files [@Cummings:2021; @Summers:2020; @GMSO:2019; @forcefield-utilities:2022].  By allowing the user to set any other dihedral in the molecule to zero, this software avoids forcing one dihedral fit to correct the inaccurate forces of another dihedral, resulting in a problematic or bad cosine series fit; thus, providing a more flexible and accurate fit by combining multiple dihedral conformational energies in a single dihedral, a strategy used in the original and modern OPLS force fields [Jorgensen:1996: @Chao:2021].  For example, a carboxylic acid with an alkyl tail has two dihedrals in the same rotation cycle; the C-C-C-O: (O: = oxygen without hydrogen) dihedral is set to zero while the C-C-O-H dihedral is fit [@Jorgensen:1996; @Chao:2021; @Ganesh:2004].  The `MoSDeF-dihedral-fit` [@Crawford:2023b] API fills the missing gap by providing a generalized and easy solution to fitting dihedrals in their commonly used forms and outputting the MM dihedral data points so users can fit other custom dihedral forms.


Commonly used force field dihedrals, such as OPLS, were originally fit when the QM simulations used to fit the dihedrals were computationally prohibitive.  Due to these limitations, scientists assumed that the dihedral fits were transferable with all the atom classes in the dihedral fit; however, this is not always an accurate assumption.  Some of the dihedrals were only fit to the first minimum and not the entire dihedral landscape, which can lead to errors in the predicted molecular conformations.  These prior assumptions in the dihedral fits may also lead to problems in reproducibility in modified force fields.  Today, with more advanced hardware and software, QM simulations can be conducted with more complex molecules, allowing for higher quality and customized dihedral fits.  The `MoSDeF-dihedral-fit` software will enable scientists to create a generalized, or molecule-specific dihedral parameters, quickly, accurately and reproducibly.

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
