[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/GOMC-WSU/MoSDeF-dihedral-fit/main.svg)](https://results.pre-commit.ci/latest/github/GOMC-WSU/MoSDeF-dihedral-fit/main)
[![CI](https://github.com/GOMC-WSU/MoSDeF-dihedral-fit/actions/workflows/CI.yml/badge.svg)](https://github.com/GOMC-WSU/MoSDeF-dihedral-fit/actions/workflows/CI.yml)

# MoSDeF-dihedral-fit
The MoSDeF-dihedral-fit: an open-source, transparent, and lightweight Python software package capable
of fitting dihedrals with QM calculations for existing forces fields. This software fits the [Optimized Potentials for Liquid Simulations (OPLS)](https://pubs.acs.org/doi/10.1021/ja9621760) style
dihedrals, then also analytically converts them to the periodic dihedral and
Ryckaert-Bellemans (RB) torsion.

```python
import unyt as u
from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import fit_dihedral_with_gomc
fit_dihedral_with_gomc(
        ["HC", "CT", "CT", "HC"], #atomclass names of the dihedral
        "molecule.mol2" # mol2 file with relevant structure
        "compound.xml", # xml file with other atomtyped parameters in foyer format
        298.15 * u.Kelvin, # relevant temperature
        gomc_binary_directory, # path to binary command from GOMC install
        {"HC_CT_CT_HC_multiplicity_1.log": []}, # log file to store info
        zero_dihedral_atom_types=None,
        qm_engine="gaussian",
        combining_rule='lorentz',
        atom_type_naming_style='general',
        gomc_cpu_cores=1,
        r_squared_min=0.99,
        r_squared_rtol=1e-03,
        opls_force_k0_zero=True
    )
import os
os.system("cat RB_torsion_k_constants_fit_energy.txt")
os.system("cat opls_torsion_k_constants_fit_energy.txt")
os.system("cat periodic_torsion_k_constants_fit_energy.txt")
```

### The plotted dihedral fits:
    - "opls_all_single_fit_dihedral_k_constants_figure.pdf"
    - "opls_all_summed_dihedrals_k_constants_figure.pdf"


## Quick Installation/Setup

**Note the mamba is used as a drop-in replacement for conda in the below installation, but conda works as well.**
```bash
mamba install -c conda-forge mosdef-dihedral-fit
conda activate mosdef-dihedral-fit
cd $CONDA_PREFIX
git clone https://github.com/GOMC-WSU/GOMC.git --branch v2.75a
cd GOMC
chmod u+x metamake.sh NVT
./metamake.sh
ln -s $CONDA_PREFIX/GOMC/bin/GOMC_CPU_NVT $CONDA_PREFIX/bin
```

GOMC will also be removed when the mosdef_dihedral_fit environment is removed. The installation can be placed anywhere though, that path will just have to manually be passed for the variable gomc_binary_directory.

#### More setup information
For more complete setup information see the [full installation documentation](https://github.com/GOMC-WSU/MoSDeF-dihedral-fit/blob/main/docs/getting_started/installation/installation.rst#installation)


## Documentation

The MoSDeF-dihedral-fit documentation can be found [here](https://mosdef-dihedral-fit.readthedocs.io/en/latest/).

## Dihedral Equations

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

$$U_{Periodic} = K_0 * (1 + cos(n_0*\theta - d_0))$$

$$+ K_1 * (1 + cos(n_1*\theta - d_1)) + K_2 * (1 + cos(n_2*\theta - d_2))$$

$$+  K_3 * (1 + cos(n_3*\theta - d_3)) +  K_4 * (1 + cos(n_4*\theta) - d_4)$$

$$where:  n_0 = 0  ;  n_1 = 1  ;  n_2 = 2  ;  n_3 = 3  ;  n_4 = 4 $$

$$d_0 = 90^o  ;  d_1 = 180^o  ;  d_2 = 0^o  ;  d_3 = 180^o  ;  d_4 = 0^o

## Examples
Some basic workflows that use this package and cover several force field (FF) types are [here](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDeF-dihedral-fit).  A summary of the FFs covered in the examples are listed below.

 - <b>OPLS-AA</b>: ethane dihedral
 - <b>OPLS-AA</b>: propanoic acid dihedral
 - <b>AMBER-AA</b>: butane dihedral
 - <b>TraPPE-UA</b>: butane dihedral
 - <b>Mie-UA</b>: butane dihedral
 - <b>Exp6-AA</b>: butane dihedral

## Resources
This package is made as an API with [MoSDeF](https://github.com/mosdef-hub) and [MoSDeF-GOMC](https://github.com/GOMC-WSU/MoSDeF-GOMC). For `mosdef_dihedral_fit` to function, the forcefield files must be in a supported MoSDeF format (preferably the GMSO force field format), and use MoSDeF-GOMC and GOMC to perform the simulation setup and simulations.

 - [MoSDeF tools](https://mosdef.org)
 - [MoSDeF-GOMC integration](https://mosdef-gomc.readthedocs.io/en/latest/index.html)
 - [GOMC Examples](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDef-GOMC)
 - [MoSDeF GMSO Sample Forcefields](https://github.com/mosdef-hub/gmso/tree/main/gmso/utils/files/gmso_xmls/test_ffstyles)

## Citations

 - Please cite MoSDeF-GOMC [here](https://mosdef-gomc.readthedocs.io/en/latest/reference/citing_mosdef_gomc_python.html)
 - Other tools used in this package can be found in the MoSDeF-Dihedral-Fit documentation.
