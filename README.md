# MoSDeF-dihedral-fit
The MoSDeF-dihedral-fit: an open-source, transparent, and lightweight Python software package capable
of fitting dihedrals with QM calculations for existing forces fields. This software fits the
Optimized Potentials for Liquid Simulations (OPLS) <https://pubs.acs.org/doi/10.1021/ja9621760>_ style
dihedrals, then also analytically converts them to the periodic/CHARMM <https://www.charmm.org> dihedral and
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
        zeroed_dihedral_atom_types=None,
        qm_engine="gaussian",
        VDWGeometricSigma=False,
        atom_type_naming_style="general",
        gomc_cpu_cores=1,
        fit_min_validated_r_squared=0.99,
        fit_validation_r_squared_rtol=1e-03,
    )
import os
os.system("cat RB_torsion_k_constants_fit_energy.txt")
```

## Installation/Setup
```
conda install -c conda-forge mosdef_dihedral_fit
git clone https://github.com/GOMC-WSU/GOMC.git
cd GOMC
chmod u+x metamake.sh
./metamake.sh
./GOMC_<CPU|GPU>_XXXX +p4 in.conf # set to 4 threads
```

## Documentation

## Examples
Some basic workflows that use this package</br>
    - [ethane dihedral](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDeF-dihedral-fit/ethane_HC_CT_CT_HC)</br>
    - [propanoic acid dihedral](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDeF-dihedral-fit/protonated_fragment_CT_CT_C_OH)</br>

## Resources
This package is made as an API with [MoSDeF](https://github.com/mosdef-hub) and [MoSDeF-GOMC](https://github.com/GOMC-WSU/MoSDeF-GOMC). For `mosdef_dihedral_fit` to function, the forcefield files must be in a supported MoSDeF format (preferably the GMSO force field format), and use MoSDeF-GOMC and GOMC to perform the simulation setup and simulations. 
    - [MoSDeF tools](https://mosdef.org)
    - [MoSDeF-GOMC integration](https://mosdef-gomc.readthedocs.io/en/latest/index.html)
    - [GOMC Examples](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDef-GOMC)
    - [MoSDeF GMSO Sample Forcefields](https://github.com/mosdef-hub/gmso/tree/main/gmso/utils/files/gmso_xmls/test_ffstyles)

## Citations
    - Please cite MoSDeF-GOMC [here.](https://mosdef-gomc.readthedocs.io/en/latest/reference/citing_mosdef_gomc_python.html)
