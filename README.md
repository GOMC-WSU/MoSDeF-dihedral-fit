# MoSDeF-dihedral-fit
MoSDeF-dihedral-fit: A simple software package to fit dihedrals via the MoSDeF software.

This package is a tool to add dihedral parameters for [Ryckaert-Bellemans](https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-ryckaert-bellemans-function) dihedral based forcefields by fitting QM calculations for the rotation about the dihedral angle to this functional form. The main functionality is found through the following code snippet.
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
pytest -v
```

## Documentation

## Examples
Some basic workflows that use this package</br>
    - [ethane dihedral](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDeF-dihedral-fit/ethane_HC_CT_CT_HC)</br>
    - [propanoic acid dihedral](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDeF-dihedral-fit/protonated_fragment_CT_CT_C_OH)</br>

## Resources
This package is made as an API with [MoSDeF](https://github.com/mosdef-hub) and [MoSDeF-GOMC](https://github.com/GOMC-WSU/MoSDeF-GOMC). In order to get the most of the dihedral-fitter, it is recommended to do your forcefielding via MoSDeF and simulations via GOMC. Included are some information to get users familiar with those pacakges.
    - [MoSDeF tools](https://mosdef.org)
    - [MoSDeF-GOMC integration](https://mosdef-gomc.readthedocs.io/en/latest/index.html)
    - [GOMC Examples](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDef-GOMC)
    -

## Citations
    - Please cite MoSDeF-GOMC [here.](https://mosdef-gomc.readthedocs.io/en/latest/reference/citing_mosdef_gomc_python.html)
