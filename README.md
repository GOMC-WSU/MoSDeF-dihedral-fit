# MoSDeF-dihedral-fit
MoSDeF-dihedral-fit: A simple software package to fit dihedrals via the MoSDeF software

This package is a tool to add dihedral parameters for [Ryckaert-Bellemans](https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-ryckaert-bellemans-function) dihedral based forcefields by fitting QM calculations for the rotation about the dihedral angle to this functional form. The main functionality is found through the following code snippet.
```python
>>> import unyt as u
>>> from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import fit_dihedral_with_gomc
>>> fit_dihedral_with_gomc(
        ["HC", "CT", "CT", "HC"],
        "molecule.mol2"
        "compound.xml",
        298.15 * u.Kelvin,
        gomc_binary_directory,
        {"HC_CT_CT_HC_multiplicity_1.log": []},
        zeroed_dihedral_atom_types=None,
        qm_engine="gaussian",
        VDWGeometricSigma=False,
        atom_type_naming_style="general",
        gomc_cpu_cores=1,
        fit_min_validated_r_squared=0.99,
        fit_validation_r_squared_rtol=1e-03,
    )
>>> import os
>>> os.system("cat RB_torsion_k_constants_fit_energy.txt")
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
Some basic workflows that use this package
    - [ethane dihedral](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDeF-dihedral-fit/ethane_HC_CT_CT_HC)
    - [propanoic acid dihedral](https://github.com/GOMC-WSU/GOMC_Examples/tree/main/MoSDeF-dihedral-fit/protonated_fragment_CT_CT_C_OH)
