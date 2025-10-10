# Frequently Asked Questions (FAQ)

## Questions

<details>
<summary>What keywords are recommended for psi4 specifications?</summary>

### Recommended Psi4 Program Specification Keywords

#### SCF Multipole Reference Point
- **QCPortal Input** `QCSpecification(keywords={"properties_origin": ["NUCLEAR_CHARGE"]})`
- **Purpose**: Ensures that SCF dipole, quadrupole, and higher-order moments are computed relative to the nuclear charge origin rather than the coordinate origin.
- **Warning**: Datasets that did not include this keyword can only reliably use the dipole moment of neutral molecules. Higher-order moments and multipole moments of charged molecules should be estimated from the partial charges.

#### Multipole Moments
- **QCPortal Input** `QCSpecification(keywords={"scf_properties": ['dipole', 'quadrupole']})`
- **QCSubmit Input** `QCSpec(scf_properties=['dipole', 'quadrupole'])`
- **Purpose**: Compute the SCF dipole, quadrupole, and possibly higher-order moments using the reference point defined by `"properties_origin"`, otherwise the coordinate origin is used.

#### Partial Charges
- **QCPortal Input** `QCSpecification(keywords={"scf_properties": ['lowdin_charges', 'mulliken_charges', 'mbis_charges']})`
- **QCSubmit Input** `QCSpec(scf_properties=['lowdin_charges', 'mbis_charges'])`
- **Purpose**: Compute atomic partial charges with a variety of methods. Lowdin is a popular option that seems to have the best balance of element coverage and passing the common sense test. MBIS charges are a popular option, but have a limited number of supported elements. Mulliken charges are often produced by default in other QC packages, but not generally recommended.
- **Warning**: MBIS charges will error for unsupported elements (e.g., I)

#### Atomic Spin Population
- **QCPortal Input** `QCSpecification(keywords={"scf_properties": ['lowdin_spin']})`
- **Purpose**: Compute the per atom fractional number of unpaired electrons. They can be used as is or used to estimate the spin angular momentum as \( S_i(S_i + 1)\hbar^2 \).

#### Bond Indices
- **QCPortal Input** `QCSpecification(keywords={"scf_properties": ['wiberg_lowdin_indices', 'mayer_indices']})`
- **QCSubmit Input** `QCSpec(scf_properties=['wiberg_lowdin_indices', 'mayer_indices'])`
- **Purpose**: These indices provide a measure of bond order derived from the density matrix and are useful for analyzing the bonding characteristics in a molecule. The Wiberg bond indices calculated using the Löwdin orthogonalized atomic orbitals. The Mayer bond indices are commonly used to quantify the strength of bonding interactions between atoms in a molecule.

#### Dipole Polarizabilities
- **QCPortal Input** `QCSpecification(keywords={'function_kwargs': {'properties': ['dipole_polarizabilities']}})`
- **QCSubmit Input** `QCSpec(scf_properties=['dipole_polarizabilities'])`

#### Other Keywords
See the [psi4 keyword documentation](https://psicode.org/psi4manual/master/oeprop.html#basic-keywords) for other options.

#### Example:
```python
spec = QCSpecification(
        program='psi4',
        driver=SinglepointDriver.gradient,
        method='b3lyp-d3bj',
        basis='dzvp',
        keywords={
            'maxiter': 500, 
            'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges', 'mulliken_charges', 'mbis_charges'],
            'function_kwargs': {'properties': ['dipole_polarizabilities']},
            'properties_origin': ['NUCLEAR_CHARGE']
        },
        protocols={'wavefunction': 'none'}
    )
```

- **Reference**: [Psi4 Documentation](https://psicode.org/psi4manual/master/oeprop.html#basic-keywords)

</details>

<details>
<summary>How to define a torsion drive?</summary>

### Defining Torsion Drive Constraints

Torsion drive constraints are used to specify the dihedral angles that should be scanned during a torsion drive. These constraints are defined in the dataset submission file.

**WARNING:** `Molecule.connectivity` must be defined so that dihedrals can be verified.

Referencing datasets in the `submissions` folder, such as `2025-04-10-OpenFF-Additional-Generated-ChEMBL-TorsionDrives-4.0`. The specification program would be `"torsiondrive"`

#### In QCPortal:
```python
TorsiondriveDatasetEntry(
    initial_molecules=[offmol],
    additional_keywords={
        "dihedrals": [[ 3, 0, 1, 2]],
        "grid_spacing": [15],
        "dihedral_ranges": [[ -165, 180]],
        "energy_upper_limit": 0.05,
})
```

#### In QCSubmit:
```python
from openff.qcsubmit.utils import get_symmetry_classes, get_symmetry_group
from openff.qcsubmit.workflow_components import TorsionIndexer

torsion_indexer = TorsionIndexer()
symmetry_classes = get_symmetry_classes(offmol)

atom_indices = (i, j, k, l)
central_bond = tuple(sorted(atom_indices[1:-1]))
symmetry_group = get_symmetry_group(central_bond, symmetry_classes)
torsion_indexer.add_torsion((i, j, k, l), symmetry_group, (-165, 180))

offmol.properties["dihedrals"] = torsion_indexer
```

- **Note**: Ensure that the indices correspond to the correct atoms in the molecule, the connectivity of the molecule must be present for validation.

</details>

<details>
<summary>How to implement constraints (in geomeTRIC)?</summary>

### Implementing Constraints with Program = geomeTRIC


Most if not all datasets in this repository that involved geometry optimization did so with geomeTRIC with psi4 as an objective function. If one were to freeze a bond, angle, or dihedral during a geometry optimization, they would have to do so here.

#### Example Dataset Reference
For a practical example of implementing constraints in geomeTRIC, refer to the dataset submission `2025-03-05-OpenFF-Protein-PDB-4mer-v4.0`

#### In QCPortal:
```python
OptimizationEntry(
    initial_molecule=offmol,
    additional_keywords={"constraints": {"freeze": [{ "type": "dihedral", "indices": [ 1, 3, 4, 5]},]}
)
```

#### In QCSubmit:
```python
offmol.add_constraint(
    constraint = 'freeze', 
    constraint_type = 'dihedral', 
    indices = constraint_index,
    bonded=True
)
```

- **Reference**: [geomeTRIC Documentation](https://geometric.readthedocs.io/en/latest/constraints.html)

</details>

<details>
<summary>Why do some directories have scaffold.json files instead of dataset.json files?</summary>

There are some molecular systems, e.g., transition metal complexes, that are not handled well by neither OpenEye nor RDKit (the toolkits leveraged in QCSubmit). The scaffold.json files are a serialization of QCPortal datasets created in a FractalSnowflake server to bypass QCSubmit, but the CI can still validate the dataset contents and submit the dataset for you. If your dataset can be submitted with QCSubmit, it MUST be prepared with QCSubmit, otherwise we will ask to update your PR.

</details>

<details>
<summary>Tagging records for compute, and tagging by molecular weight</summary>

Records can be assigned a specific compute tag with our GitHub Action CI. In the PR used to create the dataset, add a label `compute-<my tag>`.  
The `<my tag>` portion will be the updated compute tag to be used to find the records on QCArchive, such as in NRP. Commonly the PR number is chosen.

To more efficiently use compute, you can have the CI separate the molecules into molecular weight (MW) bins and tag accordingly.  

A tag like, `compute-<my tag>_200-400-600` will group the molecules into groups where:  
- `<my tag>-200` has a MW of 200 Da or less,  
- `<my tag>-400` is between 200 Da and 400 Da,  
- `<my tag>-600` is similarly between 400 Da and 600 Da, and  
- `<my tag>-large` is > 600 Da.  

Any number of sequential bin boundaries can be strung together with hyphens.  

Once the GitHub tag is in place, you must run the GitHub Action "Dataset Lifecycle - Reprioritize/Retag" to propagate the GitHub labels to QCArchive.

</details>