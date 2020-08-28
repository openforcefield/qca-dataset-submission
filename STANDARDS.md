# Dataset Standards

This file outlines the standards and requirements needed for submitting a dataset to QCArchive.
This ensures that we have a consistent data model for downstream processes.

STANDARDS version: 3-DRAFT (2020-08-23)

Implementation:

	- QCEngine 
	- QCElemental
	- QCFractal
	- QCSubmit

# Required fields 

Training a new Open Force Field requires a specific set of data. This data touches on multiple aspects, including the underlying the theoretical chemistry models, the data structures that organize the information, and the data required for downstream calculations.

## Required fields layout for Input

This layout closely follows the API and implementation provided by QCArchive. All objects can be represented as JSON, which will be provided as examples below.

### All datasets
	- Name: str
	- Metadata: Dict
		- A human description of the entry names: str
		- The changelog: Dict
		- A short description: str
		- A long description: str
		* elements: List
	* The Entry Specification:
		- name: "default"
		- description: str
		- qc_spec
			- driver: str
			- method: "b3lyp-d3bj"
			- basis: "dzvp"
			- keywords: Dict
				- maxiter: 200
				- scf_properties: List
					- dipole
					- quadrupole
					- wiberg_lowden_indices
					- mayer_indices
	* Entries: Dict
		- attributes: Dict
			- `canonical_isomeric_explicit_hydrogen_mapped_cmiles`
		- The molecule configurations: Dict
			- A valid QCSchema molecule with the follwing also required:
				- `extras`
					- `canonical_isomeric_explicit_hydrogen_mapped_cmiles`
				- `connectivity`

- Name: str
	The name of the dataset. This is the name that will be publicly visible and needed for query/access.
- Metadata: Dict
	- entry_names_description: str
	 This field helps clarify the convention used for naming the entries. Each entry is named, and this name is an arbitrary string. This is a description to help understand any information that the names are conveying. See section `Best Practices` for examples. 
	- changelog: Dict(version: str, description: str)
	A description of the changes made to the dataset for each version.
	- short_description: str
	A concise description of the dataset and the information it contains. 
	* long_description: str
	A longer description which contains all information relevent to the dataset
	* elements: List
	A list of elements that are included in the dataset molecules.

- The Entry Specification
Every dataset should have a specification named `'default'`, which includes all of the standard settings Open Force Field uses for its calculations. The complete requirements of this object are specific to each type of dataset; see the proceeding dataset sections for details. The following fields are required by all datasets.
	- name: str
	The name of the specification. A specification is analogous to the details and settings used for the calculations: all calculations grouped by a specification only differ by the input molecule.
	- description: str
	The description of this specification. See best practices for suggestions on what to put in this description
	- The Quantum Chemistry specification: QCSpecification
		- driver: str
		The type of calculation. This setting depends on the type of dataset.
		- method: str
		The method used for calculation. The `default` specification requires this to be `B3LYP-D3BJ`
		* basis: str
		The basis set used for calculation. The `default` specification requires this to be `dzvp`, and is a name specific to PSI4.
		* keywords: KeywordSet
		A collection of options used by the program. For `default` this is the settings for PSI4. Each keyword set is defined by an ID, and on the public database this is current ID '2':
			- maxiter: int
			The maximum number of SCF iterations to use in PSI4. The `default` spec requires this to be 200.
			* scf_properties: List
			A list of properties that PSI4 should report. The `default` spec requires this to be:
				- dipole
				- quadropole
				- wiberg_lowdin_indices
				- mayer_indices

- The Entries
Each entry represents a molecule. It includes the "chemical" information, such as SMILES and other metadata, as well as the "physical" information, such as positions. The information in an entry is dataset specific. Below represents the information needed by all datasets.
	* attributes: Dict
		- canonical_isomeric_explicit_hydrogen_mapped_cmiles: str
		The canonical isomeric explicity hydrogen mapped SMILES string of the molecule. The map connects the absolution position of each atom in the SMILES string to the positions defined in the molecule entries. Thus the string `[H:2][Cl:1]` would indicate that the first position in the conformations would correspond to `Cl`. The indices in the map always start at 1; maps pointing to 0 generally mean "unmapped".

### Optimizations

The optimization dataset requires the following additional details.

- The `default` EntrySpecification:

	- The Optimization Specification: OptimizationSpecification
		- program: "geometric"
		- keywords: Dict
			- coordsys: "tric"
			- enforce: 0.1
			- reset: True
			- qcconv: True
			- epsilon: 0.0

* The Optimization Specification
	The program to use for the optimization and the program-specific options it uses. The `default` spec requires the program to be `geometric`. The keywords vary depending on the dataset, for OptimizationDatasets they as listed above.

	- Conforming example:

	TODO: ADD HERE
	
- TorsionDrives

	Essentially the same as OptimizationDatasets, but adding to each Entry requires a list of torsions (4-tuples) to drive, and the constant interval spacing. 
	Additionally, each Entry contains a list of input molecules, whereas OptimizationDataset Entries contain a single molecule conformation.

	After calculation, the Optimizations will have the following information

	- Conforming example:
	TODO: ADD HERE

- Hessians
	- The qc_spec options
	- Conforming example:
	TODO: ADD HERE

- GridOptimizations
	- The qc_spec options
	- Conforming example:
	TODO: ADD HERE

* Training sets
	- Metadata
		- DOI 
	- Necessary contributed values

	- Conforming example:
	TODO: ADD HERE

* Benchmarking sets
	- Metadata
		- DOI
	- Necessary contributed values

	Conforming example:
	TODO: ADD HERE

## Required fields layout after calculation

This specifies the fields that will be present after a calculation has been performed.

## Job specifications (level of theory; settings)

OpenFF depends on a QCSpecification named "default" which corresponds to `B3LYP-D3BJ/DZVP` in Psi4. Submissions may have multiple specifications, but must include the `default`. 

# Best practices

* If any calculations are to be redone from another collection, re-use the old input (coordinates, atom ordering etc) as this will avoid running the calculation again and will just create new references in the database to the old results and should help keep the cost of the calculations down.  

# Entry naming

Entries names are not subject to any requirements other than being a string of characters. Appropriate names can give more depth and insight to a dataset, so choosing an effective naming scheme should be taken seriously. A typical example is torsion drives, where the entry names are SMILES string with all elements converted to lower case, and the 4 atoms involved in the torsion are labeled, e.g. `'c[c:1][c:2][c:3][c:4]c'`. Another example is in optimization sets, where the lower-cased SMILES with `-N` appended is used, which signifies the conformation for that given molecule. Another option is just to name each entry with a stringified number, e.g. `'1', '2', '3'`, etc, but this does not convey any information, and should be avoided if possible. One case where this may potentially be appropriate is if the structures are large (proteins) and a short, managable name is not possible. However, one may decide to use sequences (`MEGFFFKSS`) or PDB IDs instead. Overall, the hallmarks of a good name is utility and readability. Whichever naming scheme is used, the `entry_names_description` must describe the convention.

# Dataset naming and versioning

Each dataset shall be versioned.
- The naming of a dataset should have the following structure:

    `"OpenFF <descriptive and uniquely-identifying name> v<version number>"`

- The first submission of a dataset will have a version `"v3.0"`

* The major version shall indicate the STANDARDS that the dataset conforms to. Datasets which are not intended to conform to any STANDARDS should start with 0, e.g. `"v0.1"`. 

* Datasets with versions starting with `"v1.x"` and `"v2.x"` do not follow any official STANDARDS, and thus should be considered `"v0.x"`.

- A dataset with the suffix `"-beta"` is not to be used for production work.

- A minor version change (e.g. `"v3.1"`) represents a cosmetic or minor additions/problems that were addressed:
    - Cosmetic changes
	- Errors/bugs in the molecule specification
	- Changes necessary to adhere to the STANDARDS

A best-effort is made to ensure that a dataset follows its underlying STANDARDS. One must assume that the newest version of a dataset best conforms to these STANDARDS, and the same promise may not hold for earlier versions. The changelog should address any changes made to improve compliance.

Each version increment should take the information from the previous `changelog` field, and add a new entry of the form { "version": description } that explains the modifications made to the dataset. Each dataset should therefore have the complete changelog.

## Tags indicate status

A tag `"complete"` indicates that a dataset is completed as far as OpenFF is concerned.
This means that any errors remaining are known to be acceptable or impossible to fix.
It also means that no additional work is being done on the dataset to get it to completion.

A tag `"inflight"` indicates that a dataset is not completed as far as OpenFF is concerned.
This means that any errors remaining are being actively addressed.

All datasets should also feature a `"openff"` tag.

## Force Field Releases

When a new force field is released, a dataset corresponding to all results used for the force field fitting should be created.
This gives a single reference for these data instead of many references.
The format of these dataset names is:

    `"OpenFF Force Field <friendly name> <ff version>"`

## Group

The dataset's group should be set to `"OpenFF"`.

## Molecule validation

* See ["Molecule submission checklist"](https://github.com/openforcefield/qcsubmit/issues/9)

* Unique keys in `ds.data.records` must not reference the same entry.
* Unless explicitly specified in the submission descriptions, torsion drives must be on 4 connected atoms
* Torsions driving a ring will give a warning, and torsions in a a ring of  3,4,5,6 atoms is considered an error
* Warnings will be given if an atom does not have a complete valence set

# Standard functions and modules for entry preparation

* QCSubmit (https://github.com/openforcefield/qcsubmit)

## Related/ongoing discussions

### Required fields

* See ["Fields that should be required for OpenFF submissions"](https://github.com/openforcefield/qcsubmit/issues/3)


### Adding Compute to old datasets
Datasets now support multiple QC specifications and will start compute for them all simultaneously when submitted.
However in some cases you may want to add new specifications to old datasets already in the archive, to do this make a PR in the normal way with either a `dataset.json` or `compute.json` qcsubmit dataset. 
The dataset should be of the correct type and have the name set to that in the archive.   The dataset entries should be empty and only the new `qc_specifications` section should be filled in which will cause the 
CI to search the public archive for the dataset and validate the basis coverage before submitting. 
