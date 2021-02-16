This brings in examples provided by Alberto Gobbi from Genentech as also deposited in our [`open-forcefield-data`](https://github.com/openforcefield/open-forcefield-data/tree/master/pdb-examples) repo.

## Manifest
- `pubLigsChargedGoodDensity.sdf`: SDF file (see dataset description) with charged PDB ligands.
- `pubLigsNeutralGoodDensity.sdf`: SDF file (see dataset description below) with neutral PDB ligands.
- `pub*.pdf`: PDF visualizations of the same.
- `visualize.py`: Generates PDF of compounds


# Dataset description from originator

## Examples from the [PDB database](https://www.rcsb.org/)
Ligands were extracted using the in-house Proasis
software and filtered by multiple criteria:

- Electron density quality criteria such as Resolution and Rvalue and others
- MW <= 500, Rotatable Bonds <= 6, no none-organic atoms

### pubLigsNeutralGoodDensity.sdf.gz:
648 ligands from pdb. Starting with criteria mentioned above then further processed as follows:

- Computed pKa: most basic > 6, most acidic < 8 using [MoKa](http://www.moldiscovery.com/software/moka/)
- All ionized atoms where neutralized by adding or removing H when possible.
- Some cleanup of the connectivity was performed but there might still be issues.


### pubLigsChargedGoodDensity.sdf.gz
382 ligands from pdb. Starting with criteria mentioned above then further processed as follows:

- Computed pKa: most basic > 8, most acidic < 6 using [MoKa](http://www.moldiscovery.com/software/moka/)
- Ligands were also computationally protonated to their assumed state at pH 7.0 using [tauter](http://www.moldiscovery.com/software/moka/).

  Note: the computational method used is known to generate wrong protonation states and tautomers with low population.

- Some cleanup of the connectivity was performed but there might still be issues.



#### Tools used:

- [Chemalot](https://github.com/chemalot/chemalot)

    Lee, Man-Ling, Ignacio Aliagas, Jianwen A. Feng, Thomas Gabriel, T. J. O�Donnell, Benjamin D. Sellers, Bernd Wiswedel, and Alberto Gobbi.
    �Chemalot and Chemalot_knime: Command Line Programs as Workflow Tools for Drug Discovery.�
    Journal of Cheminformatics 9 (June 12, 2017): 38. https://doi.org/10.1186/s13321-017-0228-9.

- [Proasis](http://www.desertsci.com/)

- [MoKa](http://www.moldiscovery.com/software/moka/)

	New and Original pKa Prediction Method Using Grid Molecular Interaction Fields
	Francesca Milletti, Loriano Storchi, Gianluca Sforna, and Gabriele Cruciani
	J. Chem. Inf. Model., 2007, 47 (6), pp 2172-2181

	Tautomer Enumeration and Stability Prediction for Virtual Screening on Large Chemical Databases
	Francesca Milletti, Loriano Storchi, Gianluca Sforna, Simon Cross and Gabriele Cruciani
	J. Chem. Inf. Model., 2009, 49 (1), pp 68�75

	Tautomer Preference in PDB Complexes and its Impact on Structure-Based Drug Discovery
	Francesca Milletti and Anna Vulpetti
	J. Chem. Inf. Model., 2010, 50 (6), pp 1062�1074

- [OpenEye](https://www.eyesopen.com/)
