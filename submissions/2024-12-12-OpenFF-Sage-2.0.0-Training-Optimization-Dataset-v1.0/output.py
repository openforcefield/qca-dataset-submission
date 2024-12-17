Number of results:  3663
Finished converting to records
Molecules are the same?  True
# heavy atoms
  4: 2
  5: 3
  6: 4
  7: 7
  8: 20
  9: 19
 10: 46
 11: 56
 12: 81
 13: 81
 14: 94
 15: 94
 16: 59
 17: 59
 18: 47
 19: 38
 20: 36
 21: 45
 22: 30
 23: 23
 24: 21
 25: 29
 26: 6
 27: 11
 28: 16
 29: 10
 30: 8
 31: 16
 32: 16
 33: 11
 34: 11
 35: 5
 36: 15
 37: 7
 38: 8
 39: 5
* Description: A quantum chemical (QC) dataset curated to train OpenFF 2.0.0 Sage, with reparametrized Lennard-Jones (LJ) and valence parameters, the latter relevant to this dataset. This QC dataset with the OpenFF default level of theory, B3LYP-D3BJ/DZVP, is used to benchmark Sage geometries and energetics. These optimized conformer geometries where used in conjunction with the QC dataset used to train one dimensional torsional profiles. This Generation 2 dataset increases chemical diversity when compared to Generation 1, which are of value to our industry partners. Large molecules (>20 heavy atoms) were also included, including more flexible molecules and a greater degree of conformational variation which provide intramolecular interactions.

This is the complete optimization dataset used for training OpenFF 2.0.0 Sage, consisting of the following datasets: 'OpenFF Gen 2 Opt Set 1 Roche', 'OpenFF Gen 2 Opt Set 2 Coverage', 'OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy', 'OpenFF Gen 2 Opt Set 4 eMolecules  - Discrepancy', and 'OpenFF Gen 2 Opt Set 5 Bayer'.
The following filters were applied: RecordStatusFilter(status=RecordStatusEnum.complete), ConnectivityFilter(tolerance=1.2), UndefinedStereoFilter(), ConformerRMSDFilter(max_conformers=10), and ElementFilter(allowed_elements=['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'])Further information can be found in the curation scripts for the linked repositories.
* Purpose: B3LYP-D3BJ/DZVP conformers applicable to drug-like molecules for OpenFF 2.0.0 Sage
* Name: OpenFF Sage 2.0.0 Training Optimization Dataset v1.0
* Number of unique molecules: 1025
* Number of filtered molecules: 0
* Number of conformers: 3663
* Number of conformers (min, mean, max): 1.00, 3.53, 10.00
* Molecular weight (min, mean, max): 76.05, 261.38, 544.64
* Charges: -2.0, -1.0, 0.0, 1.0
* Submitter: Jennifer A Clark
## Metadata
* Elements: {S, Br, O, P, H, C, N, I, Cl, F}
* QC Specifications: default
  * basis: DZVP
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 200
  * method: B3LYP-D3BJ
  * program: psi4
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
