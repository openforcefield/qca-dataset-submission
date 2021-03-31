# OpenFF QCArchive Dataset Submission

## Dataset Lifecycle

All datasets submitted to QCArchive via this repository conform to the [Dataset Lifecycle](./LIFECYCLE.md).

See [STANDARDS.md](./STANDARDS.md) for submission standards.
Datasets must be submitted as pull requests.


## Dataset Priority

This section contains a relative ranking of dataset priority.

While queueing is still handled manually, *this* is the single place where priority is determined.

New datasets that are added should be included in this list if they are to be submitted.

Do not modify this directly in master -- Always open a PR and get at least one approving review before merging.

* ~Roche hessians~
* ~Coverage hessians~
* ~Chaya's Biphenyl set~
* ~Boron~
* 06-28-19-Nitrogen-grid-optimization
* Everything else

Look into Roche torsion set to investigate/fix failures.

# qca-dataset-submission
Data generation and submission scripts for the QCArchive ecosystem.

## `TorsionDriveDataset`
* `2019-06-03-Torsion-Stability-Benchmark`: source files for `Fragment Stability Benchmark` (`TorsionDriveDataset`) covering molecular fragments of varying size around a central torsion (@ChayaSt)
* `2019-07-25-phenyl-set`: source files for `phenyl set`
* `2019-09-07-Pfizer-discrepancy-torsion-dataset-1`: source files for `Pfizer discrepancy torsion dataset 1`

## `OptimizationDataset`
* `2019-06-25-smirnoff99Frost-coverage`: source files for `smirnoff99Frost parameter coverage` (`OptimizationDataset`) providing minimal coverage of SMIRNOFF valence parameters (@ChayaSt)
* `2019-07-02 VEHICLe optimization dataset`: source files for `Open Force Field VEHICLe optimization dataset 1.0.0` (`OptimizationDataset`) for the VEHICLe dataset (heteraromatic rings of the future) (@jchodera)
* `2019-07-05 OpenFF NCI250K Boron 1`: source files for `OpenMM NCI250K Boron 1` (`OptimizationDataset`) where small boron-containing compounds are extracted from the [NCI250K](https://cactus.nci.nih.gov/download/nci/) (@jchodera)
* `2019-09-07-Pfizer-discrepancy-optimization-dataset-1`: source files for `Pfizer discrepancy optimization dataset 1` (@jchodera)
