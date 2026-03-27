from openff.qcsubmit.datasets import OptimizationDataset
from openff.qcsubmit.common_structures import Metadata, QCSpec, MoleculeAttributes
from collections import defaultdict
import logging
# # suppress stereochemistry warnings
logging.getLogger("openff").setLevel(logging.ERROR)


dataset = OptimizationDataset.parse_file('dataset.json.bz2')

dataset.metadata.long_description_url = (
        "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
        "submissions/"
        "2024-12-09-OpenFF-NAGL2-Training-Optimization-Dataset-v4.0"
    )


dataset.export_dataset("dataset.json.bz2")
dataset.molecules_to_file('dataset.smi', 'smi')
dataset.visualize("dataset.pdf", columns=8)
