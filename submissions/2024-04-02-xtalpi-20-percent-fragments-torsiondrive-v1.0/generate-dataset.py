#!/usr/bin/env python
# coding: utf-8

# # Generation process
# 
# This notebook documents the generation of a dataset of representative torsiondrive molecules shared by XtalPi. The conformers used are the post-optimization conformers shared by XtalPi.

# ## Imports

# In[1]:


import openff.qcsubmit
import openff.toolkit
import openeye
import qcelemental
import qcportal
import pyarrow
import pyarrow.dataset as ds
import numpy as np

print("OpenFF QCSubmit:", openff.qcsubmit.__version__)
print("OpenFF Toolkit:", openff.toolkit.__version__)
print("OpenEye:", openeye.__version__)
print("QCElemental:", qcelemental.__version__)
print("QCPortal:", qcportal.__version__)
print("PyArrow:", pyarrow.__version__)


# In[2]:


import tqdm

from openff.units import unit

from openff.toolkit import Molecule
from openff.toolkit.utils import OpenEyeToolkitWrapper, ToolkitRegistry

from openff.qcsubmit.datasets import TorsiondriveDataset
from openff.qcsubmit.factories import TorsiondriveDatasetFactory


# ## Setting up dataset

# In[3]:


dataset_factory = TorsiondriveDatasetFactory()


# ## Loading input

# In[4]:


input_dataset = ds.dataset("/data/chodera/lilywang/datasets/xff/output/xff-20-percent-td-dataset")
input_dataset.schema


# In[6]:


print(f"Total #confs: {input_dataset.count_rows()}")


# In[7]:


from openff.toolkit import Molecule
from openff.units import unit
import numpy as np

from openff.qcsubmit.workflow_components import TorsionIndexer

columns = [
    "mapped_smiles", "ATOM", "ENERGY",
    "frozen_atoms", "frozen_0", "frozen_1", "frozen_2", "frozen_3",
]
df = input_dataset.to_table(columns=columns).to_pandas()
all_molecules = []
# sort by torsion scan and take lowest energy conformer
groupby = ["mapped_smiles", "frozen_0", "frozen_1", "frozen_2", "frozen_3"]
for (mapped_smiles, *_), subdf in df.groupby(by=groupby):
    mol = Molecule.from_mapped_smiles(
        mapped_smiles,
        allow_undefined_stereo=True,
    )
    subdf = subdf.sort_values("ENERGY")
    conformer = np.array(subdf.ATOM.values[0]).reshape((-1, 3))
    mol._conformers = [conformer * unit.angstrom]
    
    torsion_indexer = TorsionIndexer()
    # frozen atoms indexes from 1
    frozen_atoms = tuple(
        [atom - 1 for atom in subdf.frozen_atoms.values[0]]
    )
    assert max(frozen_atoms) < len(mol.atoms)
    torsion_indexer.add_torsion(frozen_atoms, (0, 3), (-165, 180))
    mol.properties["dihedrals"] = torsion_indexer
    
    all_molecules.append(mol)

len(all_molecules)


# In[8]:


dataset = dataset_factory.create_dataset(
    dataset_name="XtalPi 20-percent Fragments TorsiondriveDataset v1.0",
    tagline="B3LYP-D3BJ/DZVP torsion drive of 20% the fragment dataset used by XtalPi to fit XFF.",
    description=(
        "A dataset containing 20% the fragments used by XtalPi "
        "in fitting the XFF force field "
        "(DOI: 10.1021/acs.jctc.3c00920). "
        "Conformers are the post-optimization geometries shared by XtalPi. "
        "Each conformer will be converged according to the 'GAU_LOOSE' criteria."
    ),
    molecules=all_molecules
)
dataset.metadata.submitter = "lilyminium"
dataset.metadata.long_description_url = (
    "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
    "submissions/"
    "2024-04-02-xtalpi-20-percent-fragments-torsiondrive-v1.0"
)


# In[9]:


dataset.dict()


# ## Exporting dataset

# In[10]:


dataset.export_dataset("dataset.json.bz2")
dataset.molecules_to_file('dataset.smi', 'smi')
dataset.visualize("dataset.pdf", columns=8)

print(dataset.qc_specifications)


# ## Dataset information

# In[11]:


import numpy as np
from collections import Counter


# In[12]:


print("n_molecules:", dataset.n_molecules)
print("n_conformers:", dataset.n_records)


# In[13]:


n_confs = np.array(
    [mol.n_conformers for mol in dataset.molecules]
)
n_heavy_atoms = np.array(
    [mol.to_rdkit().GetNumHeavyAtoms() for mol in dataset.molecules]
)


# In[14]:


print(
    "Number of conformers (min, mean, max):",
    n_confs.min(), n_confs.mean(), n_confs.max()
)
print("# heavy atoms")
counts = Counter(n_heavy_atoms)
for n_heavy in sorted(counts):
    print(f"{str(n_heavy):>3}: {counts[n_heavy]}")


# In[15]:


unique_charges = set([
    mol.total_charge.m_as(unit.elementary_charge)
    for mol in dataset.molecules
])
unique_charges


# In[16]:


masses = np.array([
    sum([atom.mass.m for atom in mol.atoms])
    for mol in dataset.molecules
])
print("MW (min, mean, max):", masses.min(), masses.mean(), masses.max())


# In[17]:


elements = set(
    atom.symbol
    for mol in dataset.molecules
    for atom in mol.atoms
)
print(elements)


# In[ ]:




