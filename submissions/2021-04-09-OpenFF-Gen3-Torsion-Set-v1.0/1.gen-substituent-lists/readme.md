
### Short description
Step1. Generation of a list of small and rigid substituents.

### Manifest
- `gen-substituent-lists.py`: (1) reads input molecule files from `input-mol-sets`; (2) generates a pickle file, storing a list of substituents and a pdf file for each molecule set; (3) combined four lists of substituents into one. 
- `utils.py`: contains functions used in `gen-substituent-lists.py` 
- `substituents_combined.p`: pickle file storing a combined list of substituents. 2D image of the substituents can be found in `substituents_combined.pdf`
- `substituents_filtered.p`: pickle file storing a combined list of substituents, with uninteresting substituents eliminated.  (1. Hydrazines, 2. Halogen attached to noncarbon, 3. Pentavalent nitrogens, 4. Fluoren isotope, 5. Different halogen atoms attached to the same carbon)