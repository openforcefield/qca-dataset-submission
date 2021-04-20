
### Short description
Step3. Filter the molecule set from the previous step, to eliminate internal hydrogen bond forming molecules.

### Manifest
- `molecules_out.smi`: copied from `../2.gen-mol-list`
- `filter-baker-hubbard.py`: Written by SB. Generates conformers for each molecule and find molecules having internal hydrogen bond forming conformers.  
- `smiles-to-keep.smi`: Molecules not forming internal hydrogen bonds.
- `smiles-to-remove.smi`: Molecules forming internal hydrogen bonds.
- `smiles-with-errors.smi`: Molecules failed to be parsed into a molecule object.