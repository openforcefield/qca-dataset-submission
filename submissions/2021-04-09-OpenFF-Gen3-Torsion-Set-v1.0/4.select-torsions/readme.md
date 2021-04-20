
### Short description
Step4. List molecules matching to each torsion parameter and select torsions 

### Manifest
- `smiles-to-keep.smi`: copied from `../3.filter-mol-list`
- `gen_tid_mol_list.py`: lists molecules matching to each torsions and generate a dictionary    `tid_molecules_list[tid] = [{'mol_index': mol_index, 'indices': indices, 'covered_tids':covered_tids}, ...]`
- `tid_molecules_list.p`: pickle file storing `tid_molecules_list` from running `gen_tid_mol_list.py`
- `clustering.ipynb`: (1) clusters each list using DBSCAN and MACCS keys; (2)clusters each list into around 10 clusters; (3) chooses center molecule from each cluster; (4) save selected torsions into a dictionary, `tid_clusters_list[tid] = [ {'cluster_label': N, 'torsions': [...]}]`
- `tid_clusters_list.p`: pickle file storing `tid_clusters_list`. 2D images of selected torsions can be found in `elected_torsions.pdf`
