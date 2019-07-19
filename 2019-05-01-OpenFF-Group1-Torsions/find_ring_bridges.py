#!/usr/bin/env python

import os
from itertools import combinations
from forcebalance.molecule import Molecule
from bond_graph import BondGraph

mol_folder = 'processed_molecules/mol2'
total_count = 0
mol_with_bridges = []
for f in sorted(os.listdir(mol_folder)):
    fn = os.path.join(mol_folder, f)
    m = Molecule(fn)
    bg = BondGraph(m.bonds)
    rings = bg.get_rings()
    if len(rings) >= 2:
        rsets = [set(ring) for ring in rings]
        for r1, r2 in combinations(rsets, 2):
            # find all paths between the two rings
            all_paths = bg.find_all_paths(r1, r2)
            # we want only one path, and the path has len == 3 (exactly one bridge atom)
            if len(all_paths) == 1 and len(all_paths[0]) == 3:
                #print(f'found ring-bridge molecule {f}')
                mol_with_bridges.append(f)
                break
    total_count += 1


for f in mol_with_bridges:
    print(f)

print(f"Among {total_count}, found {len(mol_with_bridges)} molecules with bridge between rings")



