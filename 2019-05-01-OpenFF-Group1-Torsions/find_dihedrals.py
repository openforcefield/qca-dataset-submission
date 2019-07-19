#!/usr/bin/env python

import json
import collections
from itertools import combinations

from forcebalance.molecule import Molecule, Elements
from bond_graph import BondGraph

class DihedralSelector:
    def __init__(self, molecule, skip_straight=True):
        self.m = molecule
        self.bond_graph = BondGraph(self.m.bonds)
        self.avoid_angles_set = set()
        if skip_straight:
            self.avoid_angles_set = self.get_straight_angles()

    def get_straight_angles(self, threshold=165.0):
        straight_angles = set()
        self.m.build_topology(force_bonds=False)
        for a in self.m.find_angles():
            if self.m.measure_angles(*a)[0] >= threshold:
                straight_angles.add(a)
                straight_angles.add(a[::-1])
        return straight_angles

    def find_dihedrals(self, dihedral_filters=None):
        """ Find all dihedrals, then apply filters """
        # default empty filter
        if dihedral_filters is None: dihedral_filters = []
        # get the list of dihedrals
        dihedral_list = self.bond_graph.get_dihedrals()
        print(f'Found total {len(dihedral_list)} distinct dihedrals')
        if self.avoid_angles_set:
            dihedral_list = [d for d in dihedral_list if d[:3] not in self.avoid_angles_set and d[1:] not in self.avoid_angles_set]
            print(f'{len(dihedral_list)} dihedrals left after filter out straight angles')
        # apply filters
        available_filters = {
            'equiv_terminal': self.filter_equivalent_terminals,
            'heavy_atoms': self.filter_keep_4_heavy,
            'no_ring': self.filter_remove_ring,
            'unique_center_bond': self.filter_keep_unique_center
        }
        for filt_name in dihedral_filters:
            if filt_name not in available_filters:
                raise ValueError(f"Filter named {filt_name} not recognized, choices are {available_filters.keys()}")
            filt = available_filters[filt_name]
            dihedral_list = filt(dihedral_list)
        return dihedral_list

    def filter_equivalent_terminals(self, dihedral_list):
        eq_idxs = self.find_equivalent_terminal_atom_idxs()
        print(f"Filtering based on equivalent terminal atoms: {eq_idxs}")
        dihedrals_set = set()
        for i,j,k,l in dihedral_list:
            skipped = False
            for eq_i in eq_idxs[i]:
                for eq_l in eq_idxs[l]:
                    if (eq_i, j, k, eq_l) in dihedrals_set:
                        print(f"Filter: dihedral {i}-{j}-{k}-{l} skipped because equivalent to {eq_i}-{j}-{k}-{eq_l}")
                        skipped = True
                        break
                if skipped: break
            if not skipped:
                dihedrals_set.add((i,j,k,l))
        # resume the ordering of dihedrals
        dihedral_list = sorted(dihedrals_set, key=lambda d: (d[1], d[2], d[0], d[3]))
        return dihedral_list

    def find_equivalent_terminal_atom_idxs(self):
        elem_list = self.m.elem
        neighbor_list = self.bond_graph
        noa = self.m.na
        equal_atom_idxs = {i:{i} for i in range(noa)}
        for i in range(noa):
            for j in range(i+1, noa):
                if elem_list[i] == elem_list[j]:
                    if len(neighbor_list[i]) == len(neighbor_list[j]) == 1:
                        if neighbor_list[i] == neighbor_list[j]:
                            equal_atom_idxs[i].add(j)
                            equal_atom_idxs[j].add(i)
        return equal_atom_idxs

    def filter_keep_4_heavy(self, dihedral_list):
        """ Filter dihedrals, keep only with 4 heavy atoms """
        print("Filter: only keep dihedrals formed by 4 heavy atoms")
        elem_list = self.m.elem
        filtered_dihedral_list = [d for d in dihedral_list if not any(elem_list[idx] == 'H' for idx in d)]
        print(f"Number Left: {len(dihedral_list)} => {len(filtered_dihedral_list)}")
        return filtered_dihedral_list

    def filter_remove_ring(self, dihedral_list):
        """ Filter dihedrals, remove any dihedral that's inside a ring """
        rings = self.bond_graph.get_rings()
        print(f"Filter: Removing dihedrals that the center bond is in any of the rings {rings}")
        # build a dictionary stores the index of ring each atom belongs to
        d_rings = collections.defaultdict(set)
        for i_ring, ring in enumerate(rings):
            for atom_idx in ring:
                d_rings[atom_idx].add(i_ring)
        # go over dihedrals and check if any four belong to the same ring
        filtered_dihedral_list = []
        for i, j, k, l in dihedral_list:
            if d_rings[j] & d_rings[k]:
                print(f"Dihedral {i}-{j}-{k}-{l} skipped because bond {j}-{k} is in a ring")
            else:
                filtered_dihedral_list.append([i,j,k,l])
        print(f"Number Left: {len(dihedral_list)} => {len(filtered_dihedral_list)}")
        return filtered_dihedral_list

    def filter_keep_unique_center(self, dihedral_list):
        """ Keep only one dihedral for each unique center bond """
        print(f"Filter: Keep only one dihedral for each center bond")
        filtered_dihedral_list = []
        d_center_bond = collections.defaultdict(list)
        for i,j,k,l in dihedral_list:
            center_bond = (j,k) if j < k else (k,j)
            d_center_bond[center_bond].append([i,j,k,l])
        for center_bond, dihedral_candidates in d_center_bond.items():
            print(f"best dihedral among {dihedral_candidates}:")
            best_dihedral = self.find_best_dihedral_same_center_bond(dihedral_candidates)
            print(best_dihedral)
            filtered_dihedral_list.append(best_dihedral)
        print(f"Number Left: {len(dihedral_list)} => {len(filtered_dihedral_list)}")
        return filtered_dihedral_list

    def find_best_dihedral_same_center_bond(self, dihedral_candidates):
        """ Find the best dihedral among candidates with same center bond
        Definition of best dihedral i-j-k-l: (From Lee-Ping)
        Temporarily disconnect all i-j bonds, then check the total number of connected atoms for each i,
        Same method applies to all candidates of l.
        The dihedral angle with the maximum connected_i + connected_j wins
        Return a single dihedral as [i, j, k, l]
        """
        if len(dihedral_candidates) == 0: return
        # check center bond are all the same
        _, center_j, center_k, _ = next(iter(dihedral_candidates))
        assert all(j==center_j and k==center_k for i,j,k,l in dihedral_candidates), "all candidates should share same center"
        # build new bond graph with only heavy atoms
        heavy_atom_bonds = [[b1, b2] for b1, b2 in self.m.bonds if self.m.elem[b1] != 'H' and self.m.elem[b2] != 'H']
        # get a new bond graph with only heavy atoms
        bond_graph = BondGraph(heavy_atom_bonds)
        # find the best i among all candidates
        i_candidates = {i for i,_,_,_ in dihedral_candidates}
        # compute and store the number of connected atoms
        n_connected_i = {}
        if len(i_candidates) == 1:
            n_connected_i[i_candidates.pop()] = 0
        else:
            # temporarily remove all i-j bonds
            for i in i_candidates:
                bond_graph.remove_bond(i, center_j)
            # compare i_candidates and find the one with most connected atom
            for i in i_candidates:
                # get all atoms connect to i in the temporary graph
                n_connected_i[i] = len(bond_graph.get_connected_nodes(i))
            print(f"n_connected for each i: {n_connected_i}")
            # add back all i-j bonds
            for i in i_candidates:
                bond_graph.add_bond(i, center_j)
        # find the best_l among all candidates
        l_candidates = {l for _,_,_,l in dihedral_candidates}
        n_connected_l = {}
        if len(l_candidates) == 1:
            n_connected_l[l_candidates.pop()] = 0
        else:
            # temporarily remove all i-j bonds
            for l in l_candidates:
                bond_graph.remove_bond(center_k, l)
            # compare i_candidates and find the one with most connected atom
            for l in l_candidates:
                # get all atoms connect to i in the temporary graph
                n_connected_l[l] = len(bond_graph.get_connected_nodes(l))
            print(f"n_connected for each l: {n_connected_l}")
        # get the best dihedral
        best_dihedral = max(dihedral_candidates, key=lambda d: n_connected_i[d[0]] + n_connected_l[d[3]])
        return best_dihedral

    def write_dihedrals(self, dihedral_list, filename):
        with open(filename, 'w') as outfile:
            json.dump(dihedral_list, outfile, indent=2)

    def find_dihedral_pairs(self, pattern='ring-a-ring'):
        res = []
        if pattern == 'ring-a-ring':
            # Find all pairs of dihedrals that two center bonds are
            # (ring-a)-ring and ring-(a-ring)
            rings = self.bond_graph.get_rings()
            if len(rings) >= 2:
                rsets = [set(ring) for ring in rings]
                for r1, r2 in combinations(rsets, 2):
                    # find all paths between the two rings
                    all_paths = self.bond_graph.find_all_paths(r1, r2)
                    # we want only one path, and the path has len == 3 (exactly one bridge atom)
                    if len(all_paths) != 1 or len(all_paths[0]) != 3: continue
                    left, center, right = all_paths[0]
                    print(f"Found ring-a-ring: {r1}-{center}-{r2}")
                    # since both ends are rings, the (effective size)* of end groups are the same, randomly pick one
                    left_neighbor = (self.bond_graph[left] & r1).pop()
                    right_neighbor = (self.bond_graph[right] & r2).pop()
                    # add dihedral pair to result list
                    dihedral_pair = ((left_neighbor, left, center, right), (left, center, right, right_neighbor))
                    # filter out straight angles
                    if self.avoid_angles_set:
                        avoid_angles = {(left_neighbor, left, center), (left, center, right), (center, right, right_neighbor)} & self.avoid_angles_set
                        if avoid_angles:
                            print(f"{dihedral_pair} ignored because angle {avoid_angles} should be avoided")
                            continue
                    print(f"{dihedral_pair} added")
                    res.append(dihedral_pair)
        res.sort()
        return res



def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help='Input mol2 file for a single molecule')
    parser.add_argument("-o", "--outfile", default='dihedral_list.json', help='Output json file containing definition of dihedrals')
    args = parser.parse_args()

    molecule = Molecule(args.infile)
    selector = DihedralSelector(molecule)
    dihedral_list = selector.find_dihedrals()
    print(f"Found {len(dihedral_list)} dihedrals")
    for d in dihedral_list:
        print(d)
    selector.write_dihedrals(dihedral_list, args.outfile)

if __name__ == '__main__':
    main()
