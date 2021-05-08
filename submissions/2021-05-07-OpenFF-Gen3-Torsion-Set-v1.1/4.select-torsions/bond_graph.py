#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import collections
import numpy as np

class BondGraph(object):
    def __init__(self, bonds):
        """ Initialize with a list of bonds: [(1,2), (3,4)] """
        self._bonds = []
        for b1, b2 in bonds:
            bond = [b1, b2] if b1 < b2 else [b2, b1]
            self._bonds.append(bond)
        self._data = collections.defaultdict(set)
        for b1, b2 in bonds:
            self._data[b1].add(b2)
            self._data[b2].add(b1)

    def __getitem__(self, key):
        return self._data[key]

    def copy(self):
        return BondGraph(self._bonds)

    def add_bond(self, node1, node2):
        self._data[node1].add(node2)
        self._data[node2].add(node1)

    def remove_bond(self, node1, node2):
        assert (node1 in self._data) and (node2 in self._data), "Both nodes need to be in this bond graph!"
        self._data[node1].remove(node2)
        self._data[node2].remove(node1)

    def remove_node(self, node1):
        """ Remove node and all its bonds """
        self._data.pop(node1)
        for node2 in self._data:
            self._data[node2].discard(node1)

    def find_path(self, node1, node2):
        """ Find the shortest path from node1 to node2 """
        assert (node1 in self._data) and (node2 in self._data), "Both nodes need to be in this bond graph!"
        explored = set()
        frontier = [ [node1] ]
        while frontier:
            path = frontier.pop(0)
            last = path[-1]
            for new in self._data[last]:
                if new == node2:
                    return path + [new]
                elif new not in explored:
                    explored.add(new)
                    path2 = path + [new]
                    frontier.append(path2)
        return None # path not found

    def find_all_paths(self, node_list1, node_list2):
        """ Find all possible paths from any node in node_list1 to node_list2
        The resulting paths should not have duplicate nodes.
        Also, the resulting paths should not have any node from node_list1 in the middle
        """
        assert all(node in self._data for node in node_list1), 'All nodes in node_list1 should be in this bond graph'
        assert all(node in self._data for node in node_list2), 'All nodes in node_list2 should be in this bond graph'
        set1 = set(node_list1)
        set2 = set(node_list2)
        frontier = [[node1] for node1 in node_list1]
        res_paths = []
        # the overlaps will be valid results
        for node_overlap in set1 & set2:
            res_paths.append([node_overlap])
        # explore the graph to find all possible paths
        while frontier:
            path = frontier.pop()
            last = path[-1]
            for new in self._data[last]:
                if new in set2:
                    res_paths.append(path + [new])
                elif new not in set1 and new not in path:
                    path2 = path + [new]
                    frontier.append(path2)
        res_paths.sort(key=lambda x: (len(x), x))
        return res_paths


    def cluster_nodes(self):
        """
        Put connected nodes in to clusters.
        Return a list [cluster1, cluster2, ..], each cluster is a set of nodes {node1, node2 ...}
        """
        clusters = []
        explored = set()
        nodes_iter = iter(self._data)
        while len(explored) < len(self._data):
            node1 = next(nodes_iter)
            if node1 not in explored:
                this_cluster = set()
                new_neighbors = self._data[node1]
                while new_neighbors:
                    this_cluster.update(new_neighbors)
                    explored.update(new_neighbors)
                    just_added = new_neighbors
                    new_neighbors = set(i for node in just_added for i in self._data[node] if i not in explored)
                clusters.append(this_cluster)
        return clusters

    def get_dihedrals(self):
        """
        Iterate over each bond, find distinct side atoms
        Return a list of distinct dihedral angles [(i,j,k,l), ..]
        """
        # indices i-j-k-l
        dihedral_list = []
        for j in self._data:
            j_neighbors = self._data[j]
            if len(j_neighbors) < 2: continue
            for k in j_neighbors:
                if k <= j: continue
                k_neighbors = self._data[k]
                for i in j_neighbors:
                    if i == k: continue
                    for l in k_neighbors:
                        if l == j or l == i: continue
                        dihedral_list.append((i,j,k,l))
        return dihedral_list

    def get_rings(self):
        """
        Find distinct rings in topology
        Return a list of rings, each ring is a list of atoms
        The canonical index of a ring is defined below:
        1. The first index is the smallest index in the ring
        2. The second index is the smallest neighbor of the first
        3. The rest of the index is consistent with the atom orders in the ring
        e.g. [[1,2,3,4], [4,5,6,7], ...]
        """
        distinct_rings = []
        paths = copy.deepcopy(self._bonds)
        while len(paths) > 0:
            path = paths.pop()
            origin = path[0]
            last_idx = path[-1]
            for neighbor in self._data[last_idx]:
                if neighbor in path:
                    if neighbor == origin and len(path) > 1:
                        # found a ring
                        if path[1] < last_idx:
                            # canonical ring index
                            # 1-2-3 will be kept, 1-3-2 will be skipped
                            distinct_rings.append(path)
                elif neighbor > origin:
                    # origin should be the smallest index in canonical ring
                    new_path = path + [neighbor]
                    paths.append(new_path)
        return distinct_rings

    def get_connected_nodes(self, node):
        """
        Get all nodes that are connected to node
        Return a set of nodes like {1,2,3,4}
        """
        result = set()
        new_connected = {node}
        while new_connected:
            result.update(new_connected)
            next_new = set()
            for new_node in new_connected:
                # add new nodes to next new
                next_new.update(self._data[new_node] - result)
            new_connected = next_new
        return result
