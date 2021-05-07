from openeye.oegraphsim import *
from openeye.oechem import *
from collections import defaultdict
import numpy as np
from sklearn.cluster import DBSCAN, MeanShift
from sklearn import metrics

def gen_tid_clusters_list(tid_molecules_list, fptype=OEFPType_MACCS166, lim=10, select_option =0):
    # generate  tid_clusters_list[tid] = [ {'cluster_label': N, 'torsions': [...]}]
    tid_clusters_list = defaultdict(list)
    for tid, molecules_list in tid_molecules_list.items():
        print(f'tid: {tid}')
        tid_clusters_list[tid] = []
        if len(molecules_list) == 0: 
            cluster_dict = {}
        else: 
            sim_matrix, mol_idx_list = gen_sim_matrix(tid_molecules_list,tid, fptype)
            dis_matrix = convert_sim_matrix(sim_matrix)
            labels = cluster(dis_matrix, lim=lim) 

            selected = []
            if select_option==1:  # select core
                sim_idx_sums = sim_matrix.sum(axis=0) 
                for label in list(set(labels)):
                    indices = [idx for idx, l in enumerate(labels) if l == label]
                    sublist=[sim_idx_sums[idx] for idx in indices]
                    selected.append(indices[sublist.index(max(sublist))])

            elif select_option== 2: # select smallest
                for label in list(set(labels)):
                    indices = [idx for idx, l in enumerate(labels) if l == label]
                    sublist = [mol_idx for idx, mol_idx in enumerate(mol_idx_list) if idx in indices]
                    selected_idx = find_smallest(sublist)
                    selected.append(indices[selected_idx])
                    
            elif select_option ==0: # select off
                selected = list(range(len(labels)))

            else: 
                raise NotImplementedError

            if len(labels) > 1: 
                for label in list(set(labels)):
                    if label != -1:
                        clustered_dict = {'cluster_label': label, 'torsions':[]}
                        selected_mol_idx_list = [mol_idx for idx, (mol_idx, l) in enumerate(zip(mol_idx_list, labels)) if idx in selected and l == label]
                        for molecule in molecules_list:
                            if molecule['mol_index'] in selected_mol_idx_list: 
                                clustered_dict['torsions'].append(molecule)
                        tid_clusters_list[tid].append(clustered_dict)
            else: 
                for label in list(set(labels)):
                    clustered_dict = {'cluster_label': label, 'torsions':[]}
                    selected_mol_idx_list = [mol_idx for idx, (mol_idx, l) in enumerate(zip(mol_idx_list, labels)) if idx in selected and l == label]
                    for molecule in molecules_list:
                        if molecule['mol_index'] in selected_mol_idx_list: 
                            clustered_dict['torsions'].append(molecule)
                    tid_clusters_list[tid].append(clustered_dict)
    return tid_clusters_list     

def find_smallest(mol_idx_list):
    weights = []
    for idx, mol_idx in enumerate(mol_idx_list):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, mol_idx)
        weight = oechem.OECalculateMolecularWeight(mol)
        weights.append(weight)
    return weights.index(min(weights))

def cluster(dis_matrix, lim=10):
    # Case 1. Identical elements: Ncluster=1
    if np.count_nonzero(dis_matrix) == 0:
        labels = [0  for  i in range(len(dis_matrix))]
        n_clusters_  = 1
        n_noise_ = 0
        if len(dis_matrix) > 1:
            print('All the molecules are the same? ')
        print(f'tot: {len(dis_matrix)}, Ncluster: {n_clusters_}')

    # Case 2. len(dis_matrix)< lim +1
    elif len(dis_matrix)<lim +1: 
        labels = [i for i in range(len(dis_matrix))]
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = 0
        print(f'tot: {len(dis_matrix)}, Ncluster: {n_clusters_}')

    else: 
        # 3-1. First try eps=0.4, min_samples=2 (epsilon=0.4: reasonable value to separate distant chemistry)
        eps = 0.4
        min_samples = 2
        step_size = 0.02
        clustering  = DBSCAN(eps = eps, min_samples= min_samples, metric = 'precomputed')
        db = clustering.fit(dis_matrix)    
        labels = db.labels_
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)
        print(f'# eps: {eps:.2f}, Ncluster: {n_clusters_}')
        # 3-2. Check number of clusters and if the number is greater than 20, use larger epsilon(0.5) to decrease the number of clusters 
        if n_clusters_ > lim +1: 
            # eps = 0.5
            # clustering  = DBSCAN(eps = eps, min_samples= 2, metric = 'precomputed')
            # db = clustering.fit(dis_matrix)    
            # labels = db.labels_
            # n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            # n_noise_ = list(labels).count(-1)  
            delta_eps =  0 ## havent been tested yet
            enough_size = False
            while enough_size is False:
                delta_eps += step_size
                eps =  0.4 + delta_eps
                clustering  = DBSCAN(eps = eps, min_samples= min_samples, metric = 'precomputed')
                db = clustering.fit(dis_matrix)    
                labels = db.labels_
                n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
                n_noise_ = list(labels).count(-1)
                print(f'# eps: {eps:.2f}, Ncluster: {n_clusters_}')
                if n_clusters_ <= lim:
                    enough_size = True   
        # 3-3. Check if the number of clusters is smaller than one. If so, try epsilon from 0.5 to 0.1 until it is separated into least x clusters 
        elif n_clusters_ < lim-2: 
            delta_eps = 0
            enough_size = False
            while enough_size is False and eps > 0.1:
            # while enough_size is False and eps > 0.2:
                delta_eps += step_size
                eps =  0.4 - delta_eps
                clustering  = DBSCAN(eps = eps, min_samples= 1, metric = 'precomputed')
                db = clustering.fit(dis_matrix)    
                labels = db.labels_
                n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
                n_noise_ = list(labels).count(-1)
                print(f'# eps: {eps:.2f}, Ncluster: {n_clusters_}')
                if n_clusters_ >= lim -2:
                    enough_size = True                
        # 3-4. Print the result for checking the clustering result
        print(f'tot: {len(dis_matrix)}, eps: {eps:.2f}, Ncluster: {n_clusters_}')
    return labels

def gen_sim_matrix(tid_molecules_list, tid, fptype= OEFPType_Tree):
    # Generate similarity matrix for input torsion parameter.
    fps = [] # list of finger prints
    mols = []
    molecule_list = [] # remove duplicate SMILES 
    for rotation in tid_molecules_list[tid]:
        if rotation['mol_index'] not in molecule_list:
            molecule_list.append(rotation['mol_index'])

    for mol_index in molecule_list:
        mol = OEMol()
        OEParseSmiles(mol, mol_index)
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, fptype)
        fps.append(fp)
        mols.append(mol)
    # Generate similarity matrix
    N = len(mols)
    sim_matrix = np.zeros((N,N), np.float)
    for n in range(N):
        for m in range(N):
            sim_matrix[n,m] = OETanimoto(fps[n], fps[m])
    return sim_matrix, molecule_list

def convert_sim_matrix(sim_matrix):
    # Convert similiarity matrix to distance matrix 
    dis_matrix = np.ones_like(sim_matrix)-sim_matrix
    for i in range(len(dis_matrix)):
        dis_matrix[i][i] = 0
    return dis_matrix

import math
from typing import List, Optional, TypeVar

def draw_table(dis_matrix):
    import pandas as pd
    import matplotlib.pyplot as plt
    df = pd.DataFrame(np.array(dis_matrix), columns=[f'{i+1}' for i in range(len(dis_matrix))])
    vals = np.around(df.values, 2)
    fig, ax = plt.subplots()
    ax.axis('off')
    the_table=ax.table(cellText=vals, rowLabels=df.columns, colLabels=df.columns, loc='center', cellColours=plt.cm.RdYlGn_r(df),animated=True)
    the_table.set_fontsize(14)
    the_table.scale(1.5, 1.5)  
    plt.show()

                    
def gen_pdf(
    cluster_dic: list,
    output_path: str,
    cols: int = 8,
    cell_width: int = 200,
    cell_height: int = 200,
):
    from openeye import oechem, oedepict

    itf = oechem.OEInterface()
    PageByPage = True
    suppress_h = True
    n = sum([len(dic['torsions']) for dic in cluster_dic])
    rows = math.ceil(n / cols)

    image = oedepict.OEImage(cell_width * cols, cell_height * rows)
    grid = oedepict.OEImageGrid(image, rows, cols)

    opts = oedepict.OE2DMolDisplayOptions(
        grid.GetCellWidth(), grid.GetCellHeight(), oedepict.OEScale_AutoScale
    )
    opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
    opts.SetTitleLocation(oedepict.OETitleLocation_Bottom)

    # for i, (dic, cell) in enumerate(zip(cluster_dic, grid.GetCells())):
    count = 0
    for i, dic in enumerate(cluster_dic):
        if 'cluster_label' in dic: 
            label = dic['cluster_label']
            torsions = dic['torsions']
            for torsion in torsions: 
                cell = grid.GetCell(count//cols + 1, count%cols + 1)
                smi = torsion['mol_index']
                atom_indices = torsion['indices']
                mol = oechem.OEGraphMol()
                oechem.OESmilesToMol(mol, smi)

                title = 'label:{} {}'.format(label ,torsion['covered_tids'])
                mol.SetTitle(title)

                oedepict.OEPrepareDepiction(mol, False, suppress_h)
                disp = oedepict.OE2DMolDisplay(mol, opts)

                # Highlight element of interest
                class NoAtom(oechem.OEUnaryAtomPred):
                    def __call__(self, atom):
                        return False
                class AtomInTorsion(oechem.OEUnaryAtomPred):
                    def __call__(self, atom):
                        return atom.GetIdx() in atom_indices
                class NoBond(oechem.OEUnaryBondPred):
                    def __call__(self, bond):
                        return False
                class CentralBondInTorsion(oechem.OEUnaryBondPred):
                    def __call__(self, bond):
                        return (bond.GetBgn().GetIdx() in atom_indices[1:3]) and (bond.GetEnd().GetIdx() in atom_indices[1:3])

                atoms = mol.GetAtoms(AtomInTorsion())
                bonds = mol.GetBonds(NoBond())
                abset = oechem.OEAtomBondSet(atoms, bonds)
                oedepict.OEAddHighlighting(disp, oechem.OEColor(oechem.OEYellow), oedepict.OEHighlightStyle_BallAndStick, abset)

                atoms = mol.GetAtoms(NoAtom())
                bonds = mol.GetBonds(CentralBondInTorsion())
                abset = oechem.OEAtomBondSet(atoms, bonds)
                oedepict.OEAddHighlighting(disp, oechem.OEColor(oechem.OEMandarin), oedepict.OEHighlightStyle_BallAndStick, abset)

                oedepict.OERenderMolecule(cell, disp)
                count += 1
        elif 'mol_index' in dic: 
            cell = grid.GetCell(count//cols + 1, count%cols + 1)
            smi = dic['mol_index']
            atom_indices = dic['indices'] 
            mol = oechem.OEGraphMol()
            oechem.OESmilesToMol(mol, smi)

            oedepict.OEPrepareDepiction(mol, False, suppress_h)
            disp = oedepict.OE2DMolDisplay(mol, opts)

            # Highlight element of interest
            class NoAtom(oechem.OEUnaryAtomPred):
                def __call__(self, atom):
                    return False
            class AtomInTorsion(oechem.OEUnaryAtomPred):
                def __call__(self, atom):
                    return atom.GetIdx() in atom_indices
            class NoBond(oechem.OEUnaryBondPred):
                def __call__(self, bond):
                    return False
            class CentralBondInTorsion(oechem.OEUnaryBondPred):
                def __call__(self, bond):
                    return (bond.GetBgn().GetIdx() in atom_indices[1:3]) and (bond.GetEnd().GetIdx() in atom_indices[1:3])

            atoms = mol.GetAtoms(AtomInTorsion())
            bonds = mol.GetBonds(NoBond())
            abset = oechem.OEAtomBondSet(atoms, bonds)
            oedepict.OEAddHighlighting(disp, oechem.OEColor(oechem.OEYellow), oedepict.OEHighlightStyle_BallAndStick, abset)

            atoms = mol.GetAtoms(NoAtom())
            bonds = mol.GetBonds(CentralBondInTorsion())
            abset = oechem.OEAtomBondSet(atoms, bonds)
            oedepict.OEAddHighlighting(disp, oechem.OEColor(oechem.OEMandarin), oedepict.OEHighlightStyle_BallAndStick, abset)

            oedepict.OERenderMolecule(cell, disp)
            count += 1

        else:
            for cluster in cluster_dic[dic]:
                label = cluster['cluster_label']
                torsions = cluster['torsions']
                for torsion in torsions: 
                    cell = grid.GetCell(count//cols + 1, count%cols + 1)
                    smi = torsion['mol_index']
                    atom_indices = torsion['indices']
                    mol = oechem.OEGraphMol()
                    oechem.OESmilesToMol(mol, smi)

                    title = 'tid: {}, label:{}'.format(dic,label)
                    mol.SetTitle(title)

                    oedepict.OEPrepareDepiction(mol, False, suppress_h)
                    disp = oedepict.OE2DMolDisplay(mol, opts)

                    # Highlight element of interest
                    class NoAtom(oechem.OEUnaryAtomPred):
                        def __call__(self, atom):
                            return False
                    class AtomInTorsion(oechem.OEUnaryAtomPred):
                        def __call__(self, atom):
                            return atom.GetIdx() in atom_indices
                    class NoBond(oechem.OEUnaryBondPred):
                        def __call__(self, bond):
                            return False
                    class CentralBondInTorsion(oechem.OEUnaryBondPred):
                        def __call__(self, bond):
                            return (bond.GetBgn().GetIdx() in atom_indices[1:3]) and (bond.GetEnd().GetIdx() in atom_indices[1:3])

                    atoms = mol.GetAtoms(AtomInTorsion())
                    bonds = mol.GetBonds(NoBond())
                    abset = oechem.OEAtomBondSet(atoms, bonds)
                    oedepict.OEAddHighlighting(disp, oechem.OEColor(oechem.OEYellow), oedepict.OEHighlightStyle_BallAndStick, abset)

                    atoms = mol.GetAtoms(NoAtom())
                    bonds = mol.GetBonds(CentralBondInTorsion())
                    abset = oechem.OEAtomBondSet(atoms, bonds)
                    oedepict.OEAddHighlighting(disp, oechem.OEColor(oechem.OEMandarin), oedepict.OEHighlightStyle_BallAndStick, abset)

                    oedepict.OERenderMolecule(cell, disp)
                    count += 1
    oedepict.OEWriteImage(output_path, image)