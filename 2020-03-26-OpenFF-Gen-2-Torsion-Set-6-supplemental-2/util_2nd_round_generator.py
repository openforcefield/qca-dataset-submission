from utils_torsion_dataset_generator import *

from IPython.display import SVG 
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

def visualize(molecule_attributes, mol_idx):
    # Use RDKit, visualize input molecule
    mapped_smiles = molecule_attributes[mol_idx]['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    mol = Chem.MolFromSmiles(mapped_smiles)
    mc = Chem.Mol(mol.ToBinary())
    rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(150,  150)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    display(SVG(svg.replace('svg:','')))

from openeye.oegraphsim import *
from openeye.oechem import *

def gen_sim_matrix(tid_molecules_list, tid, fptype= OEFPType_Tree):
    # Generate similarity matrix for input torsion parameter.
    fps = [] # list of finger prints
    mols = []
    molecule_list = [] # remove duplicate SMILES 
    for rotation in tid_molecules_list[tid]:
        if len(molecule_list) == 0:
            molecule_list.append(rotation['mol_index'])
        elif rotation['mol_index'] == molecule_list[-1]:
            continue
        else:
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
            if n!= m:
                sim_matrix[n,m] = OETanimoto(fps[n], fps[m])
    return sim_matrix, molecule_list

def convert_sim_matrix(sim_matrix):
    # Convert similiarity matrix to distance matrix 
    dis_matrix = np.ones_like(sim_matrix)-sim_matrix
    for i in range(len(dis_matrix)):
        dis_matrix[i][i] = 0
    return dis_matrix


from sklearn.cluster import DBSCAN, MeanShift
from sklearn import metrics

# def clustering(dis_matrix, eps):
#     clustering  = DBSCAN(eps = eps, min_samples= 2, metric = 'precomputed')
#     # Fit clustering
#     db = clustering.fit(dis_matrix)    
#     # Pull labels
#     labels = db.labels_
#     # Number of clusters in labels, ignoring noise if present.
#     n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
#     n_noise_ = list(labels).count(-1)
#     print(f' Estimated number of clusters: {n_clusters_} (tot:{len(dis_matrix)}, noise pts: {n_noise_})' )
#     return labels

# def gen_cluster_dict(tid_molecules_list, tid, fptype=OEFPType_MACCS166, epsilon=0.5):
#     sim_matrix, molecule_list = gen_sim_matrix(tid_molecules_list,tid, fptype)
#     dis_matrix = convert_sim_matrix(sim_matrix)
#     labels = clustering(dis_matrix, epsilon) 
#     cluster_dict = {}
#     for idx, mol_idx in enumerate(molecule_list):
#         cluster_dict[mol_idx] = labels[idx]
#     return cluster_dict, set(labels)

# def gen_tid_clusters_list(tid_molecules_list, fptype=OEFPType_MACCS166, epsilon=0.5):
#     tid_clusters_list = defaultdict(list)
#     for tid, molecules_list in tid_molecules_list.items():
#         tid_clusters_list[tid] = []
#         if len(molecules_list) == 0:
#             cluster_dict = {}
#         else: 
#             cluster_dict, labels  = gen_cluster_dict(tid_molecules_list, tid, fptype=OEFPType_MACCS166, epsilon=0.5)
#             if len(labels) > 1: 
#                 for label in labels: 
#                     if label != -1:
#                         clustered_dict = {'cluster_label': label, 'torsions':[]}
#                         mol_idx_list = [mol_idx  for (mol_idx, val) in cluster_dict.items() if val == label]
#                         for molecule in molecules_list:
#                             if molecule['mol_index'] in mol_idx_list: 
#                                 clustered_dict['torsions'].append(molecule)
#                         tid_clusters_list[tid].append(clustered_dict)
#             else: 
#                 for label in labels: 
#                     clustered_dict = {'cluster_label': label, 'torsions':[]}
#                     mol_idx_list = [mol_idx  for (mol_idx, val) in cluster_dict.items() if val == label]
#                     for molecule in molecules_list:
#                         if molecule['mol_index'] in mol_idx_list: 
#                             clustered_dict['torsions'].append(molecule)
#                     tid_clusters_list[tid].append(clustered_dict)
#     return tid_clusters_list

def clustering_mod(dis_matrix):
    # input distance matrix, cluster the matrix and return list of cluster labels.
    # Case 1. Identical elements: Ncluster=1
    if np.count_nonzero(dis_matrix) == 0:
        labels = [0  for  i in range(len(dis_matrix))]
        n_clusters_  = 1
        n_noise_ = 0
        if len(dis_matrix) > 1:
            print('All the molecules are the same? ')
        print(f'tot: {len(dis_matrix)}, Ncluster: {n_clusters_}')
    # Case 2. Two different elements: Ncluster=2
    elif len(dis_matrix)<3: 
        labels = [i for i in range(len(dis_matrix))]
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = 0
        print(f'tot: {len(dis_matrix)}, Ncluster: {n_clusters_}')
    # Case 3. More than 2 different elements: Ncluster>=2
    # If needed, cluster with varying clustering parameters(eps, min_samples) so that it is clustered into at least two clusters.
    else:
        # 3-1. First try eps=0.4, min_samples=2 (epsilon=0.4: reasonable value to separate distant chemistry)
        eps = 0.4
        clustering  = DBSCAN(eps = eps, min_samples= 2, metric = 'precomputed')
        db = clustering.fit(dis_matrix)    
        labels = db.labels_
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)
        # 3-2. Check number of clusters and if the number is greater than 5, use larger epsilon(0.5) to decrease the number of clusters 
        if n_clusters_  > 5:
            eps = 0.5
            clustering  = DBSCAN(eps = eps, min_samples= 2, metric = 'precomputed')
            db = clustering.fit(dis_matrix)    
            labels = db.labels_
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            n_noise_ = list(labels).count(-1)     
        # 3-3. Check if the number of clusters is smaller than one. If so, try epsilon from 0.5 to 0.1 until it is separated into least two clusters 
        if n_clusters_ < 2: 
            delta_eps = - 0.05
            enough_size = False
            while enough_size is False and 0.5 - delta_eps > 0.1:
                delta_eps += 0.05
                eps =  0.5 - delta_eps
                clustering  = DBSCAN(eps = eps, min_samples= 1, metric = 'precomputed')
                db = clustering.fit(dis_matrix)    
                labels = db.labels_
                n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
                n_noise_ = list(labels).count(-1)
                if n_clusters_ > 1:
                    enough_size = True                
        # 3-4. Print the result for checking the clustering result
        print(f'tot: {len(dis_matrix)}, eps: {eps:.2f}, Ncluster: {n_clusters_},labels: {labels}')
    return labels

def gen_cluster_dict_mod(tid_molecules_list, tid, fptype=OEFPType_MACCS166):
    #  input list of molecules for input torsion parameters, cluster the list and return dictionary storing the clustering result
    sim_matrix, molecule_list = gen_sim_matrix(tid_molecules_list,tid, fptype)
    dis_matrix = convert_sim_matrix(sim_matrix)
    labels = clustering_mod(dis_matrix) 
    cluster_dict = {}
    for idx, mol_idx in enumerate(molecule_list):
        cluster_dict[mol_idx] = labels[idx]
    return cluster_dict, set(labels)

def gen_tid_clusters_list_mod(tid_molecules_list, fptype=OEFPType_MACCS166):
    # generate  tid_clusters_list[tid] = [ {'cluster_label': N, 'torsions': [...]}]
    tid_clusters_list = defaultdict(list)
    for tid, molecules_list in tid_molecules_list.items():
        print(f'# {tid}')
        tid_clusters_list[tid] = []
        if len(molecules_list) == 0:
            cluster_dict = {}
        else: 
            cluster_dict, labels  = gen_cluster_dict_mod(tid_molecules_list, tid, fptype=OEFPType_MACCS166,)
            if len(labels) > 1: 
                for label in labels: 
                    if label != -1:
                        clustered_dict = {'cluster_label': label, 'torsions':[]}
                        mol_idx_list = [mol_idx  for (mol_idx, val) in cluster_dict.items() if val == label]
                        for molecule in molecules_list:
                            if molecule['mol_index'] in mol_idx_list: 
                                clustered_dict['torsions'].append(molecule)
                        tid_clusters_list[tid].append(clustered_dict)
            else: 
                for label in labels: 
                    clustered_dict = {'cluster_label': label, 'torsions':[]}
                    mol_idx_list = [mol_idx  for (mol_idx, val) in cluster_dict.items() if val == label]
                    for molecule in molecules_list:
                        if molecule['mol_index'] in mol_idx_list: 
                            clustered_dict['torsions'].append(molecule)
                    tid_clusters_list[tid].append(clustered_dict)
    return tid_clusters_list

def draw_table(dis_matrix):
      df = pd.DataFrame(np.array(dis_matrix), columns=[f'{i+1}' for i in range(len(dis_matrix))])
      vals = np.around(df.values, 2)
      fig, ax = plt.subplots()
      ax.axis('off')
      the_table=ax.table(cellText=vals, rowLabels=df.columns, colLabels=df.columns, loc='center', cellColours=plt.cm.RdYlGn_r(df),animated=True)
      the_table.set_fontsize(14)
      the_table.scale(1.5, 1.5)  
      plt.show()

def  find_reusable_cluster(tid_clusters_list, tid_calculated_molecules_list): 
    # tid_clusters_list_detail[tid] = [ {'cluster_label': N, 'torsions': [...], 'reusable': False or torsion_info}, ...] 
    tid_clusters_list_detailed = defaultdict(list)
    for tid,  clusters in tid_clusters_list.items():
        tid_clusters_list_detailed[tid] = []
        for cluster in clusters:
            detailed = {'cluster_label':cluster['cluster_label'],  'torsions': cluster['torsions']}
            detailed['reusable'] = False
            for torsion in cluster['torsions']:
                if any({'mol_index': torsion['mol_index'], 'indices': torsion['indices']} in lst for lst in tid_calculated_molecules_list.values()):
                    detailed['reusable'] = torsion
            tid_clusters_list_detailed[tid].append(detailed)
    return tid_clusters_list_detailed


def gen_graph_for_2nd_round(tid_clusters_list_detailed):
    graph = defaultdict(list)
    graph_single_coverage_set = defaultdict(list)
    graph_multiple_coverage_sets = defaultdict(list)
    graph_reusable_set = defaultdict(list)
    for tid, clusters in tid_clusters_list_detailed.items():
        graph[tid] = []
        for cluster in clusters: 
            cluster_dict = {'cluster_label': cluster['cluster_label'], 'sets': []}
            if cluster['reusable'] == False:
                for torsion in cluster['torsions']:
                    if set(torsion['covered_tids']) not in cluster_dict['sets']:
                        cluster_dict['sets'].append(set(torsion['covered_tids']))
                graph[tid].append(cluster_dict)
            else: 
                cluster_dict['sets'] = [set(cluster['reusable']['covered_tids'])]
                graph_reusable_set[tid] = [cluster_dict]
    for tid, clusters in graph.items():
        if len(clusters) == 0: 
            continue
        elif len(clusters) == 1 and len(clusters[0]['sets']) == 1:
            graph_single_coverage_set[tid] = clusters
        else:
            graph_multiple_coverage_sets[tid] = clusters
        
    return graph_reusable_set, graph_single_coverage_set, graph_multiple_coverage_sets

def find_minimum_degeneracy_for_2nd_round(graph_reusable_set, graph_single_coverage_set, graph_multiple_coverage_sets):
    # selected[tid] = [{'cluster_label': cluster['cluster_label'], 'selected_set': selected set from ['sets'] }]
    final_overlap = 1000
    trial = 0 
    selected = defaultdict()
    overlap_history = []
    coverage_history = []
    final_coverage = 0
    while trial <1e+04: 
        selected_candidate = defaultdict(list)
        trial += 1
        arrows = []
        for tid, clusters in graph_multiple_coverage_sets.items():
            for cluster in clusters:
                cluster_dict = {'cluster_label': cluster['cluster_label']}
                subgroups = cluster['sets']
                subgroup = random.choice(subgroups)
                cluster_dict['selected_set'] = subgroup
                selected_candidate[tid].append(cluster_dict)
                for tid2 in subgroup:
                    arrows.append([tid, tid2])
        for tid, clusters in graph_single_coverage_set.items():
            cluster= clusters[0]
            subgroups = cluster['sets']
            selected_candidate[tid].append({'cluster_label': cluster['cluster_label'], 'selected_set': subgroups[0]})
            for tid2 in subgroups[0]:
                arrows.append([tid, tid2])
        for tid, clusters in graph_reusable_set.items():
            cluster= clusters[0]
            subgroups = cluster['sets']
            selected_candidate[tid].append({'cluster_label': cluster['cluster_label'], 'selected_set': subgroups[0]})
            for tid2 in subgroups[0]:
                arrows.append([tid, tid2])
        n_overlap = 0
        for idx, (u, v) in enumerate(arrows):
            if [v,u] in arrows[idx+1:]:
                n_overlap += 1
        coverage = set()
        for arrow in arrows:
            for i in arrow:
                coverage.add(i)
        coverage = len(coverage)
        coverage_history.append(coverage)
        overlap_history.append(n_overlap)
        if n_overlap < final_overlap:
            final_overlap = n_overlap
            selected = selected_candidate
            final_coverage = coverage
        else:
            continue
   
    return selected , final_coverage, final_overlap, coverage_history, overlap_history

def draw_graph_for_2nd_round(selected, graph_single_coverage_set, graph_reusable_set):
    plt.figure(figsize=(20,20)) 
    G = nx.MultiDiGraph()
    for tid, clusters in selected.items():
        for cluster in clusters:
            if tid not in graph_single_coverage_set.keys() and tid not in graph_reusable_set.keys():
                for tid2 in cluster['selected_set']:
                    G.add_edge(tid, tid2, weight= 0.3, direction=True)       

    for tid, clusters in graph_single_coverage_set.items():
        for cluster in clusters:
            for tid2 in cluster['sets'][0]:
                G.add_edge(tid, tid2, weight= 1, direction=True)   
    for tid, clusters in graph_reusable_set.items():
        for cluster in clusters:
            for tid2 in cluster['sets'][0]:
                G.add_edge(tid, tid2, weight= 2, direction=True) 

    fixed_edges = [(u, v) for (u, v, d) in G.edges(data=True) if 1.0>= d['weight'] > 0.5]
    selected_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= 0.5]
    reused_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] >= 2]
    
    df = pd.DataFrame(index=G.nodes(), columns=G.nodes())
    for row, data in nx.shortest_path_length(G):
        for col, dist in data.items():
            df.loc[row,col] = dist

    df = df.fillna(df.max().max())

    pos = nx.kamada_kawai_layout(G, dist=df.to_dict())

    nx.draw_networkx_nodes(G, pos, node_size=1000, alpha=0.5)

    # edges
    nx.draw_networkx_edges(G, pos, edgelist=fixed_edges,
                        alpha=1, edge_color='k', arrowsize=20, arrowstyle='fancy')  
    nx.draw_networkx_edges(G, pos, edgelist=reused_edges, alpha=1, edge_color='b', arrowsize=20, arrowstyle='fancy')
    nx.draw_networkx_edges(G, pos, edgelist=selected_edges,
                            alpha=0.5, edge_color='r', style='dashed', arrowsize=20, arrowstyle='fancy')  

    # labels
    nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')

    plt.axis('off')
    plt.show()         


def select_rotations_for_2nd_round(tid_clusters_list_detailed, selected, molecules_list_dict, tid_calculated_molecules_list=None, molecules_list_dict_from_td=None,  first_round_tid_calculated_molecules_list=None, first_round_molecules_list_dict_from_td=None):
    # selected_rotations[tid] = [{'cluster_label': cluster['cluster_label'], 'torsion': torsion_info}]
    print("\n## Selecting Torsions... ##\n" + '-'*90)    
    selected_rotations  = defaultdict(list)
    count = 0 
    saved_count  = 0
    molecules_list_dict_updated = copy.deepcopy(molecules_list_dict)
    for tid,  clusters in tid_clusters_list_detailed.items():
        for cluster in clusters: 
            # for each cluster! 
            selected_tid = next(item for item in selected[tid] if item['cluster_label'] == cluster['cluster_label'])
            pot = [molecule for molecule in cluster['torsions'] if set(molecule['covered_tids'] ) == set(selected_tid['selected_set'])]
            if len(pot) > 0:
                
                if cluster['reusable'] != False:
                    ##########
                    picked_rotation = cluster['reusable']
                    saved_count += 1
                    selected_rotations[tid].append({'cluster_label': cluster['cluster_label'], 'torsion': picked_rotation})
                    print(f'\n*{tid}({cluster["cluster_label"]}): Precalculated torsion scan is detacted from 1st round.')
                    print(f'      : {picked_rotation["mol_index"]:40s}, indices: {picked_rotation["indices"]}')
                    print(f'      : Update molecules_list_dict so that it has the same intial molecules.')
                    molecules_list_dict_updated[picked_rotation["mol_index"]] = first_round_molecules_list_dict_from_td[picked_rotation["mol_index"]]
                         
                else:
                    if tid_calculated_molecules_list == None: 
                        degenerate = False
                        for picked_molecule in pot:
                            for lst_already_selected in selected_rotations.values():
                                for already_selected_info in lst_already_selected:
                                    already_selected = already_selected_info['torsion']
                                    if already_selected['mol_index'] == picked_molecule['mol_index'] and set(already_selected['indices'][1:3]) == set(picked_molecule['indices'][1:3]):
                                        degenerate = True
                        if degenerate:
                            print(f'\n*{tid}({cluster["cluster_label"]}): Non selected since the randomly selected rotation is already included.')
                        else: 
                            count += 1
                            picked_rotation = random.choice(pot)
                            selected_rotations[tid].append({'cluster_label': cluster['cluster_label'], 'torsion': picked_rotation})
                            print(f'\n*{tid}({cluster["cluster_label"]}): {picked_rotation["mol_index"]:40s}, indices: {picked_rotation["indices"]}')

                    else:
                        # check if there's any pre-calculated torsions in pot
                        precalculated = []
                        for molecule in cluster['torsions']:
                            if any({'mol_index': molecule['mol_index'], 'indices': molecule['indices']} in lst for lst in tid_calculated_molecules_list.values()):
                                precalculated.append(molecule)
                        if len(precalculated) == 0:
                            degenerate = False
                            for picked_molecule in pot:
                                for lst_already_selected in selected_rotations.values():
                                    for already_selected_info in lst_already_selected:
                                        already_selected = already_selected_info['torsion']
                                        if already_selected['mol_index'] == picked_molecule['mol_index'] and set(already_selected['indices'][1:3]) == set(picked_molecule['indices'][1:3]):
                                            degenerate = True
                            if degenerate:
                                print(f'\n*{tid}({cluster["cluster_label"]}): Non selected since the randomly selected rotation is already included.')
                            else: 
                                count += 1
                                picked_rotation = random.choice(pot)
                                selected_rotations[tid].append({'cluster_label': cluster['cluster_label'], 'torsion': picked_rotation})
                                print(f'\n*{tid}({cluster["cluster_label"]}): {picked_rotation["mol_index"]:40s}, indices: {picked_rotation["indices"]}')
                        else: 
                            print(f'\n*{tid}({cluster["cluster_label"]}): Precalculated torsion scans are detacted. Choose one out of them.')
                            degenerate = False
                            for picked_molecule in pot:
                                for lst_already_selected in selected_rotations.values():
                                    for already_selected_info in lst_already_selected:
                                        already_selected = already_selected_info['torsion']
                                        if already_selected['mol_index'] == picked_molecule['mol_index'] and set(already_selected['indices'][1:3]) == set(picked_molecule['indices'][1:3]):
                                            degenerate = True
                            if degenerate:
                                print(f'      : Non selected since the randomly selected rotation is already included.')
                            else:
                                saved_count += 1
                                picked_rotation = random.choice(precalculated)
                                selected_rotations[tid].append({'cluster_label': cluster['cluster_label'], 'torsion': picked_rotation})
                                print(f'      : {picked_rotation["mol_index"]:40s}, indices: {picked_rotation["indices"]}')
                                print(f'      : Update molecules_list_dict so that it has the same intial molecules.')
                                molecules_list_dict_updated[picked_rotation["mol_index"]] = molecules_list_dict_from_td[picked_rotation["mol_index"]]
                            
    print(f'\nIn total, {count} new calculations are needed and {saved_count} calculations will be reused:)\n')
    print('-'*90)

    print("\n## Selected Torsion Coverage ##\n" + '-'*90)
    n_tot = len(tid_clusters_list_detailed.keys())
    covered = set()
    for tid, clusters in selected.items():
        covered.add(tid)
        for cluster in clusters:
            selected_set  = cluster['selected_set']
            for shared_tid in selected_set:
                covered.add(shared_tid)
    n_covered = len( covered )
    print(f'Coverage: {n_covered}/ {n_tot}')
    print('Uncovered tids:', set(tid_clusters_list_detailed.keys())- covered)
    print('-'*90)


    print("\n## *Final* Number of selected torsions for each torsion parameter ##\n" + '-'*90)
    print(f"{'ID':7s} {'# clusters'} {'# torsions'}")
    for idx, (tid, molecules_list) in enumerate(selected_rotations.items()):
        print(f'{tid:7s} {len(tid_clusters_list_detailed[tid]):>10} {len(molecules_list):>9}')
    print('-'*90)
    return selected_rotations, molecules_list_dict_updated

def gen_json_for_2nd_round(selected_molecules, molecule_attributes, molecules_list_dict_updated, output_json='selected_torsions.json'):
    torsions_dict = {}
    for tid, clusters  in selected_molecules.items():
        for cluster in clusters:
            selected_rotation  = cluster['torsion']
            mol_index = selected_rotation['mol_index']
            mol_attr = molecule_attributes[mol_index]
            mapped_smiles = mol_attr['canonical_isomeric_explicit_hydrogen_mapped_smiles']
            atom_indices = selected_rotation['indices']
            canonical_torsion_index = cmiles.utils.to_canonical_label(mapped_smiles, atom_indices)
            torsions_dict[canonical_torsion_index] = {
                'initial_molecules' : molecules_list_dict_updated[mol_index],
                'atom_indices'      : [atom_indices],
                'attributes'        : mol_attr,
                'tid'               : tid
            }
    with open(output_json, 'w') as jsonfile:
        json.dump(torsions_dict, jsonfile, indent=2)

