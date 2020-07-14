
import fragmenter, cmiles, json
import tarfile
from collections import Counter, defaultdict
import re, os, shutil
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import random
from openeye import oechem
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule as Off_Molecule
from openforcefield.topology import Topology as Off_Topology
from qcelemental.models import Molecule
from bond_graph import BondGraph
import pickle
import numpy as np

from forcebalance.molecule import Molecule as FBMolecule
import mdtraj as md

import qcportal as ptl

import copy

# Initialize client
client = ptl.FractalClient()
ofs = oechem.oemolostream()

def parse_input(input_file, output_json='optimization_inputs.json'):
    # Read input smi file and generate oemols
    oemols = fragmenter.chemi.file_to_oemols(input_file)
                
    optimization_input = []
    processed_canonical_smiles = []
    # Expand states. 
    for mol in oemols:
        # Filter out single atom molecules.
        if mol.GetMaxAtomIdx() == 1:
            continue
        # Expand tautomeric states and stereoisomers.
        states = fragmenter.states.enumerate_states(mol,  stereoisomers=True, tautomers=False)
        for s in states:
            # Screen out states having valence that rdkit does not accept.
            try:
                cmiles_ids = cmiles.get_molecule_ids(s)
            except:
                continue
            canonical_smiles = cmiles_ids['canonical_smiles']
            if canonical_smiles in processed_canonical_smiles:
                continue
            else: 
                processed_canonical_smiles.append(canonical_smiles)
            mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
            m = cmiles.utils.load_molecule(s)
            # Screen out states Omega fails to process.
            try:   
                conformers = fragmenter.chemi.generate_conformers(m)
            except RuntimeError:
                continue
            qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
            optimization_input.append({'initial_molecules': qcschema_molecules,
                                    'cmiles_identifiers': cmiles_ids})  
            
    # Save curated molecule set into json file.
    with open(output_json, 'w') as f:
        json.dump(optimization_input, f, indent=2, sort_keys=True)
    return optimization_input

def read_aggregate_molecules(input_json):
    molecules_list_dict = defaultdict(list)
    molecule_attributes = {}
    # open json file
    if input_json.endswith(".tar") or input_json.endswith(".tar.gz"):
        extract_file = input_json.replace(".gz", "").replace(".tar", ".json")
        with tarfile.open(input_json, 'r') as infile:
            molecule_data_list = json.load(infile.extractfile(extract_file))
    else:
        with open(input_json) as infile:
            molecule_data_list = json.load(infile)
    # put molecules and attributes into molecules_list_dict
    molecule_hash = defaultdict(set) # use a dictionary to remove duplicates
    for mdata in molecule_data_list:
        initial_molecules = mdata['initial_molecules']
        cmiles_ids = mdata['cmiles_identifiers']
        index = cmiles_ids['canonical_isomeric_smiles']
        molecule_attributes[index] = cmiles_ids
        for m_json in initial_molecules:
            m_hash = Molecule.from_data(m_json).get_hash()
            # find duplicated molecules using their hash and skip them
            if m_hash not in molecule_hash[index]:
                molecule_hash[index].add(m_hash)
                molecules_list_dict[index].append(m_json)
    return molecules_list_dict, molecule_attributes

def smirnoff_analysis_torsions(forcefield, off_mol):
    torsions_coverage = defaultdict(list)
    off_top = Off_Topology.from_molecules(off_mol)
    center_tids = defaultdict(set)
    for torsion_indices, torsion_param in forcefield.label_molecules(off_top)[0]['ProperTorsions'].items():
        i, j, k, l = torsion_indices
        # if j < k: torsion_indices = (i, j, k, l)
        # else:
        #     torsion_indices = (l, k, j, i)
        # i, j, k, l = torsion_indices
        center_tids[tuple(sorted([j,k]))].add(torsion_param.id)
        torsions_coverage[torsion_param].append(torsion_indices)
    return torsions_coverage, center_tids

def filter_torsions_coverage(torsions_coverage, oemol):
    # Collect usuful information using BondGraph
    bonds = []
    for bond in oemol.GetBonds():
        bonds.append((bond.GetBgnIdx(), bond.GetEndIdx()))
    bond_graph = BondGraph(bonds)
    rings = bond_graph.get_rings()
    d_rings = defaultdict(set)
    for i_ring, ring in enumerate(rings):
        for atom_idx in ring:
            d_rings[atom_idx].add(i_ring) 
    elem_list = []
    for atom in oemol.GetAtoms():
        elem_list.append(atom.GetAtomicNum())
    # print('elem_list',elem_list)
    # Filter out (1) unwanted in-ring rotations (2) terminal H when terminal is not specified
    filtered_torsions_coverage = defaultdict(list)
    for torsion_param, indices_list in torsions_coverage.items():
        rotatable_bond = False
        heavy_atoms = 4
        
        # Screening out unwanted in-ring rotations
        smirks_mod = re.sub(':2](\(.*\))?', ':2]', torsion_param.smirks)
        smirks_chopped = re.split('\:2\]', smirks_mod)[1]
        central_bond = re.split('\[.*:3\]', smirks_chopped)[0]  

        if central_bond in ['-;@', '-@', ':','=,:', '@']:
            rotatable_bond = False
        else: 
            rotatable_bond = True

        if re.search("[^!]#1:1", torsion_param.smirks):
            if re.search("[^!]#1:4", torsion_param.smirks):
                heavy_atoms = 2
            else: 
                heavy_atoms = 3
        elif re.search("[^!]#1:4", torsion_param.smirks):
            heavy_atoms = 3
        # validation for each indices
        for indices in indices_list:
            valid1=False
            valid2=False
            check_elem = [elem_list[idx] for idx in indices]
            if heavy_atoms == 4:
                if not any(elem_idx == 1 for elem_idx in check_elem):
                    valid1 = True
            elif heavy_atoms == 3:
                if check_elem.count(1) == 1:
                    valid1 = True
            elif heavy_atoms == 2: 
                if not any(elem_idx == 1 for elem_idx in check_elem[1:3]):
                    valid1 = True
            if rotatable_bond == False:
                valid2=True
            else: 
                i, j, k, l = indices
                if  d_rings[j] & d_rings[k]:
                    continue
                else: 
                    valid2 = True

            if valid1 and valid2:
                filtered_torsions_coverage[torsion_param.id].append(indices)

    return filtered_torsions_coverage
    
def get_torsion_definition(ff_torsion_param_list, tid):
    for torsion_param in ff_torsion_param_list:
        if torsion_param.id == tid:
            answer= torsion_param
    return answer                   

def gen_tid_molecules_list(molecule_attributes, molecules_list_dict, forcefield): 
    # gen dictionary with keys, including all tids in the input forcefield
    ff_torsion_param_list = forcefield.get_parameter_handler('ProperTorsions').parameters

    tid_molecules_list = {}
    for torsion_param in ff_torsion_param_list:
        tid_molecules_list[torsion_param.id] = []

    for idx, (mol_index, mol_attr) in enumerate(molecule_attributes.items()):
        mapped_smiles = mol_attr['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        qcjson_mol = molecules_list_dict[mol_index][0]
        oemol = cmiles.utils.load_molecule(qcjson_mol)
        off_mol = Off_Molecule.from_openeye(oemol, allow_undefined_stereo=True)
        
        torsions_coverage, center_tids = smirnoff_analysis_torsions(forcefield, off_mol)
        filtered_torsions_coverage = filter_torsions_coverage(torsions_coverage, oemol)

        for tid, indices_list in filtered_torsions_coverage.items():
            
            for indices in indices_list:
                covered_tids = []
                i,j,k,l = indices
                tids  = center_tids[tuple(sorted([j,k]))]
                for i in tids: 
                    if i not in covered_tids:
                        covered_tids.append(i)
                tid_molecules_list[tid].append({'mol_index': mol_index, 'indices': indices, 'covered_tids':covered_tids})
    print("\n## Torsion parameter: matched molecules ##\n" + '-'*90)
    print(f"{'idx':<7} {'ID':7s} {'SMIRKS Pattern':70s} {'Number of molecules matched'}")
    for idx, (tid, molecules_list) in enumerate(tid_molecules_list.items()):
        torsion_param = get_torsion_definition(ff_torsion_param_list, tid)
        print(f'{idx:<7} {torsion_param.id:7s} {torsion_param.smirks:70s} {len(molecules_list)}')
    print('-'*90)

    return tid_molecules_list

def gen_tid_molecules_list_of_interest(molecule_attributes, molecules_list_dict, forcefield, tid_list): 
    # gen dictionary with keys, including all tids in the input forcefield
    ff_torsion_param_list = forcefield.get_parameter_handler('ProperTorsions').parameters
    ff_torsion_param_list_of_interest = []
    tid_molecules_list = {}
    for torsion_param in ff_torsion_param_list:
        if torsion_param.id in tid_list:
            ff_torsion_param_list_of_interest.append(torsion_param)
            tid_molecules_list[torsion_param.id] = []

    for idx, (mol_index, mol_attr) in enumerate(molecule_attributes.items()):
        mapped_smiles = mol_attr['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        qcjson_mol = molecules_list_dict[mol_index][0]
        oemol = cmiles.utils.load_molecule(qcjson_mol)
        off_mol = Off_Molecule.from_openeye(oemol, allow_undefined_stereo=True)
        
        torsions_coverage, center_tids = smirnoff_analysis_torsions(forcefield, off_mol)
        filtered_torsions_coverage = filter_torsions_coverage(torsions_coverage, oemol)

        for tid, indices_list in filtered_torsions_coverage.items():
            if tid in tid_list:
            
                for indices in indices_list:
                    covered_tids = []
                    i,j,k,l = indices
                    tids  = center_tids[(j,k)]
                    for i in tids: 
                        if i not in covered_tids:
                            covered_tids.append(i)
                    tid_molecules_list[tid].append({'mol_index': mol_index, 'indices': indices, 'covered_tids':covered_tids})
    print("\n## Torsion parameter: matched molecules ##\n" + '-'*90)
    print(f"{'idx':<7} {'ID':7s} {'SMIRKS Pattern':70s} {'Number of molecules matched'}")
    for idx, (tid, molecules_list) in enumerate(tid_molecules_list.items()):
        torsion_param = get_torsion_definition(ff_torsion_param_list_of_interest, tid)
        print(f'{idx:<7} {torsion_param.id:7s} {torsion_param.smirks:70s} {len(molecules_list)}')
    print('-'*90)

    return tid_molecules_list

def download_torsiondrive_data(dataset_name, output_pickle='torsiondrive_data.pickle'):
    # load dataset from public qcfractal server
    ds = client.get_collection("TorsionDriveDataset", dataset_name)
    spec_name = ds.list_specifications().index[0]
    print(f"Loading TorsionDrive Scans from [ {dataset_name} ] spec [{spec_name}]")
    print(f"Found {len(ds.df)} data entries")
    # load torsiondrive record ids from the dataset
    map_record_id_entry_index = {}
    for entry_index in ds.df.index:
        data_entry = ds.get_entry(entry_index)
        td_record_id = data_entry.object_map[spec_name]
        map_record_id_entry_index[td_record_id] = entry_index, data_entry.attributes
    print(f"Found {len(map_record_id_entry_index)} torsiondrive records")
    # query all torsiondrive records at the same time
    td_record_ids = list(map_record_id_entry_index.keys())
    torsiondrive_data = {}
    for i, td_record in enumerate(client.query_procedures(id=td_record_ids), 1):
        entry_index, attributes = map_record_id_entry_index[td_record.id]
        print(f"{i:5d} : {entry_index:50s} status {td_record.status}")
        if td_record.status == 'COMPLETE':
            try:
                torsiondrive_data[entry_index] = {
                    'initial_molecules': client.query_molecules(td_record.initial_molecule),
                    'final_molecules': td_record.get_final_molecules(),
                    'final_energies': td_record.get_final_energies(),
                    'final_gradients': {gid: np.array(res.return_result) for gid, res in td_record.get_final_results().items()},
                    'keywords': td_record.keywords.dict(),
                    'attributes': attributes
                }
            except:
                print( 'Problem in storing torsiondrive_data for ', entry_index)
    print(f'Downloaded torsion drive data for {len(torsiondrive_data)} completed entries')
    # save as pickle file
    with open(output_pickle, 'wb') as pfile:
        pickle.dump(torsiondrive_data, pfile)
    return torsiondrive_data

def download_torsiondrive_data2(dataset_name, output_pickle='torsiondrive_data.pickle'):
    # load dataset from public qcfractal server
    ds = client.get_collection("TorsionDriveDataset", dataset_name)
    spec_name = ds.list_specifications().index[0]
    print(f"Loading TorsionDrive Scans from [ {dataset_name} ] spec [{spec_name}]")
    print(f"Found {len(ds.df)} data entries")
    # load torsiondrive record ids from the dataset
    map_record_id_entry_index = {}
    for entry_index in ds.df.index:
        data_entry = ds.get_entry(entry_index)
        td_record_id = data_entry.object_map[spec_name]
        map_record_id_entry_index[td_record_id] = entry_index, data_entry.attributes
    print(f"Found {len(map_record_id_entry_index)} torsiondrive records")
    # query all torsiondrive records at the same time
    td_record_ids = list(map_record_id_entry_index.keys())
    torsiondrive_data = {}
    for i, td_record in enumerate(client.query_procedures(id=td_record_ids), 1):
        entry_index, attributes = map_record_id_entry_index[td_record.id]
        print(f"{i:5d} : {entry_index:50s} status {td_record.status}")
        # if td_record.status == 'COMPLETE':
        try:
            torsiondrive_data[entry_index] = {
                'initial_molecules': client.query_molecules(td_record.initial_molecule),
                'final_molecules': td_record.get_final_molecules(),
                'final_energies': td_record.get_final_energies(),
                'final_gradients': {gid: np.array(res.return_result) for gid, res in td_record.get_final_results().items()},
                'keywords': td_record.keywords.dict(),
                'attributes': attributes
            }
        except:
            print( 'Problem in storing torsiondrive_data for ', entry_index)
    print(f'Downloaded torsion drive data for {len(torsiondrive_data)} existing entries')
    # save as pickle file
    with open(output_pickle, 'wb') as pfile:
        pickle.dump(torsiondrive_data, pfile)
    return torsiondrive_data



def test_ff_mol2(test_ff, mol2_fnm):
    """
    Test creating system with mol2 file
    """
    from openforcefield.topology import Molecule as Off_Molecule
    from openforcefield.topology import Topology as Off_Topology
    try:
        off_molecule = Off_Molecule.from_file(mol2_fnm)
        off_topology = Off_Topology.from_molecules(off_molecule)
        test_ff.create_openmm_system(off_topology)
        molecule_labels = test_ff.label_molecules(off_topology)[0]
    except Exception as e:
        return False, str(e), None
    return True, '', molecule_labels

def check_Hbond(scan_fnm, top_fnm=None):
    """
    Check if the torsion scan contains conformers with internal hydrogen bonds
    """
    
    traj = md.load(scan_fnm, top=top_fnm)
    hbonds = md.baker_hubbard(traj)
    if len(hbonds) == 0:
        return True
    else:
        return False

def gen_tid_calculated_molecules_list(torsiondrive_data, forcefield, verbose = False):

    # gen dictionary with keys, including all tids in the input forcefield
    ff_torsion_param_list = forcefield.get_parameter_handler('ProperTorsions').parameters

    tid_calculated_molecules_list = {}
    molecules_list_dict_from_td = defaultdict = {}
    for torsion_param in ff_torsion_param_list:
        tid_calculated_molecules_list[torsion_param.id] = []
    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    os.mkdir('tmp')
    os.chdir('tmp')
    for entry_index, td_data in torsiondrive_data.items():
        # pick a single initial molecule
        qcmol = td_data['initial_molecules'][0]

        # write input.mol2 file
        qcjson_mol = qcmol.dict(encoding='json')
        oemol = cmiles.utils.load_molecule(qcjson_mol)
        ofs.open(f'input.mol2')
        oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()
        # test mol2 file
        success, msg, molecule_labels = test_ff_mol2(forcefield, 'input.mol2')
        if not success: 
            if verbose ==True:
                print('Error occured while testing input.mol2. Excluded in tid_calculated_molecules_list. ')
            continue
        # check if the torsion scan contains one or more conformers forming strong internal H bonds
        if success:
            # write conf.pdb file
            fbmol = FBMolecule(f'input.mol2')
            # list of grid ids sorted
            sorted_grid_ids = sorted(td_data['final_molecules'].keys())
            # write scan.xyz 
            target_mol = FBMolecule()
            target_mol.elem = fbmol.elem
            target_mol.xyzs = []
            target_mol.qm_energies = []
            target_mol.qm_grads = []
            for grid_id in sorted_grid_ids:
                grid_qc_mol = td_data['final_molecules'][grid_id]
                # convert geometry unit Bohr -> Angstrom
                geo = grid_qc_mol.geometry * 0.529177
                target_mol.xyzs.append(geo)
                # add energy and gradient
                target_mol.qm_energies.append(td_data['final_energies'][grid_id])
                target_mol.qm_grads.append(td_data['final_gradients'][grid_id])
            target_mol.write('scan.xyz')

            no_hbonds = check_Hbond(scan_fnm ='scan.xyz', top_fnm='input.mol2')
            if not no_hbonds:
                if verbose ==True:
                    print('Internal hydrogen bond detacted. Excluded in tid_calculated_molecules_list. ')
                success = False
        if success:
            mol_index = td_data['attributes']["canonical_isomeric_smiles"]
            indices = td_data['keywords']['dihedrals'][0]
            tid = molecule_labels['ProperTorsions'][tuple(indices)].id
            
            # qcschema_molecules = [qcmol.dict(encoding='json') for qcmol in td_data['initial_molecules']]
            tid_calculated_molecules_list[tid].append({'mol_index': mol_index, 'indices': indices})

            qcschema_molecules = []
            for qcmol in td_data['initial_molecules']:
                j_dict = qcmol.dict(encoding='json')
                qcschema_molecule = {'symbols': j_dict['symbols'], 'geometry': j_dict['geometry'], 'connectivity': j_dict['connectivity'],
                    'molecular_charge': j_dict['molecular_charge'], 'molecular_multiplicity': j_dict['molecular_multiplicity']}
                qcschema_molecules.append(qcschema_molecule)

            molecules_list_dict_from_td[mol_index] = qcschema_molecules
    print("\n## Available torsion scans from QCArchive ##\n" + '-'*90)      
    print(f"{'idx':<7} {'tid':7s}  {'Number of torsion scans'}")   
    for idx, (tid, molecules_list) in enumerate(tid_calculated_molecules_list.items()):
        if len(molecules_list) > 0:
            print(f'{idx:<7} {tid:7s}  {len(molecules_list)}')      
    print('-'*90)  
    os.chdir('..')
    shutil.rmtree('tmp')
    return tid_calculated_molecules_list, molecules_list_dict_from_td


def gen_graph(tid_molecules_list):
    graph = defaultdict(list)
    graph_single_coverage_set = defaultdict(list)
    graph_multiple_coverage_sets = defaultdict(list)
    for tid, molecules in tid_molecules_list.items():
        for molecule in molecules:
            if set(molecule['covered_tids']) not in graph[tid]:
                graph[tid].append(set(molecule['covered_tids']))
    for tid, sets in graph.items():
        if len(sets) == 0: 
            continue
        elif len(sets) == 1:
            graph_single_coverage_set[tid] = sets
        else:
            graph_multiple_coverage_sets[tid] = sets
    
    return graph, graph_single_coverage_set, graph_multiple_coverage_sets

def find_minimum_degeneracy(graph_single_coverage_set, graph_multiple_coverage_sets):
    old_n_overlap = 100
    trial = 0 
    selected = defaultdict()
    while trial <1e+04: 
        selected_candidate = defaultdict()
        trial += 1
        arrows = []
        for tid,  subgroups in graph_multiple_coverage_sets.items():
            subgroup = random.choice(subgroups)
            selected_candidate[tid]  = subgroup
            for tid2 in subgroup:
                arrows.append([tid, tid2])
        for tid, subgroups in graph_single_coverage_set.items():
            selected_candidate[tid] = subgroups[0]
            for tid2 in subgroups[0]:
                arrows.append([tid, tid2])
        n_overlap = 0
        for idx, (u, v) in enumerate(arrows):
            if [v,u] in arrows[idx+1:]:
                n_overlap += 1
        if n_overlap < old_n_overlap:
            old_n_overlap = n_overlap
            selected = selected_candidate
        else:
            continue
    return selected

# may need to check coverage 
import matplotlib.pyplot as plt
def find_minimum_degeneracy_check_coverage(graph_single_coverage_set, graph_multiple_coverage_sets):
    final_overlap = 100
    trial = 0 
    selected = defaultdict()
    overlap_history = []
    coverage_history = []
    final_coverage = 0
    while trial <1e+04: 
        selected_candidate = defaultdict()
        trial += 1
        arrows = []
        for tid,  subgroups in graph_multiple_coverage_sets.items():
            subgroup = random.choice(subgroups)
            selected_candidate[tid]  = subgroup
            for tid2 in subgroup:
                arrows.append([tid, tid2])
        for tid, subgroups in graph_single_coverage_set.items():
            selected_candidate[tid] = subgroups[0]
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

def draw_graph(selected, graph_single_coverage_set):
    plt.figure(figsize=(20,20)) 
    G = nx.MultiDiGraph()
    for tid, lst_tids in selected.items():
        if tid not in graph_single_coverage_set.keys():
            for tid2 in lst_tids:
                G.add_edge(tid, tid2, weight= 0.3, direction=True)       

    for tid, lst_tids in graph_single_coverage_set.items():
        for tid2 in lst_tids[0]:
            G.add_edge(tid, tid2, weight= 1, direction=True)   

    fixed_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0.5]
    selected_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= 0.5]
    
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
    nx.draw_networkx_edges(G, pos, edgelist=selected_edges,
                            alpha=0.5, edge_color='r', style='dashed', arrowsize=20, arrowstyle='fancy')  

    # labels
    nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')

    plt.axis('off')
    plt.show()          

def select_rotations(tid_molecules_list, selected, molecules_list_dict, tid_calculated_molecules_list=None, molecules_list_dict_from_td=None):
    print("\n## Selecting Torsions... ##\n" + '-'*90)    
    selected_rotations  = defaultdict()
    count = 0 
    saved_count  = 0
    molecules_list_dict_updated = copy.deepcopy(molecules_list_dict)
    for tid, molecules in tid_molecules_list.items():
        
        pot = [molecule for molecule in molecules if set(molecule['covered_tids'] ) == selected[tid]]
        if len(pot) > 0:
            if tid_calculated_molecules_list == None: 
                degenerate = False
                for picked_molecule in pot:
                    for already_selected in selected_rotations.values():
                        if already_selected['mol_index'] == picked_molecule['mol_index'] and set(already_selected['indices'][1:3]) == set(picked_molecule['indices'][1:3]):
                            degenerate = True
                if degenerate:
                    print(f'\n*{tid}: Non selected since the randomly selected rotation is already included.')
                else: 
                    count += 1
                    picked_rotation = random.choice(pot)
                    selected_rotations[tid] = picked_rotation
                    print(f'\n*{tid}: {picked_rotation["mol_index"]:40s}, indices: {picked_rotation["indices"]}')

            else:
                # check if there's any pre-calculated torsions in pot
                precalculated = []
                for molecule in molecules:
                    if any({'mol_index': molecule['mol_index'], 'indices': molecule['indices']} in lst for lst in tid_calculated_molecules_list.values()):
                        precalculated.append(molecule)
                if len(precalculated) == 0:
                    degenerate = False
                    for picked_molecule in pot:
                        for already_selected in selected_rotations.values():
                            if already_selected['mol_index'] == picked_molecule['mol_index'] and set(already_selected['indices'][1:3]) == set(picked_molecule['indices'][1:3]):
                                degenerate = True
                    if degenerate:
                        print(f'\n*{tid}: Non selected since the randomly selected rotation is already included.')
                    else: 
                        count += 1
                        picked_rotation = random.choice(pot)
                        selected_rotations[tid] = picked_rotation
                        print(f'\n*{tid}: {picked_rotation["mol_index"]:40s}, indices: {picked_rotation["indices"]}')
                else: 
                    print(f'\n*{tid}: Precalculated torsion scans detacted. Choose one out of them.')
                    degenerate = False
                    for picked_molecule in precalculated:
                        for already_selected in selected_rotations.values():
                            if already_selected['mol_index'] == picked_molecule['mol_index'] and set(already_selected['indices'][1:3]) == set(picked_molecule['indices'][1:3]):
                                degenerate = True
                    if degenerate:
                        print(f'    : Non selected since the randomly selected rotation is already included.')
                    else:
                        saved_count += 1
                        picked_rotation = random.choice(precalculated)
                        selected_rotations[tid] = picked_rotation
                        print(f'    : {picked_rotation["mol_index"]:40s}, indices: {picked_rotation["indices"]}')
                        print(f'    : Update molecules_list_dict so that it has the same intial molecules.')
                        molecules_list_dict_updated[picked_rotation["mol_index"]] = molecules_list_dict_from_td[picked_rotation["mol_index"]]
                    
    print(f'\nIn total, {count} new calculations are needed and {saved_count} calculations will be reused:)\n')
    print('-'*90)

    print("\n## Selected Torsion Coverage ##\n" + '-'*90)
    n_tot = len(tid_molecules_list.keys())
    covered = set()
    for tid, sets in selected.items():
        covered.add(tid)
        for shared_tid in sets:
            covered.add(shared_tid)
    n_covered = len( covered )
    print(f'Coverage: {n_covered}/ {n_tot}')
    print('Uncovered tids:', set(tid_molecules_list.keys())- covered)
    print('-'*90)

    return selected_rotations, molecules_list_dict_updated

def gen_json(selected_molecules, molecule_attributes, molecules_list_dict_updated, output_json='selected_torsions.json'):
    torsions_dict = {}
    for tid, selected_rotation in selected_molecules.items():
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


