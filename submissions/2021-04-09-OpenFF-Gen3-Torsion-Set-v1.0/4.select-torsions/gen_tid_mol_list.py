from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule 
from openforcefield.topology import Topology 
from collections import defaultdict
from tqdm import tqdm
from openeye import oechem

def gen_canonical_isomeric_smiles(oemol):
    # 1. Create an OpenFF molecule from the OpenEye molecule, guessing the
    #    stereochemistry if needed.
    oe_molecule = oechem.OEMol(oemol)
    try:
        molecule = Molecule.from_openeye(oe_molecule)
    except:
        molecule = Molecule.from_openeye(oe_molecule, allow_undefined_stereo=True)
        stereoisomers = molecule.enumerate_stereoisomers(
            undefined_only=True, max_isomers=1
        )
        if len(stereoisomers) > 0:
            molecule = stereoisomers[0]
    # 2. Canonically order the molecule
    molecule = molecule.canonical_order_atoms()
    # 3. Figure out which atoms in the canonical molecule should be tagged.
    mapped_smiles = oechem.OEMolToSmiles(oe_molecule)
    torsion_match = molecule.chemical_environment_matches(
        mapped_smiles
    )[0]
    # 4. Generate a canonical isomeric mapped smiles
    molecule.properties["atom_map"] = {j: i + 1 for i, j in enumerate(torsion_match)}
    center_bond = set(molecule.properties["atom_map"].keys())
    canonical_isomeric_smiles = molecule.to_smiles(
        isomeric=True, explicit_hydrogens=True, mapped=False
    )
    return molecule, canonical_isomeric_smiles, center_bond

# input smi file, output tid molecule list 
def list_matching_torsions(smi_file, forcefield):
    from fragmenter import chemi # chemi.file_to_oemols

    # generate oemols from smi file
    oemols = chemi.file_to_oemols(smi_file)
    # list of torsion parameters 
    ff_torsion_param_list = forcefield.get_parameter_handler('ProperTorsions').parameters

    # tid_molecules_list[tid] = [{'mol_index': mol_index, 'indices': indices, 'covered_tids':covered_tids}, ...]
    tid_molecules_list = {}
    failed_smi = []
    for torsion_param in ff_torsion_param_list:
        tid_molecules_list[torsion_param.id] = []
    
    for oemol in tqdm(oemols): 
        try: 
            off_mol, mol_index, center_bond = gen_canonical_isomeric_smiles(oemol)
            oemol=Molecule.to_openeye(off_mol)
        except: 
            failed_smi.append(oechem.OEMolToSmiles(oemol))
            continue 

        torsions_coverage = defaultdict(list)
        off_top = Topology.from_molecules(off_mol)
        center_tids = defaultdict(set)
        dihedrals = []
        for torsion_indices, torsion_param in forcefield.label_molecules(off_top)[0]['ProperTorsions'].items():
            i, j, k, l = torsion_indices
            if set([j,k]) == center_bond: 
                center_tids[tuple(sorted([j,k]))].add(torsion_param.id)
                torsions_coverage[torsion_param].append(torsion_indices)
                dihedrals.append(torsion_indices)
        if not check_connectivity(dihedrals, oemol):
            print(f'## {mol_index} has diff bond info in oemol and offmol...')
            continue
        filtered_torsions_coverage = filter_torsions_coverage(torsions_coverage, oemol) # check connectivity
        
        for idx, (tid, indices_list) in enumerate(filtered_torsions_coverage.items()):
            for idxx, indices in enumerate(indices_list):
                if idxx == 0: # count once 
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

    return tid_molecules_list, failed_smi

def get_torsion_definition(ff_torsion_param_list, tid):
    for torsion_param in ff_torsion_param_list:
        if torsion_param.id == tid:
            answer= torsion_param
    return answer 

def check_connectivity(dihedrals, oemol):
    bonds = []
    for bond in oemol.GetBonds():
        bonds.append(set([bond.GetBgnIdx(), bond.GetEndIdx()]))
    bonds_connected = True
    bonds_to_check = []
    for dihedral in dihedrals: 
        i,j,k,l = dihedral
        bonds_to_check.append(set([i,j]))
        bonds_to_check.append(set([j,k]))
        bonds_to_check.append(set([k,l]))
    bonds_to_check = list(set(frozenset(bond) for bond in bonds_to_check))
    for bond_to_check in bonds_to_check:
        if bond_to_check not in bonds:
            bonds_connected  = False
    return bonds_connected

from collections import Counter, defaultdict
import re, os, shutil
from openeye import oechem
from bond_graph import BondGraph

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
import math
from typing import List, Optional, TypeVar

def gen_pdf(
    tid_molecules_list: list,
    output_path: str,
    cols: int = 8,
    cell_width: int = 200,
    cell_height: int = 200,
):
    from openeye import oechem, oedepict

    itf = oechem.OEInterface()
    PageByPage = True
    suppress_h = True
    
    n = sum([len(molecules) for tid, molecules in tid_molecules_list.items()])
    rows = math.ceil(n / cols)

    image = oedepict.OEImage(cell_width * cols, cell_height * rows)
    grid = oedepict.OEImageGrid(image, rows, cols)

    opts = oedepict.OE2DMolDisplayOptions(
        grid.GetCellWidth(), grid.GetCellHeight(), oedepict.OEScale_AutoScale
    )
    opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
    opts.SetTitleLocation(oedepict.OETitleLocation_Bottom)

    count = 0


    for tid, molecules in tid_molecules_list.items():
        for torsion in molecules: 
            cell = grid.GetCell(count//cols + 1, count%cols + 1)
            smi = torsion['mol_index']
            atom_indices = torsion['indices']
            off_mol= Molecule.from_smiles(smi, allow_undefined_stereo=True)
            off_mol = off_mol.canonical_order_atoms()

            mol=Molecule.to_openeye(off_mol)

            title = '{} ({})'.format(tid, smi)
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

def main():
    import pickle
    forcefield = ForceField('result.offxml',allow_cosmetic_attributes=True)

    tid_molecules_list, failed_smi = list_matching_torsions(smi_file='smiles-to-keep.smi', forcefield=forcefield)
    pickle.dump(tid_molecules_list, open('tid_molecules_list.p','wb'))
    gen_pdf(tid_molecules_list, output_path='tid_molecules_list.pdf')
if __name__ == '__main__':
    main()
