def is_center_bond_single(smirks):
    import re
    # Screening out unwanted in-ring rotations
    smirks_mod = re.sub(':2](\(.*\))?', ':2]', smirks)
    smirks_chopped = re.split('\:2\]', smirks_mod)[1]
    central_bond = re.split('\[.*:3\]', smirks_chopped)[0]
    if central_bond in ['-;@', '-@', ':','=,:', '@','=','#']:
        rotatable_bond = False
    else:
        rotatable_bond = True
    return rotatable_bond

def split_smirks(smirks):
    if is_center_bond_single(smirks):
        import re
        smirks_mod = re.sub(':2](\(.*\))?', ':2]', smirks)
        frag2 = '[*]'+re.split('\:2\]', smirks_mod)[1]
        frag1 = re.split(':2](\(.*\))?', smirks)[0] +']-[*]'
        return frag1, frag2
    else: 
        return None
def get_smirks(ff_torsion_param_list, tid):
    for ff_torsion_param in ff_torsion_param_list: 
        if ff_torsion_param.id == tid: 
            smirks =ff_torsion_param.smirks
    return smirks

def substructure_search_rdkit(smi, target_smiles_mod):
    import re
    smi_mod = re.sub(r"([*])", "I", smi)
    from rdkit import Chem
    rdmol = Chem.MolFromSmiles(smi_mod)
    rdmol = Chem.AddHs(rdmol)
    Chem.SanitizeMol(rdmol)
    Chem.SetAromaticity(rdmol)
    Chem.Kekulize(rdmol)
    patt = Chem.MolFromSmarts(target_smiles_mod)
    return rdmol.HasSubstructMatch(patt)

def switch_terminal_atom(target_smiles):
    import re
    if target_smiles.startswith('[*]-'):
        target_smiles_mod = re.sub("^\[([*])\]-", "I", target_smiles)
    elif target_smiles.endswith('-[*]'):
        target_smiles_mod = re.sub("-\[([*])\]$", "I", target_smiles)
    return target_smiles_mod

def get_matching_substituents(smirks, sublist):
    
    if split_smirks(smirks) is not None:
        frag1, frag2 = split_smirks(smirks)
        target_smiles_mod1 = switch_terminal_atom(frag1)
        target_smiles_mod2 = switch_terminal_atom(frag2)
        selected_substituents = {frag1: set(), frag2: set()}
        for sub in sublist:
            if substructure_search_rdkit(sub, target_smiles_mod1):
                selected_substituents[frag1].add(sub)
            if substructure_search_rdkit(sub, target_smiles_mod2):
                selected_substituents[frag2].add(sub)
    return selected_substituents

import math
from typing import List, Optional, TypeVar

def smiles_to_image_grid_mod(
    smiles: set,
    output_path: str,
    cols: int = 8,
    cell_width: int = 200,
    cell_height: int = 200,
):
    from openeye import oechem, oedepict

    itf = oechem.OEInterface()
    PageByPage = True
    suppress_h = True
    rows = math.ceil(len(smiles) / cols)

    image = oedepict.OEImage(cell_width * cols, cell_height * rows)
    grid = oedepict.OEImageGrid(image, rows, cols)

    opts = oedepict.OE2DMolDisplayOptions(
        grid.GetCellWidth(), grid.GetCellHeight(), oedepict.OEScale_AutoScale
    )
    opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
    opts.SetTitleLocation(oedepict.OETitleLocation_Bottom)

    for i, (smi, cell) in enumerate(zip(smiles, grid.GetCells())):

        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        oedepict.OEPrepareDepiction(mol, False, suppress_h)
        disp = oedepict.OE2DMolDisplay(mol, opts)
        oedepict.OERenderMolecule(cell, disp)

    oedepict.OEWriteImage(output_path, image)

def main():
    # 1. find uncovered torsions
    import pickle
    round1 = pickle.load(open('tid_clusters_list.p','rb'))
    from openff.toolkit.typing.engines.smirnoff import ForceField
    forcefield = ForceField('result.offxml',allow_cosmetic_attributes=True)
    ff_torsion_param_list = forcefield.get_parameter_handler('ProperTorsions').parameters

    uncovered = []
    poorly_covered = []
    for tid, clusters in round1.items():
        if len(clusters) == 0:
            uncovered.append(tid)
        elif len(clusters) <3 :
            poorly_covered.append(tid)

    print(f'# uncovered: {len(uncovered)}, # poorly covered: {len(poorly_covered)}')

    single_uncovered = []
    for tid in uncovered: 
        smirks = get_smirks(ff_torsion_param_list, tid)
        if is_center_bond_single(smirks):
            single_uncovered.append(tid)
    print(f'# single uncovered: {len(single_uncovered)}')

    for count, uncovered_tid in enumerate(single_uncovered):
        smirks = get_smirks(ff_torsion_param_list, uncovered_tid)
        print(f' {count}: {uncovered_tid}, {smirks}')
        
    pickle.dump(single_uncovered, open('single_uncovered.p','wb'))
    
    # 2. Load the substituent list 
    substituents_filtered=pickle.load(open('substituents_filtered.p','rb'))

    # 3. Load additional substituent list (hand-picked)
    with open('supplemental_substituents.smi','r') as file:
        supplemental_substituents = [
            line for line in file.read().split("\n") if len(line) > 0
        ]
    smiles_to_image_grid_mod(supplemental_substituents, output_path='supplemental_substituents.pdf')

    # 4. combine the lists
    substituents_filtered_new = set(substituents_filtered)
    substituents_filtered_new.update(supplemental_substituents)

    print(f'# substituents: {len(substituents_filtered)}')
    print(f'# additional substituents: {len(supplemental_substituents)}')
    print(f'# tot: {len(substituents_filtered_new)}')

    # 5. gen dict[tid]={frag1:[sub1, sub2, ...], frag2: [suba, subb, ...]}
    selected_substituents_tot = dict()
    for uncovered_tid in single_uncovered:
        smirks = get_smirks(ff_torsion_param_list, uncovered_tid)
        selected_substituents = get_matching_substituents(smirks, substituents_filtered_new)
        selected_substituents_tot[uncovered_tid] = selected_substituents 

    import pickle
    pickle.dump(selected_substituents_tot, open('selected_substituent_dict.p','wb'))

if __name__ =='__main__':
    main()