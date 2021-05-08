
def combine_substituents(sub1, sub2):
    from openeye import oechem
    import re

    sub1 = re.sub(r"([*])", r"([1*])", sub1)
    reactant = oechem.OEMol()
    oechem.OESmilesToMol(reactant, sub1)

    for atom in reactant.GetAtoms():
        if atom.GetAtomicNum() == 0:
            for bond in atom.GetBonds():
                nbor = bond.GetNbr(atom)
                nbor.SetMapIdx(1)

    for bond in nbor.GetBonds():
        nbor2 = bond.GetNbr(nbor)
        if nbor2.GetAtomicNum() != 0:
            nbor2.SetMapIdx(3)

    sub2 = re.sub(r"([*])", r"[*:1]", sub2)
    rxn = oechem.OEUniMolecularRxn(f"[*:1]-[{1}*]>>{sub2}")
    rxn(reactant)

    for atom in reactant.GetAtoms():
        if atom.GetMapIdx() == 1:
            b1 = atom.GetIdx()
            for bond in atom.GetBonds():
                nbor = bond.GetNbr(atom)
                if nbor.GetMapIdx() == 0:
                    nbor.SetMapIdx(2)
                    b2 = nbor.GetIdx()
                    break ##
        elif atom.GetMapIdx() == 3:
            atom.SetMapIdx(0)

    product_oemol = oechem.OEMol(reactant)
    print(f'{sub1:30s} + {sub2:30s} -> {oechem.OEMolToSmiles(product_oemol)}')

    return product_oemol, oechem.OEMolToSmiles(product_oemol)

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
        center_bond = []
        for atom in mol.GetAtoms():
            if atom.GetMapIdx() >0:
                center_bond.append(atom.GetIdx())
        # print(smi, center_bond)
        assert len(center_bond) == 2, print(smi)
        oedepict.OEPrepareDepiction(mol, False, suppress_h)
        disp = oedepict.OE2DMolDisplay(mol, opts)

        # Highlight element of interest
        class NoAtom(oechem.OEUnaryAtomPred):
            def __call__(self, atom):
                return False
        class NoBond(oechem.OEUnaryBondPred):
            def __call__(self, bond):
                return False
        class CentralBondInTorsion(oechem.OEUnaryBondPred):

            def __call__(self, bond):
                return (bond.GetBgn().GetIdx() in center_bond) and (bond.GetEnd().GetIdx() in center_bond)

        atoms = mol.GetAtoms(NoAtom())
        bonds = mol.GetBonds(CentralBondInTorsion())
        abset = oechem.OEAtomBondSet(atoms, bonds)
        oedepict.OEAddHighlighting(disp, oechem.OEColor(oechem.OEMandarin), oedepict.OEHighlightStyle_BallAndStick, abset)

        oedepict.OERenderMolecule(cell, disp)

    oedepict.OEWriteImage(output_path, image)

def write_smi(full_smiles_set, fnm):
    with open(fnm, 'w') as f: 
        for i, smi in enumerate(full_smiles_set):        
            f.write('{}\n'.format(smi))
            
def write_rdkit_smi(full_smiles_set, fnm):
    from rdkit import Chem
    with open(fnm, 'w') as f: 
        for i, smi in enumerate(full_smiles_set):      
            rdmol = Chem.MolFromSmiles(smi)  
            rdsmi = Chem.MolToSmiles(rdmol)
            f.write('{}\n'.format(rdsmi))
        
def smi_to_oemol(smi):
    from openeye import oechem
    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, smi)
    oechem.OEAddExplicitHydrogens(oemol)
    return oemol

def filter_list_of_substituents(list_of_substituents):
    frag_list_filtered = []
    for sub in list_of_substituents:
        if sub.startswith('*C#'):
            print(f'- exclude {sub}')
        else:
            frag_list_filtered.append(sub)
    return sorted(frag_list_filtered)

def main():
    import argparse
    import pickle
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pfile", help="List of substituents")
    parser.add_argument('--gen2D', action='store_true')
    parser.add_argument('--pdf_name', default='molecules')
    args = parser.parse_args()

    import sys
    print('# '+' '.join(sys.argv))

    from constructure.utilities.openeye import smiles_to_image_grid
    from openeye import oechem
    from collections import defaultdict

    substituents_dict = pickle.load(open(args.pfile, "rb"))

    full_smiles_set = set()
    for tid, frags in substituents_dict.items(): 
        frag1 = filter_list_of_substituents(frags[list(frags.keys())[0]])
        frag2 = filter_list_of_substituents(frags[list(frags.keys())[1]])
        for sub1 in frag1: 
            for sub2 in frag2: 
                deriv_oemol, derivative = combine_substituents(sub1, sub2)
                full_smiles_set.add(derivative)

    if args.gen2D: 

        if args.pdf_name.endswith('.pdf'):
            fnm = args.pdf_name
        else:
            fnm = args.pdf_name + '.pdf'
        smiles_to_image_grid_mod(full_smiles_set, fnm, cols=8)

    write_smi(full_smiles_set, fnm= args.pdf_name + '_out.smi')
    # write_rdkit_smi(full_smiles_set, fnm= args.pdf_name + '_out.smi')

if __name__ == "__main__":
    main()
