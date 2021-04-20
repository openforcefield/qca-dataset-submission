import math
from typing import List, Optional, TypeVar

def gen_pdf(
    tid_clusters_list: list,
    output_path: str,
    cols: int = 8,
    cell_width: int = 200,
    cell_height: int = 200,
):
    from openeye import oechem, oedepict
    from openforcefield.topology import Molecule
    itf = oechem.OEInterface()
    PageByPage = True
    suppress_h = True
    
    n = sum([len(clusters) for tid, clusters in tid_clusters_list.items()])
    rows = math.ceil(n / cols)

    image = oedepict.OEImage(cell_width * cols, cell_height * rows)
    grid = oedepict.OEImageGrid(image, rows, cols)

    opts = oedepict.OE2DMolDisplayOptions(
        grid.GetCellWidth(), grid.GetCellHeight(), oedepict.OEScale_AutoScale
    )
    opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
    opts.SetTitleLocation(oedepict.OETitleLocation_Bottom)

    count = 0


    for tid, clusters in tid_clusters_list.items():
        for cluster in clusters: 
            torsions = cluster['torsions']
            label = cluster['cluster_label']
            torsions = cluster['torsions']
            for torsion in torsions: 
                cell = grid.GetCell(count//cols + 1, count%cols + 1)
                smi = torsion['mol_index']
                atom_indices = torsion['indices']
                # mol = oechem.OEGraphMol()
                # oechem.OESmilesToMol(mol, smi)
 
                off_mol= Molecule.from_smiles(smi, allow_undefined_stereo=True)
                off_mol = off_mol.canonical_order_atoms()
                mol=Molecule.to_openeye(off_mol)

                title = '{} ({})'.format(tid ,set(torsion['covered_tids']))
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
    tid_clusters_list = pickle.load(open('tid_clusters_list.p','rb'))
    gen_pdf(tid_clusters_list, 'selected_torsions.pdf')

if __name__=='__main__':
    main()
