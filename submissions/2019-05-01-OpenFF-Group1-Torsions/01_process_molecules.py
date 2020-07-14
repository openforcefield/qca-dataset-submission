#!/usr/bin/env python
# coding: utf-8
import os
import sys
import cmiles
from openeye import oechem

def read_split_mols(filename):
    ifs = oechem.oemolistream()
    # open input file that contains many molecules
    ifs.open(filename)
    # get all molecules from sdf file
    mol_list = []
    for mol in ifs.GetOEMols():
        # explicit declaring oechem.OEGraphMol is needed
        oe_mol = oechem.OEMol(mol)
        # add explicit hydrogens
        oechem.OEAddExplicitHydrogens(oe_mol)
        mol_list.append(oe_mol)
    ifs.close()
    print(f"Loaded {len(mol_list)} molecules")
    return mol_list

def canonical_order_molecule_inplace(mol_list):
    """rearrange atoms in molecules to canonical ordering
    based on example from Chaya """
    if not mol_list: return
    toolkit = cmiles.utils._set_toolkit(mol_list[0])
    for oe_mol in mol_list:
        toolkit.canonical_order_atoms(oe_mol)
        # make sure stereochemistry is defined
        if not cmiles.utils.has_stereo_defined(oe_mol):
            oechem.OEPerceiveChiral(oe_mol)
            oechem.OE3DToAtomStereo(oe_mol)
            oechem.OE3DToBondStereo(oe_mol)

def write_split_molecules_in_folder(mol_list, folder_name):
    ofs = oechem.oemolostream()
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
    xyz_folder = os.path.join(folder_name, 'xyz')
    if not os.path.isdir(xyz_folder):
        os.mkdir(xyz_folder)
    mol2_folder = os.path.join(folder_name, 'mol2')
    if not os.path.isdir(mol2_folder):
        os.mkdir(mol2_folder)
    for idx, mol in enumerate(mol_list, 1):
        formula = oechem.OEMolecularFormula(mol)
        #filename = f'{idx:03d}_{formula}.mol2'
        filename = os.path.join(folder_name, f'{idx:03d}_{formula}.sdf')
        print(f'Writing {filename}')
        ofs.open(filename)
        oechem.OEWriteMolecule(ofs, mol)
        ofs.close()
        # write xyz files
        xyzfile = os.path.join(xyz_folder, f'{idx:03d}_{formula}.xyz')
        ofs.open(xyzfile)
        oechem.OEWriteMolecule(ofs, mol)
        ofs.close()
        # write mol2 files
        mol2file = os.path.join(mol2_folder, f'{idx:03d}_{formula}.mol2')
        ofs.open(mol2file)
        oechem.OEWriteMolecule(ofs, mol)
        ofs.close()

### functions below are for submit_torsiondrives.py

def read_sdf_to_fb_mol(filename):
    """ read sdf file and return ForceBalance.molecule.Molecule object """
    from forcebalance.molecule import Molecule, Elements
    import numpy as np
    mol_list = read_split_mols(filename)
    assert len(mol_list) == 1, 'file contains multiple molecules'
    oe_mol = mol_list[0]
    # create a new molecule
    fb_mol = Molecule()
    # load elems
    fb_mol.elem = [Elements[a.GetAtomicNum()] for a in oe_mol.GetAtoms()]
    noa = len(fb_mol.elem)
    # load coordinates
    coords_dict = oe_mol.GetCoords()
    fb_mol.xyzs = [np.array([coords_dict[i] for i in range(noa)])]
    # load bonds
    bonds, bond_orders = [], []
    for oe_bond in oe_mol.GetBonds():
        idx_a = oe_bond.GetBgnIdx()
        idx_b = oe_bond.GetEndIdx()
        bond = (idx_a, idx_b) if idx_a <= idx_b else (idx_b, idx_a)
        bonds.append(bond)
        bond_orders.append(oe_bond.GetOrder())
    fb_mol.bonds = bonds
    fb_mol.bond_orders = bond_orders
    # load atomic formal charges
    atomic_formal_charges = [a.GetFormalCharge() for a in oe_mol.GetAtoms()]
    molecular_charge = sum(atomic_formal_charges)
    fb_mol.Data['molecular_charge'] = molecular_charge
    fb_mol.Data['atomic_formal_charges'] = atomic_formal_charges
    # set the oe_mol as one attribute
    fb_mol.oe_mol = oe_mol
    # set the cmiles id
    mapped_smiles = cmiles.utils.mol_to_smiles(oe_mol)
    fb_mol.Data['cmiles_id'] = cmiles.generator.get_molecule_ids(mapped_smiles)
    return fb_mol

def generate_torsion_index(fb_mol, torsion):
    assert hasattr(fb_mol, 'oe_mol'), 'fb_mol does not have oe_mol attr, did you create it using read_sdf_to_fb_mol?'
    oe_mol = fb_mol.oe_mol
    mapped_smiles = cmiles.utils.mol_to_smiles(oe_mol)
    # create a new mol and convert mapped smiles back to mol
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, mapped_smiles)
    for a in mol.GetAtoms():
        if not a.GetMapIdx()-1 in torsion:
            a.SetMapIdx(0)
    for i, t in enumerate(torsion):
        a = mol.GetAtom(oechem.OEHasMapIdx(t+1))
        a.SetMapIdx(i+1)
    return oechem.OEMolToSmiles(mol)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help='Input sdf file')
    args = parser.parse_args()

    print(' '.join(sys.argv))

    mol_list = read_split_mols(args.infile)

    canonical_order_molecule_inplace(mol_list)

    write_split_molecules_in_folder(mol_list, 'processed_molecules')


if __name__ == '__main__':
    main()
