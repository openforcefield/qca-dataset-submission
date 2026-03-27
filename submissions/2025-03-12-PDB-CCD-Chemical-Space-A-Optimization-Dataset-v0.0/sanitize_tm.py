""" This module uses function checks that the structure has all coordinates and performs custom sanitation for ferrocene structures.
"""
import numpy as np

from rdkit import Chem

METALS_NUM = [12,21,22,23,24,25,26,27,57,28,29,30,39,40,41,42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80]

def assert_same_ring(mol, ind1, ind2, max_ring_size=6):
    """Determine whether two atoms are in the same chemical ring

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule to assess
        ind1 (int): Index of first atom of interest
        ind2 (int): Index of second atom of interest
        max_ring_size (int, optional): Maximum ring size to consider

    Returns:
        bool: True if the two indices are in the same ring
    """
    ring_info = mol.GetRingInfo()
    
    indices = []
    for ring in ring_info.AtomRings():
        if ind1 in ring and len(ring) <= max_ring_size:
            indices.extend(list(set(ring)))
    if not indices:
        return False
    else:
        return ind2 in indices


def _correct_ferrocene(mol, index):
    """Correct a ferrocene containing molecule to ensure that all hydrogens are included

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule of interest
        index (int): Index of the metal center of a ferrocene group

    Returns:
        new_mol (rdkit.Chem.rdchem.Mol): Output molecule with corrected ferrocene group
    """
    
    metal = mol.GetAtoms()[index]
    symbol = metal.GetSymbol()
    if symbol != "Fe":
        raise ValueError("Ferrocene must have an iron center.")
    if metal.GetDegree() != 10:
        raise ValueError("Fe would have 10 bonds if it's a ferrocene.")
    c_atoms = []
    for b in metal.GetBonds():
        carbon = b.GetBeginAtom() if b.GetBeginAtom().GetSymbol() != symbol else b.GetEndAtom()
        c_atoms.append(carbon.GetIdx())
        for bc in carbon.GetBonds():
            tmp_atm = bc.GetBeginAtom() if bc.GetBeginAtomIdx() != carbon.GetIdx() else bc.GetEndAtom()
            if assert_same_ring(mol, carbon.GetIdx(), tmp_atm.GetIdx()):
                bc.SetBondType(Chem.BondType.AROMATIC)
            else:
                bc.SetBondType(Chem.BondType.SINGLE)
        b.SetBondType(Chem.BondType.DATIVE)
        carbon.SetNoImplicit(False)
        if carbon.GetDegree() < 4:
            carbon.SetNumExplicitHs(1)
        carbon.UpdatePropertyCache(strict=False)
    mol = Chem.AddHs(mol, addCoords=True, explicitOnly=False, onlyOnAtoms=c_atoms)
    
    return mol


def find_missing_coords(mol, value=0):
    """Determine if an RDKit molecule has a relevant geometry

    In PDB CCD if the coordinates are missing, denoted by question marks in the cif, then the coordinate will be (0,0,0)
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule to assess
        value (float): Value used to compare to coordinates. 
        If the sum across all dimensions for one atom is equal to this value, then a coordinate is missing.

    Returns:
        bool: Whether missing coordinates were detected.
    """
    
    conf = mol.GetConformer()
    positions = conf.GetPositions()
    pos_sum = np.sum(np.abs(positions), axis=-1)

    return any(pos_sum == value)


def sanitize_complex(mol, overall_charge, verbose=False, value_missing_coord=0):
    """Separate Ligands from transition metal and determine appropriate charge

    Args:
        mol (rdkit.Chem.rdchem.Mol): Ligand molecule
        overall_charge (int): Overall charge of the complex
        verbose (bool): Print updates
        value_missing_coord (float): Value expected when the coordinates of an atom that were missing are summed.
        For example, an RDKit molecule builder may have used (0,0,0) or (NaN, NaN, NaN).

    Raises:
        ValueError: If molecule does not contain a transition metal

    Returns:
        mol (rdkit.Chem.rdchem.Mol): Customly sanitized ligand molecule
    """

    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#0]"))
    
    tmc_idx = None
    for a in mol.GetAtoms():
        a.SetIntProp("__origIdx", a.GetIdx())
        a.SetNoImplicit(True)
        if a.GetAtomicNum() in METALS_NUM:
            tmc_idx = a.GetIdx()
    if tmc_idx is None:
        raise ValueError("No transition metal found")

    mol = Chem.AddHs(mol, addCoords=True, explicitOnly=False)
    
    # Detect and correct special cases
    if mol.GetAtoms()[tmc_idx].GetDegree() == 10 and mol.GetAtoms()[tmc_idx].GetSymbol() == "Fe": # Detect ferrocene
        mol = _correct_ferrocene(mol, tmc_idx)
    
    missing_coord_indices = find_missing_coords(mol, value=value_missing_coord)
    if missing_coord_indices:
        raise ValueError("Molecule missing coordinates")

    return mol