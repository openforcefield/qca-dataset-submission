""" This module uses function architectures originally produced in [xyz2mol_tm](https://github.com/jensengroup/xyz2mol_tm/). However rather than using the Huckel method
an arrow pushing script is produced here with custom checks for ferrocene structures.
"""
import numpy as np

from rdkit import Chem
from rdkit.Chem import GetPeriodicTable, rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Geometry import Point3D

METALS_NUM = [12,21,22,23,24,25,26,27,57,28,29,30,39,40,41,
              42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80,
]
MetalNon_Hg = "[#3,#11,#12,#19,#13,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#57,#72,#73,#74,#75,#76,#77,#78,#79,#80]~[B,#6,#14,#15,#33,#51,#16,#34,#52,Cl,Br,I,#85]"
bond_order_dict = {
    "SINGLE": 1,
    "DOUBLE": 2,
    "AROMATIC": 1.5,
    "DATIVE": 0, # this is 1 in rdkit, but it should be zero to determine multiplicity
    "TRIPLE": 3,
}
bond_type_dict = {
    0: Chem.BondType.DATIVE,
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    1.5: Chem.BondType.AROMATIC,
    3: Chem.BondType.TRIPLE,
}
pt = GetPeriodicTable()
params = Chem.MolStandardize.rdMolStandardize.MetalDisconnectorOptions()
params.splitAromaticC = True
params.splitGrignards = True
params.adjustCharges = False

def get_atom_charge(atom):
    """Get effective atom charge from the atom default valance and the total bond orders.

    Args:
        atom (rdkit.Chem.rdchem.Atom): RDKit atom to assess

    Raises:
        ValueError: If provided atom could have multiple charge states given the bonding orders

    Returns:
        int: Formal charge of the atom
    """
    Ntotalbonds =  np.array(pt.GetValenceList(atom.GetAtomicNum()))
    if Ntotalbonds[0] == -1:
        Ntotalbonds = np.array([0])
    bond_orders = [bond_order_dict[b.GetBondType().name] for b in atom.GetBonds()]
    charge = sum(bond_orders) - Ntotalbonds
    charge = charge[np.where(np.min(np.abs(charge)) == np.abs(charge))[0]]
    if len(charge) > 1:
        if atom.HasProp('__origIdx'):
            raise ValueError(f"Atom {atom.GetSymbol()}, Ligand Index: {atom.GetIdx()}, Complex Index: "
                         f"{atom.GetIntProp('__origIdx')}, can have multiple charge states {charge}")            
        else:
            raise ValueError(f"Atom {atom.GetSymbol()}, Ligand Index: {atom.GetIdx()} "
                             f"can have multiple charge states {charge}")
            
    return int(charge[np.where(np.min(np.abs(charge)) == np.abs(charge))[0]])


def assess_atoms(mol):
    """Assess an RDKit molecule's atoms to determine which are charged/not satisfied.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule

    Returns:
        charge (int): Total charge of the molecule
        hanging_bonds (int): Number of handing bonds to contribute toward the oxidation state of the metal
        charged_atoms (dict): In-depth information on the atoms showing a charge/are not satisfied with full bonding
        
        - symbol (str): Atomic symbol
        - charge (int): Formal charge of the atom as defined by :func:`get_atom_charge`
        - bond_orders (list[tuple]): List of bonds to this atom. Each bond is represented by a tuple with:
        
            - bond_type (str): rdkit.Chem.BondType.name
            - bond_order (int): Custom integer value representing the bond type, the major difference from RDKit being that DATIVE=0
            - BeginAtom_Symbol (str): Atomic symbol of the "BeginAtom"
            - BeginAtomIdx (int): Molecule index of the "BeginAtom"
            - EndAtom_Symbol (str): Atomic symbol of the "EndAtom"
            - EndAtomIdx (int): Molecule index of the "EndAtom"

    """
    charged_atoms = {}
    total_charge = 0
    hanging_bonds = 0
    for a in mol.GetAtoms():
        charge = get_atom_charge(a)
        if charge != 0:
            total_charge += charge
            charged_atoms[a.GetIdx()] = {
                "symbol": a.GetSymbol(),
                "charge": charge,
                "bond_orders": [(
                    b.GetBondType().name,
                    bond_order_dict[b.GetBondType().name],
                    b.GetBeginAtom().GetSymbol(),
                    b.GetBeginAtomIdx(),
                    b.GetEndAtom().GetSymbol(),
                    b.GetEndAtomIdx(),
                ) for b in a.GetBonds()],
            }
            if charge < 0: # positively charged groups won't bind to the metal center
                hanging_bonds += abs(charge)
                
    return int(total_charge), int(hanging_bonds), charged_atoms

def get_target_atom(index, filtered_bonds, path_atoms):
    
    bonded_atoms_charges = [
        (a, get_atom_charge(a)) for b in filtered_bonds 
        for a in [b.GetBeginAtom(), b.GetEndAtom()] 
        if a.GetIdx() != index
    ]
    if len(bonded_atoms_charges) == 0: # all the neighbors are happy, let it be
        b_atom_ind = None
    elif len(bonded_atoms_charges) > 1:
        bonded_atoms_charges = [
            (a, x) for a, x in bonded_atoms_charges 
            if a.GetIdx() in path_atoms and get_atom_charge(a) != 0
        ]
        b_atom_ind = bonded_atoms_charges[0][0].GetIdx() # just take the first and maybe revise later
    else: # there is only one
        b_atom_ind = bonded_atoms_charges[0].GetIdx()
        
    return b_atom_ind
                    
def get_ligand_attributes(ligand_mol, constant_atoms=None, resolved_atoms=None, verbose=False):
    """Analyze default valence and bonds to determine ligand attributes
    
    Args:
        ligand_mol (rdkit.Chem.rdchem.Mol): Ligand molecule
        constant_atoms (list[int]): List of atom indices defining those expected to be constant
        resolved_atoms (list[int]): List of atom indices that were resolved to assess progress
        verbose (bool): Print updates. Default ``False``
        
    Returns:
        total_charge (int): Total charge of the ligand, including positive and negative
        hanging_bonds (int): Number of unused valencies that could be bonded to the 
        metal center or represent a negative charge.
        charged_atoms (dict[dict]): Dictionary of atom information by index as defined in :func:`assess_atoms`
        mol (rdkit.Chem.rdchem.Mol): Resolved Ligand molecule 
    """
    
    new_code = False
    ligand_mol = Chem.RWMol(ligand_mol)
    if not new_code:
        total_charge_after, hanging_bonds_after, charged_atoms_after = assess_atoms(ligand_mol)
        if verbose:
            print(f"Total ligand charge: {total_charge_after}, N_Hanging Bonds {hanging_bonds_after}")
            print(f"Charged atom info:")
            for x,y in charged_atoms_after.items():
                print(f"    {x}: {y}")
        return total_charge_after, hanging_bonds_after, charged_atoms_after, ligand_mol

    _, _, charged_atoms_before = assess_atoms(ligand_mol)
    ligand_mol_before = ligand_mol
    ligand_mol = Chem.RWMol(Chem.AdjustQueryProperties(ligand_mol))

    total_charge_after, hanging_bonds_after, charged_atoms_after = assess_atoms(ligand_mol)
    if not new_code:
        return total_charge_after, hanging_bonds_after, charged_atoms_after, ligand_mol

    if constant_atoms is None:
        constant_atoms = list(set(charged_atoms_before.keys()) & set(charged_atoms_after.keys()))
    unresolved_atoms = [x for x in charged_atoms_after.keys() if x not in charged_atoms_before]
    if not unresolved_atoms and constant_atoms is not None:
        unresolved_atoms = [x for x in charged_atoms_after.keys() if x not in constant_atoms]
    resolved_atoms = resolved_atoms if resolved_atoms else []
    resolved_atoms.extend([x for x in charged_atoms_before.keys() if x not in charged_atoms_after])

    if not resolved_atoms: # If no atoms were resolved, then prefer the existing structure
        total_charge_before, hanging_bonds_before, charged_atoms_before = assess_atoms(ligand_mol_before)
        if verbose:
            print(f"Total ligand charge: {total_charge_before}, N_Hanging Bonds {hanging_bonds_before}")
            print(f"Charged atom info:")
            for x,y in charged_atoms_before.items():
                print(f"    {x}: {y}")
        return total_charge_before, hanging_bonds_before, charged_atoms_before, ligand_mol_before
    
    # Presumably aromaticity was resolved, so those should be constant
    if len(unresolved_atoms) > 0: 
        atoms = ligand_mol.GetAtoms()
        neighbors = [
            x for jnd in resolved_atoms 
            for b in atoms[jnd].GetBonds() 
            for x in [b.GetBeginAtomIdx(), b.GetEndAtomIdx()]
            if x in unresolved_atoms
        ]
        if not neighbors:
            raise ValueError("It's expected that all unresolved atoms are connected to the resolved atoms."
                             f" Resolved: {resolved_atoms}, Unresolved: {unresolved_atoms}")

        for index in unresolved_atoms:
            unresolved_atom = atoms[index]
            if get_atom_charge(unresolved_atom) == 0: # could have been resolved with a previous atom
                continue
            bonds = [b for b in unresolved_atom.GetBonds()]
            bond_orders = [bond_order_dict[b.GetBondType().name] for b in bonds]
            filtered_bonds = [b for b, bo in zip(bonds, bond_orders) if bo > 1.5]
            path_atoms = list(set( # Get atoms in the path from atoms[index] to one of the constant_atoms
                y for x in constant_atoms 
                for y in Chem.GetShortestPath(ligand_mol, index, x)
            ))
            if len(filtered_bonds) == 1:
                b = filtered_bonds[0]
                b_atom_ind = b.GetBeginAtomIdx() if b.GetBeginAtomIdx() != index else b.GetEndAtomIdx()
                bond_order =  pt.GetDefaultValence(unresolved_atom.GetAtomicNum()) - sum([
                    bo for bo in bond_orders if bo <= 1.5
                ])
            elif len(filtered_bonds) == 0: # No filtered_bonds so all bonds are aromatic or single. Delete a single bond
                filtered_bonds = [b for b, bo in zip(bonds, bond_orders) if bo == 1]
                if len(filtered_bonds) == 0: # If all bonds are aromatic, and assume it's right
                    continue
                    
                b_atom_ind = get_target_atom(index, filtered_bonds, path_atoms)
                if b_atom_ind is not None:
                    bond_order =  pt.GetDefaultValence(unresolved_atom.GetAtomicNum()) - (sum([
                        bo for bo in bond_orders
                    ]) - 1) # for removed single bond
            else: # There is more than one bond with order > 1.5 (aromatic)
                b_atom_ind = get_target_atom(index, filtered_bonds, path_atoms)
                if b_atom_ind is not None:
                    bond_order =  pt.GetDefaultValence(unresolved_atom.GetAtomicNum()) - sum([
                        bo for b, bo in zip(bonds, bond_orders) if b_atom_ind not in [b.GetBeginAtomIdx(), b.GetEndAtomIdx()]
                    ])
                    
            if b_atom_ind is not None:
                ligand_mol.RemoveBond(index, b_atom_ind)
                ligand_mol.AddBond(index, b_atom_ind, bond_type_dict[bond_order])

        total_charge_after, hanging_bonds_after, charged_atoms_after, ligand_mol = get_ligand_attributes(
            ligand_mol, constant_atoms=constant_atoms, resolved_atoms=unresolved_atoms, verbose=verbose
        )
    
        total_charge_after, hanging_bonds_after, charged_atoms_after = assess_atoms(ligand_mol)
        
    if verbose:
        print(f"Total ligand charge: {total_charge_after}, N_Hanging Bonds {hanging_bonds_after}")
        print(f"Charged atom info:")
        for x,y in charged_atoms_after.items():
            print(f"    {x}: {y}")
    return total_charge_after, hanging_bonds_after, charged_atoms_after, ligand_mol


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
        new_index (int): The atomic index of the ferrocene metal center, if changed
        tm_ox (int): Oxidation state of the metal, ``tm_ox = 2`` for ferrocene
    """
    
    metal = mol.GetAtoms()[index]
    symbol = metal.GetSymbol()
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
    
    for a in mol.GetAtoms():
        a.SetIntProp("__origIdx", a.GetIdx())
        if a.GetAtomicNum() in METALS_NUM:
            # tm_atom = a.GetSymbol()
            new_index = a.GetIdx()
    
    return mol, new_index, 2


def compute_centroid_excluding(conformer, exclude_atoms):
    """Compute the centroid of a molecule while excluding specified atom indices.
    
    Args:
        conformer (rdkit.Chem.rdchem.Conformer): RDKit conformer with 3D coordinates
        exclude_atoms (list[int]): List of atom indices to exclude from centroid calculation
        
    Returns:
        Point3D: Centroid of the remaining atoms
    """
    positions = conformer.GetPositions()
    for i in range(len(positions)):
        if i in exclude_atoms:
            positions[i] = [np.nan, np.nan, np.nan]

    centroid = np.nanmean(positions, axis=0)
    return Point3D(*centroid)

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
    pos_sum = np.sum(positions, axis=-1)

    return any(pos_sum == value)
        
def fix_missing_coords(mol, tmc_idx, missing_coord_indices):
    """Add coordinates to RDKit molecule with missing coordinates

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit to be repaired
        tmc_idx (int): Atom index of the complex metal
        missing_coord_indices (list[int]): Atom indices for which to find coordinates

    """
    
    # Move bad atoms closer
    conformer = mol.GetConformer()
    center = compute_centroid_excluding(conformer, missing_coord_indices)
    for i, atm_idx in enumerate(missing_coord_indices):
        radius = 1
        tmp_coord = Point3D(*tuple(
            np.array([center.x, center.y, center.z])
            + np.random.rand(3) * 2 * radius - radius
        ))
        conformer.SetAtomPosition(atm_idx, tmp_coord)
    
    # Optimize
    ff = Chem.AllChem.UFFGetMoleculeForceField(mol)
    metal_atoms = list(set([
        x 
        for b in mol.GetAtoms()[tmc_idx].GetBonds() 
        for x in [b.GetBeginAtomIdx(), b.GetEndAtomIdx()]
    ]))
    overlap = list(set(metal_atoms) & set(missing_coord_indices))
    for atm_idx in metal_atoms:
        if atm_idx not in overlap:
            ff.AddFixedPoint(atm_idx)
    ff.Minimize(maxIts=200000)

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
        tm_info (tuple(rdkit.Chem.rdchem.Mol, int)): An RDKit Molecule of the transition metal center and it's formal charge
        ligand_info (list(tuple(rdkit.Chem.rdchem.Mol, int))): A list of ligands where each entry is composed of an
        RDKit Molecule of the ligand, its formal charge, number of handing bonds, and charged atom information as defined in
        :func:`assess_atoms`
        mol (rdkit.Chem.rdchem.Mol): Customly sanitized ligand molecule
    """
    tm_ox = None
    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#0]"))
    Chem.AddHs(mol, addCoords=True, explicitOnly=False)
    
    tmc_idx = None
    for a in mol.GetAtoms():
        a.SetIntProp("__origIdx", a.GetIdx())
        a.SetNoImplicit(True)
        if a.GetAtomicNum() in METALS_NUM:
            tmc_idx = a.GetIdx()
    if tmc_idx is None:
        raise ValueError("No transition metal found")
    
    # Detect and correct special cases
    if mol.GetAtoms()[tmc_idx].GetDegree() == 10: # Detect ferrocene
        mol, tmc_idx, tm_ox = _correct_ferrocene(mol, tmc_idx)
    
    missing_coord_indices = find_missing_coords(mol, value=value_missing_coord)
    if missing_coord_indices:
        raise ValueError("Molecule missing coordinates")
    #    mol = fix_missing_coords(mol, tmc_idx, missing_coord_indices)

    coordinating_atoms = np.nonzero(Chem.rdmolops.GetAdjacencyMatrix(mol)[tmc_idx, :])[
        0
    ]

    mdis = rdMolStandardize.MetalDisconnector(params)
    mdis.SetMetalNon(Chem.MolFromSmarts(MetalNon_Hg))
    frags = mdis.Disconnect(mol)
    frag_mols = rdmolops.GetMolFrags(frags, asMols=True, sanitizeFrags=False)
    if verbose:
        print(f"Along with the metal, there are {len(frag_mols)-1} ligands")

    total_lig_bonds = 0
    total_lig_charge = 0
    flag_tm = False
    lig_info = []
    for i, f in enumerate(frag_mols):
        m = Chem.Mol(f)
        atoms = m.GetAtoms()
        for atom in atoms:
            if atom.GetAtomicNum() in METALS_NUM:
                if len(atoms) > 1:
                    raise ValueError("Not all ligands were separated.")
                flag_tm = True
                tm_mol = Chem.RWMol(frag_mols[i])
                break
        else:
            if verbose:
                print(f"Ligand {i+1} of {len(frag_mols)-1}")
            total_charge, hanging_bonds, charged_atoms, m = get_ligand_attributes(m, verbose=verbose)
            
            lig_coordinating_atom_bonds = [
                charged_atoms[a.GetIdx()]["charge"]
                for a in m.GetAtoms()
                if (
                    a.GetIntProp("__origIdx") in coordinating_atoms
                    and a.GetIdx() in charged_atoms
                )
            ]

            total_lig_bonds += abs(sum(lig_coordinating_atom_bonds))
            total_lig_charge += total_charge
            lig_info.append((m, total_charge, hanging_bonds, charged_atoms))

    if not flag_tm:
        raise ValueError("No transition metal found")
    
    # Set complex charge
    elif tm_ox is None and overall_charge is None:
        tm_ox = total_lig_bonds
    elif tm_ox is None:
        # In case there is an ammonium or something changing the charge of the ligand
        tm_ox = overall_charge - (total_lig_charge + total_lig_bonds) + total_lig_bonds

    for a in tm_mol.GetAtoms():
        if a.GetAtomicNum() in METALS_NUM:
            a.SetFormalCharge(tm_ox)
            
    tmc_mol = reform_metal_complex((tm_mol, tm_ox), lig_info, coordinating_atoms, verbose=verbose)

    return (tm_mol, tm_ox), lig_info, tmc_mol

def reform_metal_complex(tm_info, lig_info, coordinating_atoms, verbose=True):
    """Reconnect ligands with the metal complex

    Args:
        tm_info (tuple(rdkit.Chem.rdchem.Mol, int)): An RDKit Molecule of the transition metal center and it's formal charge
        ligand_info (list(tuple(rdkit.Chem.rdchem.Mol, int))): A list of ligands where each entry is composed of an
        RDKit Molecule of the ligand, its formal charge, number of handing bonds, and charged atom information as defined in
        :func:`assess_atoms`
        coordinating_atoms (list[int]): List of atom indices that were connected to the metal center in the original complex,
        as determined from a custom set atom property, "__origIdx" 
        verbose (bool, optional): Print updates. Defaults to True.
    """
    # Reform complex with ligands
    tm = tm_info[0]
    tm_symbol = tm.GetAtoms()[0].GetSymbol()
    for lmol, _, _, _ in lig_info:
        tm = Chem.CombineMols(tm, lmol)

   # Add bonds
    tmc_mol = Chem.RWMol(tm)
    coordinating_atoms_idx = [
        a.GetIdx()
        for a in tmc_mol.GetAtoms()
        if a.GetIntProp("__origIdx") in coordinating_atoms
    ]
    tm_idx = [
        a.GetIdx() for a in tmc_mol.GetAtoms() if a.GetSymbol() == tm_symbol
    ][0]
    
    for i in coordinating_atoms_idx:
        bond = tmc_mol.GetBondBetweenAtoms(i, tm_idx)
        a = tmc_mol.GetAtoms()[i]
        bonds = [b for b in a.GetBonds()]
        bond_orders = [bond_order_dict[b.GetBondType().name] for b in bonds]
        bond_order =  pt.GetDefaultValence(a.GetAtomicNum()) - sum([
            bo for b, bo in zip(bonds, bond_orders) if tm_idx not in [b.GetBeginAtomIdx(), b.GetEndAtomIdx()]
        ])
        if bond:
            old_bo = bond_order_dict[bond.GetBondType().name]
            if bond_order != old_bo:
                print(f"Changing bond between {tm_idx} and {i} from {bond_type_dict[old_bo]}"
                      f" to {bond_type_dict[bond_order]}")
            tmc_mol.RemoveBond(i, tm_idx)
        tmc_mol.AddBond(i, tm_idx, bond_type_dict[bond_order])
        
#    Chem.SanitizeMol(tmc_mol) # This is breaking things
    
    return tmc_mol