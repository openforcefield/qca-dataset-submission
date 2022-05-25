#!/usr/bin/env python

"""
Utility functions for SMARTS-based querying of GeoOpt datasets on QCArchive
Created by:
- @pavankum
"""

from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from rdkit.Chem import rdMolTransforms
from simtk.unit import bohrs
from simtk import openmm, unit
import numpy as np
import pickle
from openeye import oequacpac



smilesdict={'t130b' : ['C1=CC=C(C=C1)S(=O)(=O)NN', 'CS(=O)(=O)NN(CC)S(=O)(=O)C', 'CCOC1=CC=C(C=C1)S(=O)(=O)NN', 'CN(C)C1=CC=CC2=C1C=CC=C2S(=O)(=O)NN', 'CS(=O)(=O)NN(CC)S(=O)(=O)C'],
't132b' : ['CC(C)NC(=O)C1=CC=C(C=C1)CNNC', 'CCNNC', 'CCNNCS(=O)(=O)C', 'C1=CC=C(C=C1)CCNN', 'c1ccc(cc1)C(=O)CNN'],
't133':['C1C[N+](CN1)([N+]2(CCNC2)[O-])[O-]',
'C1C[NH+]([NH+]1[O-])[O-]',
'C1CC[N+](CC1)([N+]2(CCCCC2)[O-])[O-]', 'CC(=O)C(C)(C)[NH+]([N+](C)(CCNC(C)(C)C(=O)CN1C=CN=C1[N+](=O)[O-])[O-])[O-]', 'CC(=O)C(C)(C)NCC[N+](C)([NH+](C(C)(C)C(=O)C)[O-])[O-]'],
't133a':['CN(CCC(=O)[O-])[N+](C)(C)C', 'C[N+](C)(C)NCCC(=O)[O-]', 'C[N+](C)(C=C)NCCC(=O)[O-]', 'CN(C)[N+](C)(C)C', 'C[N+](C)(C)NCCS(=O)(=O)[O-]'],
't142b':['CNSN(C)C(=O)ON', 'CN1C(=O)C=CS1', 'CCCCCCCCN1C(=O)C=CS1', 'C1=CC=C2C(=C1)C(=O)NS2', 'C1CCC(CC1)NSC2=NC3=CC=CC=C3S2'],
't142d':['C1=N[S+](NC1=N)[O-]', 'CC(=O)O[S+](C)C', 'CCN=C1C(=N[S+](N1)[O-])N', 'C[S+](C)OS(=O)(=O)c1ccccc1', 'c1ccc(cc1)[S+](c2ccccc2)O'],
't142e':['CN1C(=O)c2cc(c(cc2[S+]1[O-])OC)OC', 'C1(=N[S+](NC1=N)[O-])N', 'CN=C1C(=NCCS)N[S+](N1)[O-]', 'CSCC(N=C1C(=N[S+](N1)[O-])N)O', 'c1cc(oc1)CN=C2C(=N[S+](N2)[O-])N'],
't142f':['c1ccc(cc1)[S+](c2ccccc2)S(=O)(=O)O', 'C[SH+]S(=O)(=O)c1ccccc1', 'COc1ccccc1S(=O)(=O)[S+](C)S(=O)(=O)C', 'CC[S+](CC)S(=O)(=O)c1ccc[nH]1', 'CC[S+](C)S(=O)(=O)O'],
't142c':['CC(C(=O)OC)SP(=S)(OC)OC', 'COP(=S)(OC)SCN1C(=O)c2ccccc2C1=O', 'CCOP(=S)(OCC)S', 'CCSP(=S)(OC(C)C)OC(C)C', 'COP(=S)(OC)Sc1ccccc1'],
't143d':['c1ccc(cc1)C2=NS(=O)ON2', 'C1=NS(=O)ON1', 'Cc1ccs(=O)n1', 'CC(C(=O)C1=NS(=O)ON1)O', 'CCc1cs(=O)nc1N'],
't143e':['CN(C)S(=O)N(C)c1ccccc1', 'CCOS(=O)NC', 'CCOS(=O)Nc1ccccc1', 'CN(C)S(=O)N(C)C', 'NS(=O)O'],
't143f':['C[N+]1(CCCC1)S(=O)O', 'C[N+](C)(C)S(=O)O', 'C[N+](C)(c1ccccc1)S(=O)O', '[NH3+]S(=O)O', 'CC(=O)[NH2+]S(=O)O'],
't122d':['CC1=CC=C(C=C1)S(=O)C#C', 'CC(C)(C)OC#CS(=O)C1=CC=CC=C1', 'C1=CC=C(C=C1)C2=CN=[C-]S2=O', 'CCCC#CS(=O)C1=CC=C(C=C1)C', 'c1ccc(cc1)S(=O)C#N'],
't122e':['c1ccc(cc1)S(=O)c2ccccc2', 'c1ccc(cc1)S(=O)c2ccc3c(c2)[nH]c(n3)N', 'c1ccc2c(c1)Sc3ccccc3S2=O', 'Cc1ccc2c(c1)N(C(=O)c3ccccc3S2=O)C', 'c1cc(cc(c1)S(=O)c2cc(cc(c2)N)N)N'],
't122a': ['CC#CS(=O)(=O)C1=CC=CC=C1', 'C1=CC=C(C=C1)S(=O)(=O)C#N', 'CC1=CC=C(C=C1)S(=O)(=O)C#CC2=CC=CC=C2', 'CC#CS(=O)(=O)c1cc(cc(c1)N)N', 'c1c(cc(cc1O)S(=O)(=O)C#N)N'],
't116b':['C=C[S+]1CCCC1', 'C[S+](C)C=CC(=O)[O-]', 'CCC(C=C[S+](C)C)N', 'C[S+](C)C=CC(=O)Nc1ccccc1', 'C=C[S+]1CC(C(C1)N)N'],
't116c':['CN1C2C[S+]3CCCC3C2N(C1=O)C', 'C[S+](C)CCC(C(=O)O)N', 'C[S+](C)CCC(=O)[O-]', 'C[S+](C)CCC(=O)Nc1ccccc1', 'C[S+](C)CCc1ccccc1'],
't133':['C1C[N+](CN1)([N+]2(CCNC2)[O-])[O-]', 'C1C[NH+]([NH+]1[O-])[O-]', 'C1CC[N+](CC1)([N+]2(CCCCC2)[O-])[O-]', 'CC[N+](C)([NH+](C(C)(C)C(=O)C)[O-])[O-]', 'CC[NH+]([N+](C)(CC)[O-])[O-]'],
't33a':['CN(CCC(=O)[O-])[N+](C)(C)C', 'C[N+](C)(C)NCCC(=O)[O-]', 'C[N+](C)(C=C)NCCC(=O)[O-]', 'CN(C)[N+](C)(C)C', 'C[N+](C)(C)NCCS(=O)(=O)[O-]'],
't74a':['CC(=O)[N-]NC', 'CCC(=O)[N-]C(=O)C', 'CCC(=O)[N-]C', 'CS(=O)(=O)[N-]C(=O)c1ccc(cc1)N', 'CC[N-]C(=O)c1ccncc1'],
't157':['NS(=O)(=O)O', 'C(CS(=O)(=O)O)N', 'CS(=O)(=O)OCCCCOS(=O)(=O)C', 'OS(=O)(=O)O', 'CC1=CC(=O)NS(=O)(=O)O1'],
't157a': ['OS(=O)O', 'C(CS(=O)O)N', 'C(=N)(N)S(=O)O', 'CCOS(=O)OCC', 'C1COS(=O)O1']
}


count_dict={}
match_check={}
for key, item in smilesdict.items():
    count_dict[key]=0
    match_check[key]=0
def get_gopt_matching_improper(smiles, ff, hess=None):
    """Get all Optimization entries from the given QCArchive instance datasets that
    match the given SMARTS pattern.

    Parameters
    ----------
    smarts : str
        SMARTS to query optimization datasets against optimzation datasets

    Returns
    -------
    opt_recs:
        List of optimization records matching the given smarts
    ind_set
        Dict of a list of matching indices per record id
    wbos
        Dict of wbos corresponding to the improper 2nd and 4th atoms

    Examples
    --------
    Get back TDEntries corresponding to a C-C bond torsiondrive:

    """
    from collections import defaultdict

    # first, we want to grab all datasets into memory
    opt_recs = []
    ind_set = defaultdict(list)
    dihs_set = []
    wbos = defaultdict(list)
    failed_smiles = []
    imp_dict=defaultdict(list)
    imp_dict_mm=defaultdict(list)
    plot_dict={}
    pccharge=defaultdict(list)
    imp_params=defaultdict(list)
    angle_params=defaultdict(list)
    torsion_params=defaultdict(list)
    mols_dict=defaultdict(list)
    molecules=0
    torsions_dict={}
    allmols=[]
    dictmols=[]
    for mol_smiles in smiles:
        try:
            #pavan, this is where I am creating the Molecule from smiles
            offmol = Molecule.from_smiles(mol_smiles, allow_undefined_stereo=True)
        except:
            failed_smiles.append(mol_smiles)
            print("Failed: ", mol_smiles)
            continue

        # apply SMARTS to each Optimization object
        #matching_indices = offmol.chemical_environment_matches(smarts, allow_cosmetic_atttributes=True)
        forcefield = ForceField(ff, allow_cosmetic_attributes=True)
        topology = Topology.from_molecules(offmol)
        molecule_force_list = forcefield.label_molecules(topology)
        imp_parameters={}
        print(molecule_force_list)
        params=[]
        molecules+=1
        allmols.append(mol_smiles)
        for mol_idx, mol_forces in enumerate(molecule_force_list):
            # print(f'Forces for molecule {mol_idx}')
            for force_tag, force_dict in mol_forces.items():
                if force_tag == "ProperTorsions":
                    for (atom_indices, parameter) in force_dict.items():
                        #params.append(parameter.id)
                        params.append(parameter.id)
                        if parameter.id in count_dict.keys():
                            print('match')
                            match_check[parameter.id]+=1
                        if parameter.id in count_dict.keys() and mol_smiles in smilesdict[parameter.id]:
                            count_dict[parameter.id]+=1
                            print(list(atom_indices))
                            #mols_dict[mol_smiles].append(list(atom_indices))
                            dictmols.append(mol_smiles)
                            if parameter.id not in torsions_dict.keys():
                                torsions_dict[parameter.id] = {}
                            #pavan - this is where I am storing the mapped smiles and the dihedral indices for the dictionary with the torsion drive data
                            if mol_smiles not in torsions_dict[parameter.id].keys():
                                torsions_dict[parameter.id][offmol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)] = []
                            torsions_dict[parameter.id][offmol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)].append(list(atom_indices))
                        if parameter.id[0:3] == 't73':
                            torsion_params[parameter.id].append(mol_smiles)
                        if parameter.id[0:3] == 't74':
                            torsion_params[parameter.id].append(mol_smiles)
                        if parameter.id[0:4] == 't116':
                            torsion_params[parameter.id].append(mol_smiles)
                        if parameter.id[0:4] == 't122':
                            torsion_params[parameter.id].append(mol_smiles)
                        if parameter.id[0:4] == 't130':
                            torsion_params[parameter.id].append(mol_smiles)
                        if parameter.id[0:4] == 't132':
                            torsion_params[parameter.id].append(mol_smiles)
                        if parameter.id[0:4] == 't133':
                            torsion_params[parameter.id].append(mol_smiles)
                        if parameter.id[0:4] == 't142':
                            torsion_params[parameter.id].append(mol_smiles)
                        if parameter.id[0:4] == 't143':
                            torsion_params[parameter.id].append(mol_smiles)
        print(set(params))
        print('complete')
        with open('torsions_check_new_dataset_8.pickle', 'wb') as handle:
            pickle.dump(torsion_params, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print(mols_dict)
        #torsions_dict
        #pavan
        #this is where I am generating the dictionary of data
        #for the torsion drive dataset
        with open('torsion_drive_data.pickle', 'wb') as handle:
            pickle.dump(torsions_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(count_dict)
    count=0
    for key, item in smilesdict.items():
        for x in item:
            count+=1
    print(count)
    print(molecules)
    #print(failed_smiles)
    print(len(allmols))
    print(len(dictmols))
    #print(torsions_dict)
    for key, item in torsions_dict.items():
        print(key)
        print(len(item))
    for x in allmols:
        if x not in dictmols:
            print(x)
            for key, item in smilesdict.items():
                if x in item:
                    print(key)

    print(match_check)
    return



def show_oemol_struc(oemol, torsions=False, atom_indices=[], width=500, height=300):
    from IPython.display import Image
    from openeye import oechem, oedepict

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

    class BondInTorsion(oechem.OEUnaryBondPred):
        def __call__(self, bond):
            return (bond.GetBgn().GetIdx() in atom_indices) and (
                bond.GetEnd().GetIdx() in atom_indices
            )

    class CentralBondInTorsion(oechem.OEUnaryBondPred):
        def __call__(self, bond):
            return (bond.GetBgn().GetIdx() in atom_indices[1:3]) and (
                bond.GetEnd().GetIdx() in atom_indices[1:3]
            )

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())

    oedepict.OEPrepareDepiction(oemol)
    img = oedepict.OEImage(width, height)
    display = oedepict.OE2DMolDisplay(oemol, opts)
    if torsions:
        atoms = oemol.GetAtoms(AtomInTorsion())
        bonds = oemol.GetBonds(NoBond())
        abset = oechem.OEAtomBondSet(atoms, bonds)
        oedepict.OEAddHighlighting(
            display,
            oechem.OEColor(oechem.OEYellow),
            oedepict.OEHighlightStyle_BallAndStick,
            abset,
        )

    oedepict.OERenderMolecule(img, display)
    png = oedepict.OEWriteImageToString("png", img)
    return Image(png)


def calc_improper(atom0, atom1, atom2, atom3, translate=False):
    """
    Calculate the improper dihedral angle of a set of given four atoms.

    Parameters
    ----------
    atom0 : numpy array
        CENTRAL atom coordinates
    atom1 : numpy array
        outer atom coordinates
    atom2 : numpy array
        outer atom coordinates
    atom3 : numpy array
        outer atom coordinates
    translate : bool
        True to translate central atom to origin, False to keep as is.
        This should not affect the results of the calculation.

    Returns
    -------
    float
        Angle in degrees.
    """
    if translate:
        atom1 = atom1 - atom0
        atom2 = atom2 - atom0
        atom3 = atom3 - atom0
        atom0 = atom0 - atom0 # central must be moved last

    # calculate vectors
    v0 = atom0-atom1
    v1 = atom2-atom1
    v2 = atom2-atom3
    w1 = np.cross(v0, v1)
    w2 = np.cross(v1, v2)
    angle = angle_between(w1,w2) # this angle should be in range [0,90]

    # compute distance from plane to central atom
    # eq 6 from http://mathworld.wolfram.com/Point-PlaneDistance.html
    # here I'm using atom1 for (x,y,z), but could also use atom2 or atom3
    numer = w2[0]*(atom0[0]-atom1[0]) + w2[1]*(atom0[1]-atom1[1]) + w2[2]*(atom0[2]-atom1[2])
    denom = np.sqrt(w2[0]**2 + w2[1]**2 + w2[2]**2)
    dist = numer/denom
    # set reference so that if central atom is above plane, angle -> [90,180]
    if dist > 0:
        #angle = 180-angle
        angle = -1*angle

    return angle

def angle_between(v1, v2):
    """
    Calculate the angle in degrees between vectors 'v1' and 'v2'.
    Modified from: https://tinyurl.com/yb89sstz

    Parameters
    ----------
    v1 : tuple, list, or numpy array
    v2 : tuple, list, or numpy array

    Returns
    -------
    float
        Angle in degrees.

    """
    v1_u = v1/np.linalg.norm(v1)
    v2_u = v2/np.linalg.norm(v2)
    return np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))





