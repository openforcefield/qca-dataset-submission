#!/usr/bin/env python

"""
Generate a PDF of all small molecules in the JSON dataset.
"""

import gzip
import json
import cmiles
from openeye import oechem

# Highlight element of interest
#subs = oechem.OESubSearch("[#6]") # carbon
subs = None

# Read the input SMILES
def read_molecules(input_json):
    """
    Parameters
    ----------
    input_json: str,
        JSON file name to the output json of generate.py (prepared as if for an OptimizationDataset)
        The data in the json file should be a list of {'initial_molecules': [..], 'cmiles_identifiers':{}}.

    Returns
    -------
    molecules_list_dict: dict
        The dictionary maps the index of a molecule to a Molecule object. e.g.
        {
            index1: [Molecule_json1a, Molecule_json1b, ..],
            index2: [Molecule_json2a, Molecule_json2b, ..],
        }

    """
    if input_json.endswith(".tar") or input_json.endswith(".tar.gz"):

        extract_file = input_json.replace(".gz", "").replace(".tar", ".json")
        with tarfile.open(input_json, 'r') as infile:

            molecule_data_list = json.load(infile.extractfile(extract_file))

    if input_json.endswith(".gz"):

        import gzip
        with gzip.open(input_json, 'r') as infile:

            molecule_data_list = json.loads(infile.read().decode('utf-8'))

    else:
        with open(input_json) as infile:
            molecule_data_list = json.load(infile)

    print(f'{len(molecule_data_list)} molecules read')
    return molecule_data_list

json_molecules = read_molecules('selected_torsions.json')

# Generate a PDF of all molecules in the set
pdf_filename = 'selected_torsions.pdf'

from openeye import oedepict
itf = oechem.OEInterface()
PageByPage = True
suppress_h = True
rows = 10
cols = 6
ropts = oedepict.OEReportOptions(rows, cols)
ropts.SetHeaderHeight(25)
ropts.SetFooterHeight(25)
ropts.SetCellGap(2)
ropts.SetPageMargins(10)
report = oedepict.OEReport(ropts)
cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
opts = oedepict.OE2DMolDisplayOptions(cellwidth, cellheight, oedepict.OEScale_Default * 0.5)
opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_On, 1.0)
opts.SetDefaultBondPen(pen)
oedepict.OESetup2DMolDisplayOptions(opts, itf)
for json_molecule in json_molecules.values():
    # Create oemol
    oemol = cmiles.utils.load_molecule(json_molecule['initial_molecules'][0])

    # Get atom indices
    atom_indices = json_molecule['atom_indices'][0]

    # Render molecule
    cell = report.NewCell()
    mol = oechem.OEMol(oemol)
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
    class BondInTorsion(oechem.OEUnaryBondPred):
        def __call__(self, bond):
            return (bond.GetBgn().GetIdx() in atom_indices) and (bond.GetEnd().GetIdx() in atom_indices)
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
    oedepict.OEAddHighlighting(disp, oechem.OEColor(oechem.OEOrange), oedepict.OEHighlightStyle_BallAndStick, abset)

    oedepict.OERenderMolecule(cell, disp)
    #oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

oedepict.OEWriteReport(pdf_filename, report)
