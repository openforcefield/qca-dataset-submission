#imports
import time
import pprint
import re
import numpy as np

from openeye import oechem
from openeye import oeomega
import qcportal as ptl
import cmiles
import qcelemental as qcel

# Custom exception for the case when there is no nitrogen
class NoNitrogenException(Exception): pass

#identifies the invertible nitrogen that the grid optimization will occur around
def find_nitrogen(mol):
    """Returns the trivalent nitrogen atom in a molecule"""
    for atom in mol.GetAtoms():
        if oechem.OEIsInvertibleNitrogen()(atom):
            return atom, atom.GetIdx()
    raise NoNitrogenException()

# Initialize Omega
omega = oeomega.OEOmega()

omega.SetMaxConfs(1)
omega.SetIncludeInput(True)
omega.SetCanonOrder(True)
omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
omega.SetStrictStereo(True)
omega.SetStrictAtomTypes(True)
omega.SetIncludeInput(False) # don't include input


client = ptl.FractalClient("https://localhost:7777/", verify=False)


def make_ptl_mol(oemol):
    """Builds a QCPortal Molecule from an OpenEye molecule"""
    coords = oemol.GetCoords()
    symbols_list = [oechem.OEGetAtomicSymbol(atom.GetAtomicNum()) for atom in mol.GetAtoms()]

    #convert to bohr
    print(coords)
    for key, item in coords.items():
        coords[key] = (item[0]*qcel.constants.conversion_factor('angstrom', 'Bohr'), item[1]*qcel.constants.conversion_factor('angstrom', 'Bohr'), item[2]*qcel.constants.conversion_factor('angstrom', 'Bohr'))


    coord_list = [c for atom in mol.GetAtoms() for c in coords[atom.GetIdx()] ]
    conn_list = np.array([[bond.GetBgnIdx(),
                           bond.GetEndIdx(),
                           bond.GetOrder()] for bond
            in mol.GetBonds()])
    ptl_mol = ptl.Molecule.from_data(
        {'geometry':coord_list,
        'symbols':symbols_list,
        'connectivity':conn_list})

    return ptl_mol

def send_qm_job(ptl_mol, nitrogen, nitrogen_i,  mol):
    """Sends a job to the QM Client - returns a submitted object"""
    indices = [nitrogen_i] + [nbor.GetIdx() for nbor in list(nitrogen.GetAtoms())]

    AtomsAroundNit = list(nitrogen.GetAtoms())
    print(AtomsAroundNit)
    for atom in AtomsAroundNit:
        print(atom.GetAtomicNum())
    #make a list of the valence indices for the restrained optimization using the
    #newlist = sorted(AtomsAroundNit, key=lambda x: x.GetAtomicNum(), reverse=True)
    AtomsAroundNit.sort(key=lambda x: x.GetAtomicNum(), reverse=True)
    print(AtomsAroundNit)


    for atom in AtomsAroundNit:
        print(atom.GetIdx())

    try:
        valenceIdx=[AtomsAroundNit[0].GetIdx(), nitrogen_i, AtomsAroundNit[1].GetIdx()]
    except:
        pass

    #print(f"indices: {indices}")
    keywords = ptl.models.KeywordSet(values={"scf_properties":["wiberg_lowdin_indices"]})
    try:
        #keywords_id = (client.add_keywords([keywords])[0])

        keywords_id = str(client.add_keywords([keywords])[0])



        smiles=cmiles.utils.mol_to_smiles(mol, mapped=False, explicit_hydrogen=False)
        mol_id = cmiles.get_molecule_ids(smiles, toolkit='openeye', strict=False)

        connectivity=np.array(ptl_mol.connectivity).tolist()
        geometry=np.array([[ptl_mol.geometry]]).ravel().tolist()
        symbols=np.array([[ptl_mol.symbols]]).ravel().tolist()
        jsonDict={
                "cmiles_ids":mol_id,
                "keywords": {
                    "preoptimization": True,
                    "scans": [{
                        "type": "dihedral",
                        "indices": list(indices),
                        "steps": [-52 ,-48,-44,-40, -36, -32, -28, -24, -20, -16, -12, -8, -4, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52],
                        "step_type": "absolute"
                    }
                   , {
                "type": "angle",
                "indices": list(valenceIdx),
                "steps": [100, 105, 110, 115, 120, 125, 130, 135, 140],
                "step_type": "absolute"}
                    ]
                },
                "optimization_spec": {
                    "program": "geometric",
                    "keywords": {
                        "coordsys": "tric",
                    }
                },
                "qc_spec": {
                    "driver": "gradient",
                    "method": "mp2",
                    "basis": "def2-SV(P)",
                    "keywords": keywords_id,
                    "program": "psi4",
                },
                "initial_molecule":{
                    "geometry":geometry,
                    "symbols":symbols,
                    "connectivity":connectivity
                    }}
        return jsonDict, smiles

    except:
        pass
    return



#This is where we submit the job.
#The molecule we ran this example with stored as a smile string in the tiny.smi file.
#This should be adapted for the directory "Molecules_to_run" for the .sdf files

first = True

results = [] # {"molecule": <OEMol>, "nitrogen": <OEAtom>, "nitrogen_i": <int>,
             #  "ptl_molecule": <PtlMol>, submitted": <submitted object>,
             #  "res": <result object> from QCPortal}

import glob
file_list = glob.glob('./Molecules_to_run/*.*')
jobsDict={}
for f in file_list:
    tmp_mol = oechem.OEMol()
    ifs = oechem.oemolistream(f)
    oechem.OEReadMolecule(ifs, tmp_mol)
    mol = oechem.OEMol(tmp_mol)
    status = omega(mol)
    nitrogen, nitrogen_i = find_nitrogen(mol)
    ptl_mol = make_ptl_mol(mol)
    subDict = send_qm_job(ptl_mol, nitrogen, nitrogen_i, mol)
    try:
        jobsDict[subDict[1]]=subDict[0]
    except:
        pass
import json
with open('nitrogen_Jobs_2dscans.json', 'w') as fp:
        json.dump(jobsDict, fp, indent=2, sort_keys=True)
