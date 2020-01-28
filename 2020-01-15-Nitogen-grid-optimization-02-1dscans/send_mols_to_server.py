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
            if atom.IsInRing()==False:
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
                    }]
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

jobsDict={}


#mol_smiles=['CN(c1ccccc1)O','CNc1ccccc1','CN(C)c1ccccc1','c1ccc(cc1)N(O)O','c1ccc(cc1)N(c2ccccc2)c3ccccc3','CN(c1ccccc1)c2ccccc2','c1ccc(cc1)N(c2ccccc2)O','c1ccc(cc1)Nc2ccccc2','CN(c1ccccc1)S','c1ccc(cc1)N(S)S','c1ccc(cc1)N(c2ccccc2)S','c1ccc(cc1)N(S(=O)=O)S(=O)=O','c1ccc(cc1)NS(=O)=O','c1ccc(cc1)N(c2ccccc2)S(=O)=O','c1ccc(cc1)NO','CN(c1ccccc1)S(=O)=O','CN(c1ccccc1)S(=O)=O','c1ccc(cc1)N(O)S(=O)=O','N(S(=O)=O)(S(=O)=O)S(=O)=O','CN(C)S(=O)=O','CN(S(=O)=O)S(=O)=O','CNS(=O)=O','NS(=O)=O','CNS(=O)=O','N(O)(S(=O)=O)S(=O)=O','N(S(=O)=O)S(=O)=O','CN(c1ccccc1)S(=O)=O','c1ccc(cc1)N(c2ccccc2)S(=O)=O','CNC','CN','N','c1ccc(cc1)Nc2ccccc2','c1ccc(cc1)N(c2ccccc2)c3ccccc3','CN(c1ccccc1)c2ccccc2','CN(C)c1ccccc1','c1ccc(cc1)N','CNc1ccccc1','CN(C)C','N(O)O','NO','N(O)(O)O','c1ccc(cc1)NO','c1ccc(cc1)N(O)O','c1ccc(cc1)N(c2ccccc2)O','CNC','c1cnc(cn1)Nc2cnccn2','CNc1cnccn1','c1cnc(cn1)N','CN(c1cnccn1)c2cnccn2','c1cnc(cn1)N(c2cnccn2)c3cnccn3','c1cnc(cn1)NS(=O)=O','c1cnc(cn1)N(S(=O)=O)S(=O)=O','c1cnc(cn1)NS(=O)=O','c1cnc(cn1)N(c2cnccn2)S(=O)=O','CN(C)c1cnccn1','CNc1ccnnc1','c1cnncc1Nc2ccnnc2','c1cnncc1N(c2ccnnc2)c3ccnnc3','CN(c1ccncc1)c2ccncc2','CNc1ccncc1','CN(C)c1ccncc1','c1cnccc1N','c1cnccc1N(c2ccncc2)c3ccncc3','c1cnccc1N(c2ccncc2)S(=O)=O','c1cnccc1Nc2ccncc2','CN(c1cncnc1)c2cncnc2','c1c(cncn1)Nc2cncnc2','c1c(cncn1)N(c2cncnc2)c3cncnc3','CNc1cncnc1','c1c(cncn1)N','CN(C)c1cncnc1','CN(c1cncnc1)S(=O)=O','c1c(cncn1)N(S(=O)=O)S(=O)=O','CN(C)C=O','C(=O)N','CNC=O','C(=O)NC=O','CN(C=O)C=O','C(CN)C=O','CNCCC=O','CN(C)CCC=O','CN(C=O)(C=O)C=O','C(C=O)NCC=O','CN(CC=O)CC=O','C(C=O)N(CC=O)CC=O','c1([nH]nnn1)Nc2[nH]nnn2','CN(c1[nH]nnn1)c2[nH]nnn2','CNc1[nH]nnn1','CN(C)c1[nH]nnn1','c1([nH]nnn1)N','c1([nH]nnn1)N(c2[nH]nnn2)c3[nH]nnn3','CN(c1cnccn1)c2cnccn2','c1cnc(cn1)Nc2cnccn2','c1cnc(cn1)N(c2cnccn2)c3cnccn3','CN(C)c1cnccn1','CNc1cnccn1','c1cnc(cn1)N','c1cnncc1N','CNc1ccnnc1','CN(c1ccnnc1)c2ccnnc2','c1cnncc1Nc2ccnnc2','c1cnncc1N(c2ccnnc2)c3ccnnc3','CN(C)c1ccnnc1','c1c[nH]cc1N','CN(C)c1cc[nH]c1','CNc1cc[nH]c1','CN(c1cc[nH]c1)c2cc[nH]c2','c1c[nH]cc1Nc2cc[nH]c2','c1c[nH]cc1N(c2cc[nH]c2)c3cc[nH]c3','CNc1[nH]nnn1','CN(C)c1[nH]nnn1','c1([nH]nnn1)N','CN(c1[nH]nnn1)c2[nH]nnn2','c1([nH]nnn1)Nc2[nH]nnn2','c1([nH]nnn1)N(c2[nH]nnn2)c3[nH]nnn3','CNc1cnsn1']


mol_smiles=['CN(c1ccccc1)O','CNc1ccccc1','CN(C)c1ccccc1','c1ccc(cc1)N(O)O','c1ccc(cc1)N(c2ccccc2)c3ccccc3','CN(c1ccccc1)c2ccccc2','c1ccc(cc1)N(c2ccccc2)O','c1ccc(cc1)Nc2ccccc2','CN(c1ccccc1)S','c1ccc(cc1)N(S)S','c1ccc(cc1)N(c2ccccc2)S','CS(=O)(=O)N(c1ccccc1)S(=O)(=O)C','CS(=O)(=O)Nc1ccccc1','CS(=O)(=O)N(c1ccccc1)c2ccccc2','c1ccc(cc1)NO','CN(c1ccccc1)S(=O)(=O)C','CN(c1ccccc1)S(=O)(=O)C','CS(=O)(=O)N(c1ccccc1)O','CS(=O)(=O)N(S(=O)(=O)C)S(=O)(=O)C','CN(C)S(=O)(=O)C','CN(S(=O)(=O)C)S(=O)(=O)C','CNS(=O)(=O)C','CS(=O)(=O)N','CNS(=O)(=O)C','CS(=O)(=O)N(O)S(=O)(=O)C','CS(=O)(=O)NS(=O)(=O)C','CN(c1ccccc1)S(=O)(=O)C','CS(=O)(=O)N(c1ccccc1)c2ccccc2','CNC','CN','N','c1ccc(cc1)Nc2ccccc2','c1ccc(cc1)N(c2ccccc2)c3ccccc3','CN(c1ccccc1)c2ccccc2','CN(C)c1ccccc1','c1ccc(cc1)N','CNc1ccccc1','CN(C)C','N(O)O','NO','N(O)(O)O','c1ccc(cc1)NO','c1ccc(cc1)N(O)O','c1ccc(cc1)N(c2ccccc2)O','CNC','c1cnc(cn1)Nc2cnccn2','CNc1cnccn1','c1cnc(cn1)N','CN(c1cnccn1)c2cnccn2','c1cnc(cn1)N(c2cnccn2)c3cnccn3','CS(=O)(=O)Nc1cnccn1','CS(=O)(=O)N(c1cnccn1)S(=O)(=O)C','CS(=O)(=O)Nc1cnccn1','CS(=O)(=O)N(c1cnccn1)c2cnccn2','CN(C)c1cnccn1','CNc1ccnnc1','c1cnncc1Nc2ccnnc2','c1cnncc1N(c2ccnnc2)c3ccnnc3','CN(c1ccncc1)c2ccncc2','CNc1ccncc1','CN(C)c1ccncc1','c1cnccc1N','c1cnccc1N(c2ccncc2)c3ccncc3','CS(=O)(=O)N(c1ccncc1)c2ccncc2','c1cnccc1Nc2ccncc2','CN(c1cncnc1)c2cncnc2','c1c(cncn1)Nc2cncnc2','c1c(cncn1)N(c2cncnc2)c3cncnc3','CNc1cncnc1','c1c(cncn1)N','CN(C)c1cncnc1','CN(c1cncnc1)S(=O)(=O)C','CS(=O)(=O)N(c1cncnc1)S(=O)(=O)C','CN(C)C=O','C(=O)N','CNC=O','C(=O)NC=O','CN(C=O)C=O','C(CN)C=O','CNCCC=O','CN(C)CCC=O','C(=O)N(C=O)C=O','C(C=O)NCC=O','CN(CC=O)CC=O','C(C=O)N(CC=O)CC=O','c1([nH]nnn1)Nc2[nH]nnn2','CN(c1[nH]nnn1)c2[nH]nnn2','CNc1[nH]nnn1','CN(C)c1[nH]nnn1','c1([nH]nnn1)N','c1([nH]nnn1)N(c2[nH]nnn2)c3[nH]nnn3','CN(c1cnccn1)c2cnccn2','c1cnc(cn1)Nc2cnccn2','c1cnc(cn1)N(c2cnccn2)c3cnccn3','CN(C)c1cnccn1','CNc1cnccn1','c1cnc(cn1)N','c1cnncc1N','CNc1ccnnc1','CN(c1ccnnc1)c2ccnnc2','c1cnncc1Nc2ccnnc2','c1cnncc1N(c2ccnnc2)c3ccnnc3','CN(C)c1ccnnc1','c1c[nH]cc1N','CN(C)c1cc[nH]c1','CNc1cc[nH]c1','CN(c1cc[nH]c1)c2cc[nH]c2','c1c[nH]cc1Nc2cc[nH]c2','c1c[nH]cc1N(c2cc[nH]c2)c3cc[nH]c3','CNc1[nH]nnn1','CN(C)c1[nH]nnn1','c1([nH]nnn1)N','CN(c1[nH]nnn1)c2[nH]nnn2','c1([nH]nnn1)Nc2[nH]nnn2','c1([nH]nnn1)N(c2[nH]nnn2)c3[nH]nnn3','CNc1cnsn1']


added_mols=['CCN(CC)CC','CCNCC','CCN','CCCNCCC','CCCN(CCC)CCC','CCCN','C1CC1NC2CC2','C1CC1N(C2CC2)C3CC3','C1CC1N','C1CCC(C1)NC2CCCC2','C1CCC(C1)N(C2CCCC2)C3CCCC3','C1CCC(C1)N','C1CNCCC1NC2CCNCC2','C1CNCCC1N','C1CNCCC1N(C2CCNCC2)C3CCNCC3','N(S)S','NS','N(S)(S)S','N(F)F','N(F)(F)F','NF','C1CCC(CC1)NC2CCCCC2','C1CCC(CC1)N(C2CCCCC2)C3CCCCC3','C1CCC(CC1)N','N(Cl)Cl','NCl','N(Cl)(Cl)Cl','C1CCCC(CC1)NC2CCCCCC2','N(Br)Br','N(Br)(Br)Br','NBr','C(N(CS(=O)=O)CS(=O)=O)S(=O)=O','C(N)S(=O)=O','C1CCCC(CC1)N','C1CCCC(CC1)N(C2CCCCCC2)C3CCCCCC3','C(CS(=O)=O)NCCS(=O)=O','C(CS(=O)=O)N','C(NCS(=O)=O)S(=O)=O','C(CCS(=O)=O)CN','C(CN)CS(=O)=O','C(CCN)CCS(=O)=O','CN(C)C','C(CS(=O)=O)N(CCS(=O)=O)CCS(=O)=O','C(CN(CCCS(=O)=O)CCCS(=O)=O)CS(=O)=O','C(CCNCCCCCS(=O)=O)CCS(=O)=O','C(CCS(=O)=O)CNCCCCS(=O)=O','C(CCS(=O)=O)CN(CCCCS(=O)=O)CCCCS(=O)=O','C(CCN(CCCCCS(=O)=O)CCCCCS(=O)=O)CCS(=O)=O','CCCCN(CCCC)CCCC','CCCCNCCCC','CCCCN','CCCCCN(CCCCC)CCCCC','CCCCCNCCCCC','CCCCCN','c1ccc(cc1)CNCc2ccccc2','c1ccc(cc1)CCNCCc2ccccc2','c1ccc(cc1)CCCN(CCCc2ccccc2)CCCc3ccccc3','c1ccc(cc1)CCN(CCc2ccccc2)CCc3ccccc3','c1ccc(cc1)CN(Cc2ccccc2)Cc3ccccc3','c1ccc(cc1)CCCN','c1ccc(cc1)CCN','c1ccc(cc1)CN','c1ccc(cc1)CCCNCCCc2ccccc2']

error_type=[]
smile_errors=[]
for smiles in added_mols:
    #tmp_mol = oechem.OEMol()
    #ifs = oechem.oemolistream(f)
    #oechem.OEReadMolecule(ifs, tmp_mol)
    #mol = oechem.OEMol(tmp_mol)
    #status = omega(mol)

    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)

    # Initialize Omega
    omega = oeomega.OEOmega()

    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(True)
    omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
    omega.SetStrictStereo(True)
    omega.SetStrictAtomTypes(True)
    omega.SetIncludeInput(False) # don't include input
    status = omega(mol)

    try:
        nitrogen, nitrogen_i = find_nitrogen(mol)
        ptl_mol = make_ptl_mol(mol)
        subDict = send_qm_job(ptl_mol, nitrogen, nitrogen_i, mol)
    except Exception as e:
        smile_errors.append(smiles)
        error_type.append(e)

        ofs = oechem.oemolostream(smiles+".mol2")
        oechem.OEWriteMolecule(ofs, mol)

        print("this smile caused issues:" + smiles)
    try:
        jobsDict[subDict[1]]=subDict[0]
    except:
        pass








import json
with open('additional_1d_scan_jobs.json', 'w') as fp:
        json.dump(jobsDict, fp, indent=2, sort_keys=True)



print("This is the final list of smiles that did not get incldued in the .json")
print(smile_errors)
print(error_type)
print(len(added_mols))
