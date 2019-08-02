#imports
import time
import pprint
import re
import numpy as np

from openeye import oechem
from openeye import oeomega
import qcportal as ptl


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
    coord_str = '\n'.join(
        (f"{oechem.OEGetAtomicSymbol(atom.GetAtomicNum())}   "
                   f"{'   '.join(str(c) for c in coords[atom.GetIdx()])}")
                  for atom in mol.GetAtoms())
    print(coord_str)
    conn = np.array([[bond.GetBgnIdx(), bond.GetEndIdx(), bond.GetOrder()] for bond
            in mol.GetBonds()])
    return ptl.Molecule.from_data(coord_str, connectivity=conn)

def send_qm_job(ptl_mol, nitrogen, nitrogen_i):
    """Sends a job to the QM Client - returns a submitted object"""
    indices = [nitrogen_i] + [nbor.GetIdx() for nbor in list(nitrogen.GetAtoms())]
    print(f"indices: {indices}")
    keywords = ptl.models.KeywordSet(values={"scf_properties":["wiberg_lowdin_indices"]})
    keywords_id = client.add_keywords([keywords])[0]
    service = ptl.models.GridOptimizationInput(**{
            "keywords": {
                "preoptimization": True,
                "scans": [{
                    "type": "dihedral",
                    "indices": indices,
                    "steps": [-40, -36, -32, -28, -24, -20, -16, -12, -8, -4, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40],
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
            "initial_molecule": ptl_mol,
        })
    submitted = client.add_service([service])
    print("Job submitted")
    return submitted



#This is where we submit the job.
#The molecule we ran this example with stored as a smile string in the tiny.smi file.
#This should be adapted for the directory "Molecules_to_run" for the .sdf files

tmp_mol = oechem.OEMol()
ifs = oechem.oemolistream("tiny.smi")
first = True

results = [] # {"molecule": <OEMol>, "nitrogen": <OEAtom>, "nitrogen_i": <int>,
             #  "ptl_molecule": <PtlMol>, submitted": <submitted object>,
             #  "res": <result object> from QCPortal}

while oechem.OEReadMolecule(ifs, tmp_mol):
    # Separate outputs by a line
    if first: first = False
    else: print()

    mol = oechem.OEMol(tmp_mol)
    status = omega(mol)
    nitrogen, nitrogen_i = find_nitrogen(mol)
    print(f"Nitrogen found at index {nitrogen_i}")
    print(f"Generating conformer: {'done' if status else 'failed'}")

    print(f"Nitrogen is at: {mol.GetCoords()[nitrogen_i]}")

    ptl_mol = make_ptl_mol(mol)
    sub = send_qm_job(ptl_mol, nitrogen, nitrogen_i)

    results.append({
        "molecule": mol,
        "nitrogen": nitrogen,
        "nitrogen_i": nitrogen_i,
        "ptl_molecule": ptl_mol,
        "submitted": sub,
    })
    break

