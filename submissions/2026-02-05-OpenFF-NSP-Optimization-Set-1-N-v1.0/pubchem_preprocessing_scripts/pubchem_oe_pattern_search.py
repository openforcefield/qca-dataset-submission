import pandas as pd
import json
from openeye import oechem
from collections import defaultdict
from openeye import oegraphsim
from openeye import oemolprop

def IsBetween(min, max, val):
    if min <= val <= max:
        return True
    return False


def IsMoleculeInHeavyAtomCountRange(min, max, mol):
    count = oechem.OECount(mol, oechem.OEIsHeavy())
    return IsBetween(min, max, count)

def has_allowed_elements(mol, allowed_elements):
    elements = set()
    for atom in mol.GetAtoms():
        # Get atomic number (int)
        atomic_num = atom.GetAtomicNum()
        # Get atomic symbol (str)
        atomic_symbol = oechem.OEGetAtomicSymbol(atomic_num)
        elements.add(atomic_symbol)
    if elements <= allowed_elements:
        return True
    else:
        return False


def get_tree_fp(oemol):
    fp = oegraphsim.OEFingerPrint()
    oegraphsim.OEMakeFP(fp, oemol, oegraphsim.OEFPType_Tree)
    return fp


n_df = pd.read_csv('nitrogen_summary_updated.md', delimiter=',')
n_df = n_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

s_df = pd.read_csv('sulfur_summary_updated.md', delimiter=',')
s_df = s_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

p_df = pd.read_csv('phosphorous_summary_updated.md', delimiter=',')
p_df = p_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

n_smarts = n_df['SMARTS'].to_list()
s_smarts = s_df['SMARTS'].to_list()
p_smarts = p_df['SMARTS'].to_list()

all_smarts = n_smarts + s_smarts + p_smarts
fname = 'file_name'
ifs = oechem.oemolistream(f'/dfs6/pub/pbehara/check_chembl_parameter_coverage/pubchem/Compounds/{fname}')

allowed_elements = set(['H','C','N','O','F','P','S','Cl','Br','I'])
n_heavy = 40
n_mols = 20
# 1. Create an input stream from your custom filter file
filter_file_path = "FILTER.OE"
filter_istream = oechem.oeifstream(filter_file_path)

# Check if the file was successfully opened
if not filter_istream.IsValid():
    print(f"Error: Could not open custom filter file {filter_file_path}")
    exit()

# 2. Initialize the OEFilter object with the custom filter stream
filter_obj = oemolprop.OEFilter(filter_istream)
print(filter_obj.GetTypeCheck())
filter_obj.SetTypeCheck(True)

mol = oechem.OEGraphMol()
smarts_dict = defaultdict(list)
while oechem.OEReadMolecule(ifs, mol):
    added_to_pattern = False
    oechem.OEDeleteEverythingExceptTheFirstLargestComponent(mol)
    if not (has_allowed_elements(mol, allowed_elements) and IsMoleculeInHeavyAtomCountRange(5, n_heavy, mol) and filter_obj(mol)): 
        continue
    mol_name = mol.GetTitle()
    for pattern in all_smarts:
        # skip already added molecules
        if added_to_pattern:
            break
        if len(smarts_dict[pattern]) < n_mols:
            # create a substructure search object
            ss = oechem.OESubSearch(pattern)
            oechem.OEPrepareSearch(mol, ss)
        
            if ss.SingleMatch(mol):
                smiles = oechem.OEMolToSmiles(mol)
                if len(smarts_dict[pattern]) == 0:
                    smarts_dict[pattern].append((f"{smiles} Pubchem_CID_{mol_name}", get_tree_fp(mol)))
                    added_to_pattern = True
                else:
                    fpmol = get_tree_fp(mol)
                    sim_with_rest = [oegraphsim.OETanimoto(fpmol, fp) for _,fp in smarts_dict[pattern]]
                    # If tanimoto sim is less than 0.4 then add molecule to the list
                    if all(x < 0.4 for x in sim_with_rest):
                        smarts_dict[pattern].append((f"{smiles} Pubchem_CID_{mol_name}", fpmol))
                        added_to_pattern = True
        else:
            all_smarts.remove(pattern)
            break
print(all_smarts)
print('###########################\n')
print(set(all_smarts) - set(smarts_dict.keys()))

for key, value in smarts_dict.items():
    smarts_dict[key] = [x[0] for x in value] 
with open(f'{fname[:-7]}_smarts_dict.json', "w") as f:
    json.dump(smarts_dict, f)
