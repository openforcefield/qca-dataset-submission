import json
import pandas as pd
from openeye import oemolprop
from openeye import oechem
from collections import defaultdict
from openff.toolkit.topology import Molecule

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


new_dict1 = defaultdict(list)
new_dict2 = defaultdict(list)

with open('pubchem_NPS_search.json', "r") as f:
    data = json.load(f)

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

for k, v in data.items():
    nd_val = []
    for smi in v:
        if len(nd_val) == n_mols:
            break
        try:
            offmol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
            mol = offmol.to_openeye()
        except:
            continue
        
        if has_allowed_elements(mol, allowed_elements) and IsMoleculeInHeavyAtomCountRange(5, n_heavy, mol) and filter_obj(mol):
            nd_val.append(smi)

    new_dict1[k] = nd_val[:10]
    new_dict2[k] = nd_val[10:]


def subset_dict(d, key_list):
    return {k: d[k] for k in key_list if k in d}

set1_N = subset_dict(new_dict1, n_smarts)
set1_S = subset_dict(new_dict1, s_smarts)
set1_P = subset_dict(new_dict1, p_smarts)

set2_N = subset_dict(new_dict2, n_smarts)
set2_S = subset_dict(new_dict2, s_smarts)
set2_P = subset_dict(new_dict2, p_smarts)


# ---- write out ----
def dump_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=4)

dump_json(set1_N, f"set1-N_upto_{n_heavy}_hac_{n_mols}.json")
dump_json(set1_S, f"set1-S_upto_{n_heavy}_hac_{n_mols}.json")
dump_json(set1_P, f"set1-P_upto_{n_heavy}_hac_{n_mols}.json")

dump_json(set2_N, f"set2-N_upto_{n_heavy}_hac_{n_mols}.json")
dump_json(set2_S, f"set2-S_upto_{n_heavy}_hac_{n_mols}.json")
dump_json(set2_P, f"set2-P_upto_{n_heavy}_hac_{n_mols}.json")

# ---- WRITE JSON SETS ----

with open(f'pubchem_NPS_search_upto_{n_heavy}_hac_{int(n_mols/2)}_mols_set1.json', "w") as f:
    json.dump(new_dict1, f, indent=4)

with open(f'pubchem_NPS_search_upto_{n_heavy}_hac_{int(n_mols/2)}_mols_set2.json', "w") as f:
    json.dump(new_dict2, f, indent=4)


# ---- WRITE SMILES FILES ----

# --- SET 1 ---
all_smiles_set1 = []
for k, v in new_dict1.items():
    all_smiles_set1.extend(v)

with open(f"all_smiles_upto_{n_heavy}_and_{int(n_mols/2)}mols_set1.smi", "w") as f:
    f.write("\n".join(all_smiles_set1) + "\n")


# --- SET 2 ---
all_smiles_set2 = []
for k, v in new_dict2.items():
    all_smiles_set2.extend(v)

with open(f"all_smiles_upto_{n_heavy}_and_{int(n_mols/2)}mols_set2.smi", "w") as f:
    f.write("\n".join(all_smiles_set2) + "\n")
