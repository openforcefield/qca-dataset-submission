import json

for file in ['set1-N_upto_40_hac_20.json', 'set1-S_upto_40_hac_20.json', 'set1-P_upto_40_hac_20.json']:
    with open(file, "r") as f:
        data = json.load(f)
    
    all_smiles = []
    for k, v in data.items():
        print(k, len(v))
        all_smiles.extend(v)
    
    with open(f"{file[:6]}-smiles.smi", mode="w") as f:
        # Joins list items with a newline, then writes the complete string
        f.write("\n".join(all_smiles) + "\n") 
