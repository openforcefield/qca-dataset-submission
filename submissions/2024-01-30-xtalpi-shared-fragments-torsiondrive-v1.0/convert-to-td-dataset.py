import pathlib
import re

import click
import tqdm
import numpy as np

def parse_opt_mol2(filename: str):
    from rdkit import Chem
    from openff.toolkit import Molecule
    
    rdmol = Chem.MolFromMol2File(str(filename), sanitize=False)
    if rdmol is None:
        raise ValueError(f"Could not read {filename}")

    filename = pathlib.Path(filename)
    
    entry = {}
    
    header_pattern = "@<.+>"
    with filename.open("r") as f:
        contents = f.read()
    sections = [
        section.split("\n")
        for section in re.split(header_pattern, contents)
        if section
    ]
    
    _to_array = lambda x: np.array(list(map(float, x)))
    
    # parse sections manually
    IGNORE = ["MOLECULE", "BOND"]
    COORDINATES = ["ATOM", "INPUT_MOLECULE"]
    FLOAT_TYPES = ["CORE_TIME", "ENERGY"]
    STRING_TYPES = ["CA_TYPE", "FINGERPRINT"]
    for section in sections:
        header = section[0]
        if header in IGNORE:
            continue
        lines = [x for x in section[1:] if x]
        if header in COORDINATES:
            coord_col = [
                x for line in lines
                for x in line.strip().split()[2:5]
            ]
            coords_ = _to_array(coord_col)
            entry[header] = coords_
        # we have to manually assign formal charges
        elif header == "UNITY_ATOM_ATTR":
            while lines:
                atom_index = int(lines[0].strip().split()[0]) - 1
                field, charge_ = lines[1].strip().split()
                if field == "charge":
                    print(atom_index, charge_)
                    rdatom = rdmol.GetAtomWithIdx(atom_index)
                    rdatom.SetFormalCharge(int(float(charge_)))
                lines = lines[2:]
        elif header in FLOAT_TYPES:
            assert len(lines) == 1
            entry[header] = float(lines[0].strip())
        elif header in STRING_TYPES:
            assert len(lines) == 1
            entry[header] = lines[0].strip()
        elif header == "FROZEN":
            values = list(map(int, lines[0].strip().split()))
            entry["frozen_atoms"] = values[:-1]
            entry["frozen_0"] = values[0]
            entry["frozen_1"] = values[1]
            entry["frozen_2"] = values[2]
            entry["frozen_3"] = values[3]
            entry["angle"] = values[-1]

    Chem.SanitizeMol(rdmol)
    offmol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    atomic_numbers = sorted(set([a.atomic_number for a in offmol.atoms]))

    entry.update(
        {
            "filename": filename.stem,
            "directory": filename.parent.name,
            "parent": filename.parent.parent.name,
            "smiles": Chem.MolToSmiles(rdmol),
            "mapped_smiles": offmol.to_smiles(mapped=True),
            "n_atoms": len(offmol.atoms),
            "atomic_numbers": atomic_numbers,
        }
    )
        
    return entry
            

def convert_all(
    input_directory: str = "input/XtalPi_shared_fragments",
    output_directory: str = "output/xff-td-dataset",
):
    import pyarrow as pa
    import pyarrow.dataset as ds

    output_directory = pathlib.Path(output_directory)
    output_directory.parent.mkdir(exist_ok=True, parents=True)

    input_directory = pathlib.Path(input_directory)
    combine_mol2s = list(input_directory.glob("*/*/scan*.mol2"))
    print(f"Found {len(combine_mol2s)} mol2 files")

    entries = []
    for filename in tqdm.tqdm(combine_mol2s):
        entry = parse_opt_mol2(filename)
        entries.append(entry)
    
    table = pa.Table.from_pylist(entries)
    print(table.schema)
    ds.write_dataset(table, output_directory, format="parquet")
    print(f"Wrote dataset to {output_directory}")


if __name__ == "__main__":
    convert_all()