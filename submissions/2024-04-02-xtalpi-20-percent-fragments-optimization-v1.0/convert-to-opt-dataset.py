import pathlib
import re

import tqdm
import numpy as np

def parse_opt_mol2(filename: str):
    """Parse Mol2 file with custom XtalPi tags"""
    from rdkit import Chem
    from openff.toolkit import Molecule
    from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper, RDKitToolkitWrapper
    from openff.toolkit.utils.exceptions import RadicalsNotSupportedError
    
    # Use RDKit to do heavy lifting in constructing molecule
    # *However*, we have to disable sanitization because it doesn't
    # parse the UNITY_ATOM_ATTR section correctly to assign charges
    rdmol = Chem.MolFromMol2File(str(filename), sanitize=False)
    if rdmol is None:
        raise ValueError(f"Could not read {filename}")
    
    filename = pathlib.Path(filename)
    
    entry = {}
    
    # parse file into dataset
    header_pattern = "@<.+>"
    with filename.open("r") as f:
        contents = f.read()
    sections = [
        section.split("\n")
        for section in re.split(header_pattern, contents)
        if section
    ]
    
    try:
        _to_array = lambda x: np.array(list(map(float, x)))
    except Exception:
        raise Exception(filename)

    # parse sections manually
    IGNORE = ["MOLECULE", "BOND"]
    CHARGES = ["AM1BCC_CHARGE", "RESP_CHARGE"]
    COORDINATES = ["ATOM", "INPUT_MOLECULE"]
    ARRAY_TYPES = ["HESSIAN", "GRADIENT"]
    FLOAT_TYPES = ["CORE_TIME", "ENERGY"]
    STRING_TYPES = ["CA_TYPE", "FINGERPRINT"]
    seen_unity_atom_attr = False
    for section in sections:
        header = section[0]
        if header in IGNORE:
            continue
        lines = [x for x in section[1:] if x]
        if header in CHARGES:
            charge_col = [line.strip().split()[1] for line in lines]
            charges_ = _to_array(charge_col)
            entry[header] = charges_
        # we have to manually assign formal charges
        elif header == "UNITY_ATOM_ATTR":
            seen_unity_atom_attr = True
            while lines:
                atom_index = int(lines[0].strip().split()[0]) - 1
                field, charge_ = lines[1].strip().split()
                if field == "charge":
                    rdatom = rdmol.GetAtomWithIdx(atom_index)
                    rdatom.SetFormalCharge(int(float(charge_)))
                lines = lines[2:]
        elif header in COORDINATES:
            coord_col = [
                x for line in lines
                for x in line.strip().split()[2:5]
            ]
            coords_ = _to_array(coord_col)
            entry[header] = coords_
        elif header in ARRAY_TYPES:
            entry[header] = _to_array(" ".join(lines).split())
        elif header in FLOAT_TYPES:
            assert len(lines) == 1
            entry[header] = float(lines[0].strip())
        elif header in STRING_TYPES:
            assert len(lines) == 1
            entry[header] = lines[0].strip()
    
    Chem.SanitizeMol(rdmol)
    toolkit_registry = RDKitToolkitWrapper()
    try:
        offmol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    except RadicalsNotSupportedError:
        if not seen_unity_atom_attr:
            # can parse these incorrectly at times, it seems?
            offmol = Molecule.from_file(str(filename), "MOL2")
            toolkit_registry = OpenEyeToolkitWrapper()
        else:
            raise

    atomic_numbers = sorted(set([a.atomic_number for a in offmol.atoms]))

    entry.update(
        {
            "filename": filename.stem,
            "directory": filename.parent.name,
            "parent": filename.parent.parent.name,
            "smiles": Chem.MolToSmiles(rdmol),
            "mapped_smiles": offmol.to_smiles(mapped=True, toolkit_registry=toolkit_registry),
            "n_atoms": len(offmol.atoms),
            "atomic_numbers": atomic_numbers,
        }
    )
        
    return entry
            

def convert_all(
    input_directory: str = "input/XtalPi_20_percent_training_set", # "input/XtalPi_shared_fragments",
    output_directory: str = "output/xff-20-percent-opt-dataset",
):
    import pyarrow as pa
    import pyarrow.dataset as ds

    output_directory = pathlib.Path(output_directory)
    output_directory.parent.mkdir(exist_ok=True, parents=True)

    input_directory = pathlib.Path(input_directory)
    combine_mol2s = list(input_directory.glob("*/*/combine*.mol2"))[::-1]
    print(f"Found {len(combine_mol2s)} mol2 files")

    entries = []
    for filename in tqdm.tqdm(combine_mol2s):
        try:
            entry = parse_opt_mol2(filename)
        except Exception as e:
            print(f"Error parsing {filename}: {e}")
            raise e
        entries.append(entry)
    
    table = pa.Table.from_pylist(entries)
    print(table.schema)
    ds.write_dataset(table, output_directory, format="parquet")
    print(f"Wrote dataset to {output_directory}")


if __name__ == "__main__":
    convert_all()
