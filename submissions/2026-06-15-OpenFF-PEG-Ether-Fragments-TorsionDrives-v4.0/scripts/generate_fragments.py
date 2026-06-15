"""Generate the PEG ether fragment library and write ``input.json``.

This is the fragment-generation code used to build this dataset's molecule set
for a polymer fitting demo. It is self-contained (only RDKit is needed).
Running it reproduces ``../input.json``, the input the notebook loads.

    python scripts/generate_fragments.py
"""

import json
import pathlib

from rdkit import Chem

# --- targets ----------------------------------------------------------------
# The first N_TORSIONS molecules of the library are the ones torsion-driven.
N_TORSIONS = 50

# Molecules are *fragments* of a long PEG chain -- T_left-(O-CH2-CH2)_n-O-T_right
# -- sampled with variable length, end-group registration (hydroxyl OR sp3
# carbon, independently on each end) and the occasional methyl branch.
TERMINI = ["", "C", "CC", "CCC", "C(C)C"]  # "" = hydroxyl end; else methyl/ethyl/...
N_UNITS = [1, 2, 3, 4]  # ethylene-glycol repeats in the fragment

# Driven dihedral preference: a CX4-flanked O-CH2-CH2-O if one exists,
# else any glycol O-CH2-CH2-O (covers hydroxyl-terminated fragments).
_CORE_STRICT = Chem.MolFromSmarts("[CX4][OX2][CX4H2][CX4H2][OX2][CX4]")
_CORE_ANY = Chem.MolFromSmarts("[OX2][CX4H2][CX4H2][OX2]")


def _peg_smiles(left: str, n: int, right: str, branched: bool = False) -> str:
    """A PEG-chain fragment T_left-(O-CH2-CH2)_n-O-T_right. Empty terminus
    leaves a hydroxyl. `branched` puts a methyl on the first unit's second
    carbon (a substituent away from the driven dihedral)."""
    first = "OCC(C)" if branched else "OCC"
    return left + first + "OCC" * (n - 1) + "O" + right


def central_dihedral(rdmol: Chem.Mol):
    """Return the (O, C, C, O) indices of a glycol dihedral to drive,
    preferring one flanked by CX4 on both sides, else any O-CH2-CH2-O."""
    match = rdmol.GetSubstructMatch(_CORE_STRICT)
    if match:
        # match = (C, O, CH2, CH2, O, C); the driven dihedral is O-C-C-O.
        return match[1], match[2], match[3], match[4]
    match = rdmol.GetSubstructMatch(_CORE_ANY)  # (O, CH2, CH2, O)
    return tuple(match) if match else None


def build_library() -> list[str]:
    """Canonical SMILES for PEG-chain fragments across lengths/termini, plus
    singly methyl-branched variants. Kept iff a glycol O-CH2-CH2-O torsion
    exists to drive (see central_dihedral)."""
    smis, seen = [], set()
    for n in N_UNITS:
        for left in TERMINI:
            for right in TERMINI:
                cands = [_peg_smiles(left, n, right)]
                if n >= 2:  # keep an unbranched dihedral available elsewhere
                    cands.append(_peg_smiles(left, n, right, branched=True))
                for cand in cands:
                    rdmol = Chem.MolFromSmiles(cand)
                    if rdmol is None or central_dihedral(rdmol) is None:
                        continue
                    canon = Chem.MolToSmiles(rdmol)
                    if canon not in seen:
                        seen.add(canon)
                        smis.append(canon)
    return smis


def _prepare(smiles: str):
    """Return (rdkit_mol_with_explicit_Hs, mapped_smiles).

    Atoms are tagged with map numbers = idx+1 so the mapped SMILES encodes
    this exact atom order; downstream `Molecule.from_mapped_smiles` then
    recovers the same order as the stored torsion indices.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    return mol, Chem.MolToSmiles(mol)


def main() -> None:
    here = pathlib.Path(__file__).resolve().parents[1]
    json_out = here / "input.json"
    smi_out = here / "input.smi"

    records = []
    smiles_list = []
    for smiles in build_library()[:N_TORSIONS]:
        rdmol, mapped = _prepare(smiles)
        idxs = central_dihedral(rdmol)
        if idxs is None:
            continue
        records.append(
            {
                "mapped_smiles": mapped,
                "dihedral_indices": [int(i) for i in idxs],
            }
        )
        # build_library returns canonical SMILES; keep input.smi (a readable,
        # unmapped companion) in lockstep with input.json.
        smiles_list.append(smiles)

    with open(json_out, "w") as fh:
        json.dump(records, fh, indent=2)
    smi_out.write_text("\n".join(smiles_list) + "\n")

    print(f"Wrote {len(records)} records to {json_out}")
    print(f"Wrote {len(smiles_list)} SMILES to {smi_out}")


if __name__ == "__main__":
    main()
