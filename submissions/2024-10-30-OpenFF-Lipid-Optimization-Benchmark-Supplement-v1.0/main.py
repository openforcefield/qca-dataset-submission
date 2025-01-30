from cura.query import symbols_to_bits
from cura.store import Store
from cura.utils import mol_from_smiles
from openff.toolkit import Molecule
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.SimDivFilters import rdSimDivPickers
from tqdm import tqdm


def store_heavy_atoms(store, tag):
    res = store.cur.execute(
        "SELECT id, smiles, natoms FROM molecules WHERE tag = ?", (tag,)
    )
    data = []
    for id, smiles, natoms in tqdm(res.fetchall()):
        mol = mol_from_smiles(smiles)
        natoms_ = mol.GetNumAtoms()
        nheavy = mol.GetNumHeavyAtoms()
        assert natoms == natoms_

        data.append((nheavy, id))

    store.cur.executemany(
        """
    UPDATE molecules
    SET nheavy = ?
    WHERE id = ?
    """,
        data,
    )
    store.con.commit()


def store_fingerprints(store):
    res = store.cur.execute(
        "SELECT id, smiles FROM molecules WHERE tag = ?", (tag,)
    )
    fpgen = AllChem.GetMorganGenerator()

    for id, smiles in tqdm(res.fetchall()):
        mol = mol_from_smiles(smiles)
        fp = fpgen.GetFingerprint(mol)
        store.cur.execute(
            """
        UPDATE molecules
        SET morgan = ?
        WHERE id = ?
        """,
            (fp.ToBitString(), id),
        )
        store.con.commit()


def load_fingerprints(store, tag, max_nheavy, min_nheavy, elements):
    res = store.cur.execute(
        """
        SELECT id, morgan, elements
        FROM molecules
        WHERE tag = ?
        AND nheavy < ?
        AND nheavy > ?
        """,
        (tag, max_nheavy, min_nheavy),
    )

    # element filtering adapted from cura ElementFilter and Store.get_molecules
    # (for int conversion)
    mask = symbols_to_bits(elements)

    ret = [
        (id, DataStructs.CreateFromBitString(fp))
        for id, fp, elements in res.fetchall()
        if (int.from_bytes(elements, "big") | mask) == mask
    ]

    ids, fps = zip(*ret)
    return ids, fps


def get_smiles(store, mol_id, table):
    res = store.cur.execute(
        f"SELECT smiles FROM {table} WHERE id = ?", (mol_id,)
    )
    return res.fetchone()[0]


# this database is ~14 GB, so a bit difficult to include here
store = Store("/home/brent/omsf/projects/curato/store.sqlite")
tag = "data/lipidmaps.smi"

# this should contain all the inchis from the opt and td datasets I've been
# using for refitting the torsion multiplicity force field, plus the industry
# dataset
with open("inchis.dat") as inp:
    old_inchis = {inchi.strip() for inchi in inp}

# update previous inchis with those from the fragment dataset
with open(
    "../2024-10-08-OpenFF-Lipid-Optimization-Training-Supplement-v1.0/output.smi"
) as inp:
    for s in inp:
        mol = Molecule.from_smiles(s.strip(), allow_undefined_stereo=True)
        old_inchis.add(mol.to_inchikey())

lp = rdSimDivPickers.LeaderPicker()

ids, fps = load_fingerprints(
    store,
    tag,
    max_nheavy=100,
    min_nheavy=3,
    elements=["Cl", "P", "Br", "I", "H", "C", "O", "N", "F", "S"],
)
print(f"processing {len(fps)} fingerprints")
thresh = 0.708
picks = lp.LazyBitVectorPick(fps, len(fps), thresh)
print(f"found {len(picks)} clusters")

molecules = (
    mol_from_smiles(get_smiles(store, ids[p], "molecules")) for p in picks
)

smiles = set()
for m in tqdm(molecules, total=len(picks)):
    s = Chem.MolToSmiles(m)
    inchi = Molecule.from_smiles(s, allow_undefined_stereo=True).to_inchikey()
    if inchi not in old_inchis:
        smiles.add(s)

print(f"found {len(smiles)} deduplicated smiles")

with open("input.smi", "w") as out:
    for s in sorted(smiles, key=lambda s: len(s)):
        print(s, file=out)
