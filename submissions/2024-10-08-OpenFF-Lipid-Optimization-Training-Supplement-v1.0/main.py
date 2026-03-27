from cura.query import symbols_to_bits
from cura.store import Store
from cura.utils import mol_from_smiles
from openff.toolkit import Molecule
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.SimDivFilters import rdSimDivPickers
from tqdm import tqdm


def store_fragment_heavy_atoms(store, tag):
    res = store.cur.execute(
        "SELECT id, smiles, natoms FROM fragments WHERE tag = ?", (tag,)
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
    UPDATE fragments
    SET nheavy = ?
    WHERE id = ?
    """,
        data,
    )
    store.con.commit()


def store_fragment_fingerprints(store, tag):
    res = store.cur.execute(
        "SELECT id, smiles FROM fragments WHERE tag = ?", (tag,)
    )
    fpgen = AllChem.GetMorganGenerator()

    data = []
    for id, smiles in tqdm(res.fetchall()):
        mol = mol_from_smiles(smiles)
        fp = fpgen.GetFingerprint(mol)
        data.append((fp.ToBitString(), id))

    store.cur.executemany(
        """
    UPDATE fragments
    SET morgan = ?
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


def load_fragment_fingerprints(store, tag, max_nheavy, min_nheavy, elements):
    res = store.cur.execute(
        """
        SELECT id, morgan, elements
        FROM fragments
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


def load_fingerprints(store):
    res = store.cur.execute(
        "SELECT morgan FROM molecules WHERE tag = ?", (tag,)
    )
    fps = [DataStructs.CreateFromBitString(fp) for (fp,) in res.fetchall()]
    return fps


def get_smiles(store, mol_id, table="fragments"):
    res = store.cur.execute(
        f"SELECT smiles FROM {table} WHERE id = ?", (mol_id,)
    )
    return res.fetchone()[0]


def replace_dummies(mol: Chem.Mol) -> Chem.Mol:
    """Replace dummy atoms in ``mol`` with hydrogen.

    Also handle a special case where S(=O)(=O)* is replaced with a S(=O)([O-]).

    Adapted from some fragmentation code I got from Lily.
    """
    rd_dummy_replacements = [
        # Handle the special case of -S(=O)(=O)[*] -> -S(=O)(-[O-])
        (Chem.MolFromSmiles("S(=O)(=O)*"), Chem.MolFromSmiles("S(=O)([O-])")),
        # Handle the general case
        (Chem.MolFromSmiles("*"), Chem.MolFromSmiles("[H]")),
    ]

    for pat, rep in rd_dummy_replacements:
        mol = AllChem.ReplaceSubstructs(mol, pat, rep, replaceAll=True)[0]

    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    return mol


# this database is ~14 GB, so a bit difficult to include here
store = Store("/home/brent/omsf/projects/curato/store.sqlite")
tag = "data/lipidmaps.smi"

# this should contain all the inchis from the opt and td datasets I've been
# using for refitting the torsion multiplicity force field, plus the industry
# dataset
with open("inchis.dat") as inp:
    old_inchis = {inchi.strip() for inchi in inp}

lp = rdSimDivPickers.LeaderPicker()

ids, fps = load_fragment_fingerprints(
    store,
    tag,
    max_nheavy=70,
    min_nheavy=3,
    # include X to handle fragments
    elements=["Cl", "P", "Br", "I", "H", "C", "O", "N", "F", "S", "X"],
)
print(f"processing {len(fps)} fingerprints")
thresh = 0.6559
picks = lp.LazyBitVectorPick(fps, len(fps), thresh)
print(f"found {len(picks)} clusters")

fragments = (mol_from_smiles(get_smiles(store, ids[p])) for p in picks)

smiles = set()
for frag in tqdm(fragments, total=len(picks)):
    m = replace_dummies(frag)
    s = Chem.MolToSmiles(m)
    inchi = Molecule.from_smiles(s, allow_undefined_stereo=True).to_inchikey()
    if inchi not in old_inchis:
        smiles.add(s)

print(f"found {len(smiles)} deduplicated smiles")

with open("input.smi", "w") as out:
    for s in sorted(smiles, key=lambda s: len(s)):
        print(s, file=out)
