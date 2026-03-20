"""Fix and dock the two compounds that failed in the initial screen."""

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

BASE = Path(__file__).resolve().parent.parent
PDBQT_DIR = BASE / "ligands" / "pdbqt"
DOCKED_DIR = BASE / "results" / "docked_poses"
RECEPTOR_PDBQT = BASE / "receptor" / "2VBA_receptor.pdbqt"
BOX_CENTER = [40.94, 69.56, 64.78]
BOX_SIZE   = [25.0, 25.0, 25.0]

compounds = [
    ("Cerulenin",      "C/C=C/C/C=C/CCC(=O)[C@@H]1[C@@H](O1)C(=O)N"),
    ("Betulinic_acid", "CC(=C)[C@@H]1CC[C@]2([C@H]1[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)C(=O)O"),
]

def smiles_to_pdbqt(smiles, name):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    if not mol_setups:
        return None
    pdbqt_string, is_ok, err = PDBQTWriterLegacy.write_string(mol_setups[0])
    if not is_ok:
        return None
    safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
    out_path = PDBQT_DIR / f"{safe_name}.pdbqt"
    out_path.write_text(pdbqt_string)
    return str(out_path)

v = Vina(sf_name="vina", cpu=0, verbosity=0)
v.set_receptor(str(RECEPTOR_PDBQT))
v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)

for name, smiles in compounds:
    print(f"\nProcessing {name}...")
    pdbqt_path = smiles_to_pdbqt(smiles, name)
    if not pdbqt_path:
        print(f"  FAILED to prepare ligand")
        continue
    print(f"  Ligand prepared: {pdbqt_path}")
    v.set_ligand_from_file(pdbqt_path)
    v.dock(exhaustiveness=8, n_poses=5)
    affinity = float(v.energies(n_poses=1)[0][0])
    print(f"  Affinity: {affinity:.2f} kcal/mol")
    out_pose = DOCKED_DIR / f"{name}_docked.pdbqt"
    v.write_poses(str(out_pose), n_poses=1, overwrite=True)
    print(f"  Pose saved: {out_pose}")
