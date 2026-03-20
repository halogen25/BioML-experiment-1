"""
Step 3: Prepare ligands for docking.

Converts SMILES → 3D SDF (RDKit) → PDBQT (Meeko) for AutoDock Vina.
Skips compounds that fail to generate valid 3D conformers.
"""

import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTMolecule

BASE = Path(__file__).parent.parent
LIGAND_DIR = BASE / "ligands"
PDBQT_DIR = LIGAND_DIR / "pdbqt"
PDBQT_DIR.mkdir(exist_ok=True)


def smiles_to_pdbqt(smiles: str, name: str, max_mw: float = 600) -> str | None:
    """Convert a SMILES string to a PDBQT file. Returns path or None on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    if mw > max_mw:
        return None

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        return None

    AllChem.MMFFOptimizeMolecule(mol)

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    if not mol_setups:
        return None

    safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
    out_path = PDBQT_DIR / f"{safe_name}.pdbqt"

    from meeko import PDBQTWriterLegacy
    pdbqt_string, is_ok, err = PDBQTWriterLegacy.write_string(mol_setups[0])
    if not is_ok:
        return None

    out_path.write_text(pdbqt_string)
    return str(out_path)


def prepare_library(smi_file: str) -> list[dict]:
    """Process an SMI file and return list of successfully prepared compounds."""
    smi_path = Path(smi_file)
    lines = smi_path.read_text().strip().splitlines()

    # Handle TSV (seed library) or plain SMI (ZINC)
    has_header = lines[0].startswith("smiles") or "\t" in lines[0]
    records = []

    if has_header:
        import csv
        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            records.append({
                "smiles": row.get("smiles", row.get("SMILES", "")),
                "name": row.get("name", row.get("Name", "UNK")),
            })
    else:
        for line in lines:
            parts = line.strip().split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else f"cpd_{len(records)}"
            records.append({"smiles": smiles, "name": name})

    prepared = []
    failed = 0
    print(f"Preparing {len(records)} compounds...")

    for i, rec in enumerate(records):
        path = smiles_to_pdbqt(rec["smiles"], rec["name"])
        if path:
            prepared.append({"name": rec["name"], "smiles": rec["smiles"], "pdbqt": path})
        else:
            failed += 1
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(records)} processed, {len(prepared)} OK, {failed} failed")

    print(f"\nDone: {len(prepared)} ligands prepared, {failed} failed.")
    return prepared


if __name__ == "__main__":
    smi_file = sys.argv[1] if len(sys.argv) > 1 else str(LIGAND_DIR / "seed_natural_products.smi")
    prepared = prepare_library(smi_file)
    print(f"\nPDBQT files written to: {PDBQT_DIR}")
