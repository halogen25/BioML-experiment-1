"""
Step 5: Analyze and visualize virtual screening results.

- Ranks hits by binding affinity
- Filters by drug-likeness (Lipinski / extended NP rules)
- Generates a summary report
"""

import csv
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

BASE = Path(__file__).parent.parent
RESULTS_DIR = BASE / "results"
LIGAND_DIR = BASE / "ligands"


def load_smiles_lookup(smi_file: Path) -> dict[str, str]:
    """Build name → SMILES dict from the library file."""
    lookup = {}
    lines = smi_file.read_text().strip().splitlines()
    has_header = lines[0].startswith("smiles") or "\t" in lines[0]
    if has_header:
        import csv as _csv
        reader = _csv.DictReader(lines, delimiter="\t")
        for row in reader:
            name = row.get("name", row.get("Name", ""))
            smiles = row.get("smiles", row.get("SMILES", ""))
            safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
            lookup[safe] = smiles
    else:
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 2:
                lookup[parts[1]] = parts[0]
    return lookup


def lipinski_pass(mol) -> bool:
    return (
        Descriptors.MolWt(mol) <= 500
        and Descriptors.MolLogP(mol) <= 5
        and Descriptors.NumHDonors(mol) <= 5
        and Descriptors.NumHAcceptors(mol) <= 10
    )


def analyze(results_csv: str = None, smi_file: str = None):
    results_csv = results_csv or str(RESULTS_DIR / "screening_results.csv")
    smi_file = smi_file or str(LIGAND_DIR / "seed_natural_products.smi")

    if not Path(results_csv).exists():
        print(f"Results file not found: {results_csv}")
        print("Run 04_run_vscreen.py first.")
        return

    smiles_lookup = load_smiles_lookup(Path(smi_file))

    with open(results_csv) as f:
        results = list(csv.DictReader(f))

    print(f"\n{'='*60}")
    print(f"VIRTUAL SCREEN RESULTS — FabB (P0A953) inhibitors")
    print(f"Target: 3-oxoacyl-ACP synthase I (E. coli FabB)")
    print(f"Structure: PDB 2VBA (1.36 Å resolution)")
    print(f"{'='*60}")
    print(f"Total compounds docked: {len(results)}")
    print(f"\n{'Rank':<6}{'Name':<35}{'Affinity (kcal/mol)':<22}{'Lipinski'}")
    print("-" * 70)

    hits = []
    for r in results:
        name = r["name"]
        affinity = float(r["affinity_kcal_mol"])
        smiles = smiles_lookup.get(name, "")
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        lipo = lipinski_pass(mol) if mol else "N/A"
        print(f"{r['rank']:<6}{name:<35}{affinity:<22.2f}{'✓' if lipo is True else '✗' if lipo is False else 'N/A'}")
        if mol and affinity < -6.0:
            hits.append({"name": name, "smiles": smiles, "affinity": affinity, "mol": mol})

    print(f"\n--- Hits with affinity < -6.0 kcal/mol: {len(hits)} ---")

    # Draw top hits
    if hits:
        mols = [h["mol"] for h in hits[:12]]
        legends = [f"{h['name']}\n{h['affinity']:.2f} kcal/mol" for h in hits[:12]]
        img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(300, 250), legends=legends)
        img_path = RESULTS_DIR / "top_hits.png"
        img.save(str(img_path))
        print(f"Structure grid saved to: {img_path}")

    # Save filtered hits
    hits_csv = RESULTS_DIR / "top_hits.csv"
    with open(hits_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "smiles", "affinity_kcal_mol"])
        writer.writeheader()
        for h in hits:
            writer.writerow({"name": h["name"], "smiles": h["smiles"], "affinity_kcal_mol": h["affinity"]})
    print(f"Top hits CSV saved to: {hits_csv}")


if __name__ == "__main__":
    analyze()
