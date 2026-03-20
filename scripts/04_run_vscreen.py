"""
Step 4: Run virtual screen with AutoDock Vina.

Docks all prepared PDBQT ligands against FabB (P0A953, PDB: 2VBA).

Binding box is centered on the co-crystallized P4T ligand active site:
  Center: x=40.94, y=69.56, z=64.78
  Size:   25 x 25 x 25 Å (covers full active site + tunnel)

Results are ranked by best binding affinity (kcal/mol, lower = better).
"""

import subprocess
import json
import csv
from pathlib import Path
from vina import Vina

BASE = Path(__file__).parent.parent
RECEPTOR_DIR = BASE / "receptor"
LIGAND_DIR = BASE / "ligands" / "pdbqt"
RESULTS_DIR = BASE / "results"
RESULTS_DIR.mkdir(exist_ok=True)
DOCKED_DIR = RESULTS_DIR / "docked_poses"
DOCKED_DIR.mkdir(exist_ok=True)

# Binding box: centered on P4T co-crystal ligand in 2VBA
BOX_CENTER = [40.94, 69.56, 64.78]
BOX_SIZE   = [25.0, 25.0, 25.0]

RECEPTOR_PDBQT = RECEPTOR_DIR / "2VBA_receptor.pdbqt"


def convert_receptor_to_pdbqt(receptor_pdb: Path, out_pdbqt: Path):
    """Use obabel to convert prepared PDB → PDBQT (remove H, add Gasteiger charges)."""
    if out_pdbqt.exists():
        print(f"Receptor PDBQT already exists: {out_pdbqt}")
        return
    print("Converting receptor PDB → PDBQT...")
    cmd = [
        "obabel", str(receptor_pdb),
        "-O", str(out_pdbqt),
        "-xr",        # rigid receptor
        "--partialcharge", "gasteiger",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"obabel failed:\n{result.stderr}")
    print(f"Receptor PDBQT written to {out_pdbqt}")


def dock_ligand(v: Vina, ligand_pdbqt: Path) -> float | None:
    """Dock a single ligand. Returns best affinity (kcal/mol) or None on failure."""
    try:
        v.set_ligand_from_file(str(ligand_pdbqt))
        v.dock(exhaustiveness=8, n_poses=5)
        energies = v.energies(n_poses=1)
        return float(energies[0][0])  # best pose affinity
    except Exception as e:
        print(f"  Docking failed for {ligand_pdbqt.stem}: {e}")
        return None


def run_screen():
    receptor_pdb = RECEPTOR_DIR / "2VBA_prepared.pdb"
    if not receptor_pdb.exists():
        raise FileNotFoundError(
            f"Prepared receptor not found: {receptor_pdb}\n"
            "Run 01_prepare_receptor.py first."
        )

    convert_receptor_to_pdbqt(receptor_pdb, RECEPTOR_PDBQT)

    ligand_files = sorted(LIGAND_DIR.glob("*.pdbqt"))
    if not ligand_files:
        raise FileNotFoundError(
            f"No ligand PDBQT files found in {LIGAND_DIR}\n"
            "Run 03_prepare_ligands.py first."
        )

    print(f"\nScreening {len(ligand_files)} ligands against FabB (2VBA)...")
    print(f"Box center: {BOX_CENTER}, size: {BOX_SIZE} Å\n")

    v = Vina(sf_name="vina", cpu=0, verbosity=0)  # cpu=0 → use all cores
    v.set_receptor(str(RECEPTOR_PDBQT))
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)

    results = []
    for i, lig in enumerate(ligand_files):
        affinity = dock_ligand(v, lig)
        if affinity is not None:
            results.append({"rank": None, "name": lig.stem, "affinity_kcal_mol": affinity})
            print(f"  [{i+1}/{len(ligand_files)}] {lig.stem:40s} {affinity:.2f} kcal/mol")

            # Save best pose
            out_pose = DOCKED_DIR / f"{lig.stem}_docked.pdbqt"
            v.write_poses(str(out_pose), n_poses=1, overwrite=True)

    # Rank by affinity (most negative = best)
    results.sort(key=lambda r: r["affinity_kcal_mol"])
    for i, r in enumerate(results):
        r["rank"] = i + 1

    # Save CSV summary
    csv_out = RESULTS_DIR / "screening_results.csv"
    with open(csv_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "name", "affinity_kcal_mol"])
        writer.writeheader()
        writer.writerows(results)

    print(f"\n--- Top 10 Hits ---")
    for r in results[:10]:
        print(f"  #{r['rank']:2d}  {r['name']:40s}  {r['affinity_kcal_mol']:.2f} kcal/mol")

    print(f"\nFull results saved to: {csv_out}")
    print(f"Docked poses saved to: {DOCKED_DIR}")
    return results


if __name__ == "__main__":
    run_screen()
