"""
Step 1: Prepare the FabB (P0A953) receptor for docking.

- Removes heteroatoms (ligands, water) from the crystal structure
- Fixes missing residues/atoms with PDBFixer
- Saves a clean receptor PDB ready for AutoDock Vina
"""

import os
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

BASE = Path(__file__).parent.parent
RECEPTOR_DIR = BASE / "receptor"

def prepare_receptor(input_pdb: str, output_pdb: str):
    print(f"Loading {input_pdb}...")
    fixer = PDBFixer(filename=input_pdb)

    print("Removing heterogens (keeping water: False)...")
    fixer.removeHeterogens(keepWater=False)

    print("Finding missing residues...")
    fixer.findMissingResidues()

    print("Finding missing atoms...")
    fixer.findMissingAtoms()

    print("Adding missing atoms...")
    fixer.addMissingAtoms()

    print("Adding missing hydrogens (pH 7.4)...")
    fixer.addMissingHydrogens(7.4)

    print(f"Saving cleaned receptor to {output_pdb}...")
    with open(output_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print("Done.")

if __name__ == "__main__":
    prepare_receptor(
        input_pdb=str(RECEPTOR_DIR / "2VBA.pdb"),
        output_pdb=str(RECEPTOR_DIR / "2VBA_prepared.pdb"),
    )
