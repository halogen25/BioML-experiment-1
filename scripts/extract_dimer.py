"""
Extract a single FabB A+B dimer plus bound ligands from a multi-chain PDB.
Saves a clean PDB with only the requested chains and any heteroatoms in those chains.
"""

import sys
from pathlib import Path

def extract_dimer(input_pdb: str, output_pdb: str, keep_chains: list[str], keep_hetatm_resnames: list[str] = None):
    keep_chains = set(keep_chains)
    keep_hetatm = set(keep_hetatm_resnames) if keep_hetatm_resnames else None

    kept = []
    with open(input_pdb) as f:
        for line in f:
            rec = line[:6].strip()

            if rec == "ATOM":
                chain = line[21]
                if chain in keep_chains:
                    kept.append(line)

            elif rec == "HETATM":
                chain = line[21]
                resname = line[17:20].strip()
                if chain in keep_chains:
                    # Always skip water; optionally filter by resname
                    if resname == "HOH":
                        continue
                    if keep_hetatm is None or resname in keep_hetatm:
                        kept.append(line)

            elif rec in ("TER", "END"):
                kept.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(kept)

    n_atom   = sum(1 for l in kept if l[:4] == "ATOM")
    n_hetatm = sum(1 for l in kept if l[:6].strip() == "HETATM")
    print(f"Wrote {output_pdb}: {n_atom} ATOM, {n_hetatm} HETATM records")


if __name__ == "__main__":
    BASE = Path(__file__).resolve().parent.parent

    # 2VBA: chains A+B (receptor for our docked poses — no ligand needed, handled separately)
    extract_dimer(
        input_pdb=str(BASE / "receptor" / "2VBA_prepared.pdb"),
        output_pdb=str(BASE / "receptor" / "2VBA_AB.pdb"),
        keep_chains=["A", "B"],
    )

    # 2AQB: chains A+B + TL6 ligand
    extract_dimer(
        input_pdb=str(BASE / "receptor" / "2AQB.pdb"),
        output_pdb=str(BASE / "receptor" / "2AQB_AB.pdb"),
        keep_chains=["A", "B"],
        keep_hetatm_resnames=["TL6"],
    )
