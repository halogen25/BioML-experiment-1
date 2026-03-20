"""
Extract a single FabB A+B dimer plus bound ligands from a multi-chain PDB.
Saves a clean PDB with only the requested chains, ligands, and CONECT records.
"""

from pathlib import Path


def extract_dimer(input_pdb: str, output_pdb: str, keep_chains: list, keep_hetatm_resnames: list = None):
    keep_chains = set(keep_chains)
    keep_hetatm = set(keep_hetatm_resnames) if keep_hetatm_resnames else None

    kept = []
    kept_serials = set()

    with open(input_pdb) as f:
        for line in f:
            rec = line[:6].strip()

            if rec == "ATOM":
                chain = line[21]
                if chain in keep_chains:
                    kept.append(line)
                    kept_serials.add(line[6:11].strip())

            elif rec == "HETATM":
                chain = line[21]
                resname = line[17:20].strip()
                if chain in keep_chains and resname != "HOH":
                    if keep_hetatm is None or resname in keep_hetatm:
                        kept.append(line)
                        kept_serials.add(line[6:11].strip())

            elif rec == "CONECT":
                # Keep CONECT if all referenced atoms are in our selection
                serials = [line[6:11].strip(), line[11:16].strip(),
                           line[16:21].strip(), line[21:26].strip(),
                           line[26:31].strip()]
                serials = [s for s in serials if s]
                if serials[0] in kept_serials:
                    kept.append(line)

            elif rec in ("TER", "END"):
                kept.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(kept)

    n_atom   = sum(1 for l in kept if l[:4] == "ATOM")
    n_hetatm = sum(1 for l in kept if l[:6].strip() == "HETATM")
    n_conect = sum(1 for l in kept if l[:6].strip() == "CONECT")
    print(f"Wrote {output_pdb}: {n_atom} ATOM, {n_hetatm} HETATM, {n_conect} CONECT records")


if __name__ == "__main__":
    BASE = Path(__file__).resolve().parent.parent

    extract_dimer(
        input_pdb=str(BASE / "receptor" / "2VBA_prepared.pdb"),
        output_pdb=str(BASE / "receptor" / "2VBA_AB.pdb"),
        keep_chains=["A", "B"],
    )

    extract_dimer(
        input_pdb=str(BASE / "receptor" / "2AQB.pdb"),
        output_pdb=str(BASE / "receptor" / "2AQB_AB.pdb"),
        keep_chains=["A", "B"],
        keep_hetatm_resnames=["TL6"],
    )

    extract_dimer(
        input_pdb=str(BASE / "receptor" / "1FJ8.pdb"),
        output_pdb=str(BASE / "receptor" / "1FJ8_AB.pdb"),
        keep_chains=["A", "B"],
        keep_hetatm_resnames=["CER"],
    )
