# BioML Experiment 1 — Virtual Screen for FabB Inhibitors

Virtual screening for inhibitors of **FabB** (3-oxoacyl-ACP synthase I, UniProt [P0A953](https://www.uniprot.org/uniprot/P0A953)) — an essential enzyme in *E. coli* fatty acid biosynthesis and validated antibiotic target.

## Target
| Field | Value |
|---|---|
| Protein | FabB (3-oxoacyl-ACP synthase I) |
| UniProt | P0A953 |
| Organism | *Escherichia coli* K-12 |
| Function | Catalyzes the elongation step in bacterial fatty acid synthesis |
| Structure | PDB [2VBA](https://www.rcsb.org/structure/2VBA) — 1.36 Å resolution |
| Drug target class | Antibacterial — fatty acid synthesis |

## Pipeline

```
01_prepare_receptor.py   — Fix/clean 2VBA crystal structure (PDBFixer)
02_fetch_natural_products.py — Download natural product library (ZINC / seed set)
03_prepare_ligands.py    — SMILES → 3D conformers → PDBQT (RDKit + Meeko)
04_run_vscreen.py        — AutoDock Vina docking (binding site from P4T co-crystal)
05_analyze_results.py    — Rank hits, filter by Lipinski, draw top structures
```

## Docking Box
Centered on the co-crystallized inhibitor P4T in PDB 2VBA:
- Center: x=40.94, y=69.56, z=64.78
- Size: 25 × 25 × 25 Å

## Environment Setup
```bash
conda create -n bioml_vscreen -c conda-forge -c bioconda vina rdkit meeko pdbfixer openbabel python=3.10
conda activate bioml_vscreen
```

## Run
```bash
bash scripts/run_pipeline.sh
```

## Tools Used
- [AutoDock Vina 1.2.6](https://github.com/ccsb-scripps/AutoDock-Vina) — molecular docking
- [RDKit](https://www.rdkit.org/) — ligand 3D generation
- [Meeko](https://github.com/forlilab/Meeko) — ligand PDBQT preparation
- [PDBFixer](https://github.com/openmm/pdbfixer) — receptor preparation
- [ZINC Natural Products](https://zinc22.docking.org/) — ligand library
- [LabDAO PLEX](https://docs.labdao.xyz/) — P2P compute orchestration
