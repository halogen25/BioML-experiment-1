#!/bin/bash
# Full virtual screening pipeline for FabB (P0A953) inhibitor discovery
# Run from the BioML-experiment-1 directory

set -e
conda activate bioml_vscreen

echo "=== Step 1: Prepare receptor ==="
python scripts/01_prepare_receptor.py

echo "=== Step 2: Fetch natural product library ==="
python scripts/02_fetch_natural_products.py

echo "=== Step 3: Prepare ligands ==="
python scripts/03_prepare_ligands.py ligands/seed_natural_products.smi

echo "=== Step 4: Run virtual screen ==="
python scripts/04_run_vscreen.py

echo "=== Step 5: Analyze results ==="
python scripts/05_analyze_results.py

echo "=== Done! Check results/ directory for outputs ==="
