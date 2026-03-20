"""
Step 2: Fetch a natural product ligand library.

Uses the ZINC Natural Products subset — a well-curated collection of
purchasable natural product-derived compounds in ready-to-dock SMILES format.

Starts with a manageable subset (~500 compounds) for the first screen.
For a full screen, set FULL_LIBRARY = True to download ~170k compounds.
"""

import requests
import csv
import gzip
from pathlib import Path

BASE = Path(__file__).parent.parent
LIGAND_DIR = BASE / "ligands"
LIGAND_DIR.mkdir(exist_ok=True)

# ZINC22 Natural Products subset (tranches: in-stock, natural products)
# Starting with a focused antibacterial-relevant subset from ChEMBL
FULL_LIBRARY = False

def fetch_zinc_natural_products():
    """Fetch natural product SMILES from ZINC15 API (250 compounds)."""
    output_file = LIGAND_DIR / "natural_products_zinc.smi"

    if output_file.exists():
        lines = output_file.read_text().strip().splitlines()
        print(f"Library already exists: {len(lines)} compounds in {output_file}")
        return str(output_file)

    # ZINC15 REST API — natural products subset, drug-like, in-stock
    url = (
        "https://zinc15.docking.org/substances/subsets/natural-products/"
        "?count=250&output_fields=smiles,zinc_id,name&format=csv"
    )
    print(f"Downloading ZINC15 natural products library (250 compounds)...")
    try:
        resp = requests.get(url, timeout=60, headers={"Accept": "text/csv"})
        resp.raise_for_status()
        # Validate it's actually CSV, not HTML
        if resp.text.strip().startswith("<"):
            raise ValueError("Received HTML instead of CSV — API unavailable")
        output_file.write_text(resp.text)
        lines = resp.text.strip().splitlines()
        print(f"Downloaded {len(lines) - 1} compounds.")
        return str(output_file)
    except Exception as e:
        print(f"ZINC download failed ({e}), falling back to curated seed set...")
        return fetch_seed_library()

def fetch_seed_library():
    """
    Fallback: a hand-curated set of ~30 natural products with known
    antibacterial activity (fatty acid synthesis inhibitors and related).
    Useful for validating the pipeline before a full screen.
    """
    compounds = [
        # Known FabB/FabF inhibitors and related natural products
        ("Cerulenin",           "O=C1CC(=O)/C=C/[C@@H]1CC/C=C/CCCC",              "NP001"),
        ("Platensimycin",       "O=C(O)c1cc(NC(=O)[C@@H]2CC[C@]3(C)CC[C@@H]4CC(=O)C[C@H]4[C@@H]3C2)ccc1O", "NP002"),
        ("Thiolactomycin",      "CC1=C(C(=O)O)C(C)=C(C)S1",                       "NP003"),
        ("Platencin",           "O=C(O)c1cc(NC(=O)[C@@H]2CC[C@]3(CC[C@@H]4CC(=O)C[C@H]43)C2)ccc1O", "NP004"),
        ("Epigallocatechin",    "Oc1cc(cc(O)c1O)[C@@H]1OC2=CC(=O)CC(O)C2=[C@@H]1O", "NP005"),
        ("Quercetin",           "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",    "NP006"),
        ("Berberine",           "COc1ccc2cc3c(cc2c1OC)N(C)CCC3",                  "NP007"),
        ("Carvacrol",           "Cc1ccc(C(C)C)cc1O",                               "NP008"),
        ("Thymol",              "Cc1cc(O)c(C(C)C)cc1",                             "NP009"),
        ("Eugenol",             "C=CCc1ccc(O)c(OC)c1",                             "NP010"),
        ("Resveratrol",         "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",                 "NP011"),
        ("Curcumin",            "O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1", "NP012"),
        ("Apigenin",            "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",           "NP013"),
        ("Luteolin",            "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",       "NP014"),
        ("Naringenin",          "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c21",           "NP015"),
        ("Kaempferol",          "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",       "NP016"),
        ("Catechin",            "Oc1cc2c(cc1O)C(O)Cc1cc(O)c(O)cc1-2",            "NP017"),
        ("Capsaicin",           "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",            "NP018"),
        ("Piperin",             "O=C(/C=C/C=C/c1ccc2c(c1)OCO2)N1CCCCC1",        "NP019"),
        ("Allicin",             "O=S(CC=C)SCC=C",                                  "NP020"),
        ("Cinnamaldehyde",      "O=C/C=C/c1ccccc1",                               "NP021"),
        ("Farnesol",            "CC(C)=CCC/C(C)=C/CCC(C)=CCO",                   "NP022"),
        ("Geraniol",            "CC(=CCC/C(C)=C/CO)C",                             "NP023"),
        ("Linalool",            "CC(=C)CCC(O)(C)C=C",                             "NP024"),
        ("Terpinen-4-ol",       "CC1CC(O)(C(C)C)CC=C1",                           "NP025"),
        ("Oleanolic acid",      "CC1(C)CCC2(C(=O)O)CCC3(C)C(=CC4C3CC(O)(CC4C)C(C)(C)C2)C1", "NP026"),
        ("Ursolic acid",        "CC1CCC2(C(=O)O)CCC3(C)C(=CC4C3CC(O)(CC4C)C(C)(C)C2)C1C", "NP027"),
        ("Betulinic acid",      "CC(=C)[C@@H]1CC[C@@]2(C(=O)O)CC[C@]3(C)[C@H](=C[C@@H]4[C@@]3(CC[C@@H]4[C@@H]2C1)C)C", "NP028"),
        ("Coumarin",            "O=c1ccc2ccccc2o1",                                "NP029"),
        ("Caffeic acid",        "OC(=O)/C=C/c1ccc(O)c(O)c1",                      "NP030"),
    ]

    output_file = LIGAND_DIR / "seed_natural_products.smi"
    with open(output_file, "w") as f:
        f.write("smiles\tname\tzinc_id\n")
        for name, smiles, zinc_id in compounds:
            f.write(f"{smiles}\t{name}\t{zinc_id}\n")

    print(f"Wrote {len(compounds)} seed natural product compounds to {output_file}")
    return str(output_file)

if __name__ == "__main__":
    # Try ZINC first, fall back to seed library
    library_path = fetch_zinc_natural_products()
    print(f"\nLibrary ready at: {library_path}")
