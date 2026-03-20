"""
Step 6: Compile and export the complete screening results table.

Merges all 30 compounds (including the 2 fixed ones) into a single ranked CSV
with SMILES, compound class, Lipinski pass/fail, and affinity.
"""

import csv
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors

BASE = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE / "results"

# Full compound data: name, SMILES, class, corrected affinity
COMPOUNDS = [
    ("Platensimycin",   "O=C(O)c1cc(NC(=O)[C@@H]2CC[C@]3(C)CC[C@@H]4CC(=O)C[C@H]4[C@@H]3C2)ccc1O",  "Terpenoid / FabB inhibitor",       -8.963),
    ("Quercetin",       "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",                                  "Flavonoid",                        -8.837),
    ("Luteolin",        "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",                                     "Flavonoid",                        -8.657),
    ("Catechin",        "Oc1cc2c(cc1O)C(O)Cc1cc(O)c(O)cc1-2",                                          "Flavonoid",                        -8.411),
    ("Resveratrol",     "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",                                              "Stilbenoid",                       -8.395),
    ("Kaempferol",      "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",                                     "Flavonoid",                        -8.365),
    ("Apigenin",        "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",                                        "Flavonoid",                        -8.339),
    ("Piperin",         "O=C(/C=C/C=C/c1ccc2c(c1)OCO2)N1CCCCC1",                                      "Alkaloid",                         -8.303),
    ("Platencin",       "O=C(O)c1cc(NC(=O)[C@@H]2CC[C@]3(CC[C@@H]4CC(=O)C[C@H]43)C2)ccc1O",          "Terpenoid / FabB inhibitor",       -8.302),
    ("Capsaicin",       "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",                                         "Vanilloid alkaloid",               -8.249),
    ("Naringenin",      "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c21",                                        "Flavanone",                        -8.151),
    ("Curcumin",        "O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1",                       "Diarylheptanoid",                  -8.024),
    ("Berberine",       "COc1ccc2cc3c(cc2c1OC)N(C)CCC3",                                              "Isoquinoline alkaloid",            -7.509),
    ("Betulinic_acid",  "CC(=C)[C@@H]1CC[C@]2([C@H]1[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)C(=O)O", "Lupane triterpenoid", -7.940),
    ("Cerulenin",       "C/C=C/C/C=C/CCC(=O)[C@@H]1[C@@H](O1)C(=O)N",                                "Epoxy amide / FabB inhibitor",     -7.010),
    ("Farnesol",        "CC(=CCC/C(C)=C/CCC(C)=CCO)C",                                                "Sesquiterpene",                    -7.245),
    ("Caffeic_acid",    "OC(=O)/C=C/c1ccc(O)c(O)c1",                                                  "Hydroxycinnamic acid",             -7.086),
    ("Coumarin",        "O=c1ccc2ccccc2o1",                                                            "Benzopyrone",                      -6.727),
    ("Geraniol",        "CC(=CCC/C(C)=C/CO)C",                                                        "Monoterpene",                      -6.206),
    ("Carvacrol",       "Cc1ccc(C(C)C)cc1O",                                                          "Monoterpenoid phenol",             -6.188),
    ("Thymol",          "Cc1cc(O)c(C(C)C)cc1",                                                        "Monoterpenoid phenol",             -6.070),
    ("Cinnamaldehyde",  "O=C/C=C/c1ccccc1",                                                           "Phenylpropanoid",                  -6.053),
    ("Terpinen-4-ol",   "CC1CC(O)(C(C)C)CC=C1",                                                       "Monoterpenoid",                    -5.986),
    ("Eugenol",         "C=CCc1ccc(O)c(OC)c1",                                                        "Phenylpropanoid",                  -5.918),
    ("Thiolactomycin",  "CC1=C(C(=O)O)C(C)=C(C)S1",                                                  "Thiolactone / FabB inhibitor",     -5.523),
    ("Linalool",        "CC(=C)CCC(O)(C)C=C",                                                         "Monoterpene alcohol",              -5.160),
    ("Oleanolic_acid",  "CC1(C)CCC2(C(=O)O)CCC3(C)C(=CC4C3CC(O)(CC4C)C(C)(C)C2)C1",                 "Oleanane triterpenoid",            -4.286),
    ("Allicin",         "O=S(CC=C)SCC=C",                                                             "Organosulfur",                     -4.052),
    ("Ursolic_acid",    "CC1CCC2(C(=O)O)CCC3(C)C(=CC4C3CC(O)(CC4C)C(C)(C)C2)C1C",                   "Ursane triterpenoid",              -3.607),
    ("Thiolactomycin",  "CC1=C(C(=O)O)C(C)=C(C)S1",                                                  "Thiolactone / FabB inhibitor",     -5.523),
]

# Deduplicate (Thiolactomycin listed twice above — remove)
seen = set()
unique = []
for row in COMPOUNDS:
    if row[0] not in seen:
        seen.add(row[0])
        unique.append(row)
COMPOUNDS = unique

def lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "N/A"
    passes = (
        Descriptors.MolWt(mol) <= 500
        and Descriptors.MolLogP(mol) <= 5
        and Descriptors.NumHDonors(mol) <= 5
        and Descriptors.NumHAcceptors(mol) <= 10
    )
    return "Pass" if passes else "Fail"

def mw(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return f"{Descriptors.MolWt(mol):.1f}" if mol else "N/A"

def logp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return f"{Descriptors.MolLogP(mol):.2f}" if mol else "N/A"

# Sort by affinity
COMPOUNDS.sort(key=lambda x: x[3])

rows = []
for rank, (name, smiles, cls, aff) in enumerate(COMPOUNDS, 1):
    rows.append({
        "rank":               rank,
        "name":               name,
        "class":              cls,
        "affinity_kcal_mol":  f"{aff:.3f}",
        "mol_weight":         mw(smiles),
        "logP":               logp(smiles),
        "lipinski":           lipinski(smiles),
        "smiles":             smiles,
    })

out_csv = RESULTS_DIR / "all_compounds_ranked.csv"
fields = ["rank", "name", "class", "affinity_kcal_mol", "mol_weight", "logP", "lipinski", "smiles"]

with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

# Print table
col_w = [5, 22, 30, 20, 12, 8, 10]
header = f"{'Rank':<{col_w[0]}} {'Name':<{col_w[1]}} {'Class':<{col_w[2]}} {'Affinity (kcal/mol)':<{col_w[3]}} {'MW (Da)':<{col_w[4]}} {'logP':<{col_w[5]}} {'Lipinski':<{col_w[6]}}"
print(header)
print("-" * len(header))
for r in rows:
    print(f"{r['rank']:<{col_w[0]}} {r['name']:<{col_w[1]}} {r['class']:<{col_w[2]}} {r['affinity_kcal_mol']:<{col_w[3]}} {r['mol_weight']:<{col_w[4]}} {r['logP']:<{col_w[5]}} {r['lipinski']:<{col_w[6]}}")

print(f"\nExported to: {out_csv}")
