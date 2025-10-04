"""
PAINS (Pan-Assay Interference Compounds) filter patterns.

These patterns identify compounds that frequently show up as hits in high-throughput screens
but are likely to be false positives due to various interference mechanisms.
"""

PAINS_PATTERNS = {
    # Quinones and derivatives
    "quinone_A": "O=C1C=CC(=O)C=C1",
    "quinone_B": "O=c1ccc(=O)cc1",
    # Catechols and related structures
    "catechol": "c1cc(O)c(O)cc1",
    "catechol_A": "[OH]c1ccccc1[OH]",
    # Michael acceptors
    "acyl_hydrazone": "C=CC(=O)N",
    "chalcone": "O=C(C=C)c1ccccc1",
    # Rhodanines and derivatives
    "rhodanine": "O=C1CSC(=S)N1",
    "thiazolidinone": "O=C1CSC(=O)N1",
    # Alkylidene barbiturates
    "barbiturate": "O=C1CC(=O)NC(=O)N1",
    # Hydroxyphenyl hydrazones
    "phenol_hydrazone": "c1cc(O)ccc1N=N",
    # Ene_rhodanine
    "ene_rhodanine": "S=C1SC(=O)N(C)C1=C",
    # Phenolic Mannich bases
    "mannich_phenol": "c1cc(O)c(CN)cc1",
    # Curcumin-like structures
    "curcumin": "O=C(C=C)C=C(C=C)C(=O)",
    # Beta-lactams and related
    "beta_lactam": "C1C(=O)NC1",
    # Styrene-like structures
    "styrene": "C=Cc1ccccc1",
    # Nitro aromatics
    "nitro_aromatic": "c1ccc([N+](=O)[O-])cc1",
    # Azo compounds
    "azo": "cN=Nc",
    # Thiourea
    "thiourea": "NC(=S)N",
    # Isothiocyanate
    "isothiocyanate": "N=C=S",
}
