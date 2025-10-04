"""
Reactive group filter patterns.

These patterns identify functional groups that are chemically reactive and may cause
issues in biological assays or drug development.
"""

REACTIVE_PATTERNS = {
    # Electrophiles
    "acyl_halide": "[C;!$(C-[OH])][C;X3](=[OX1])[F,Cl,Br,I]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "alkyl_halide": "[CX4][F,Cl,Br,I]",
    "anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
    "epoxide": "C1OC1",
    "isocyanate": "N=C=O",
    "sulfonyl_halide": "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]",
    # Michael acceptors
    "alpha_beta_unsaturated_carbonyl": "C=CC(=O)",
    "vinyl_sulfone": "C=CS(=O)(=O)",
    "vinyl_sulfonamide": "C=CS(=O)(=O)N",
    # Peroxides
    "peroxide": "[OX2][OX2]",
    # Diazo compounds
    "diazo": "[N-]=[N+]=C",
    "diazonium": "[N+]#N",
    # Acid halides
    "sulfonyl_chloride": "S(=O)(=O)Cl",
    "phosphoryl_halide": "P(=O)[F,Cl,Br,I]",
    # Metal chelators (potential issues)
    "phosphonate": "P(=O)(O)(O)",
    "hydroxamic_acid": "C(=O)N(O)",
    # Thiols (can be reactive)
    "thiol": "[SH]",
    # Azides
    "azide": "N=[N+]=[N-]",
    # Crown ethers (binding issues)
    "crown_ether": "C1COCCOCCOCCOC1",
    # Hydrazines
    "hydrazine": "N-N",
    # Peroxy acids
    "peroxy_acid": "C(=O)OO",
}
