"""
Heterocycle filter patterns.

These patterns identify common heterocyclic scaffolds in medicinal chemistry.
"""

HETEROCYCLE_PATTERNS = {
    # Five-membered rings with one heteroatom
    "pyrrole": "c1cc[nH]c1",
    "furan": "c1ccoc1",
    "thiophene": "c1ccsc1",
    
    # Five-membered rings with two heteroatoms
    "imidazole": "c1cnc[nH]1",
    "pyrazole": "c1cn[nH]c1",
    "oxazole": "c1cnco1",
    "isoxazole": "c1cnoc1",
    "thiazole": "c1cncs1",
    
    # Six-membered rings with one heteroatom
    "pyridine": "c1ccncc1",
    "pyran": "C1CCOCC1",
    
    # Six-membered rings with two heteroatoms
    "pyrimidine": "c1cncnc1",
    "pyrazine": "c1cnccn1",
    "pyridazine": "c1cnncc1",
    "oxazine": "C1COCNC1",
    
    # Six-membered rings with three heteroatoms
    "triazine": "c1ncncn1",
    
    # Fused bicyclic systems
    "indole": "c1ccc2[nH]ccc2c1",
    "benzofuran": "c1ccc2occc2c1",
    "benzothiophene": "c1ccc2sccc2c1",
    "benzimidazole": "c1ccc2[nH]cnc2c1",
    "quinoline": "c1ccc2ncccc2c1",
    "isoquinoline": "c1ccc2cnccc2c1",
    "quinazoline": "c1ccc2ncncc2c1",
    "quinoxaline": "c1ccc2nccnc2c1",
    "purine": "c1nc2[nH]cnc2n1",
    
    # Other important heterocycles
    "piperidine": "C1CCNCC1",
    "piperazine": "C1CNCCN1",
    "morpholine": "C1COCCN1",
    "thiomorpholine": "C1CSCCN1",
    "pyrrolidine": "C1CCNC1",
    
    # Seven-membered rings
    "azepine": "C1=CC=CNCC1",
    "oxepine": "C1=CC=COCC1",
}
