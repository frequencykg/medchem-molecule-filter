Filter Rules Reference
======================

This page documents all the filter patterns used in the library.

PAINS Patterns
--------------

Pan-Assay Interference Compounds (PAINS) are molecules that frequently show up as hits in high-throughput screens but are likely false positives. These patterns identify problematic structures.

Quinones and Derivatives
~~~~~~~~~~~~~~~~~~~~~~~~

* **quinone_A**: ``O=C1C=CC(=O)C=C1`` - p-Benzoquinone pattern
* **quinone_B**: ``O=c1ccc(=O)cc1`` - Alternative quinone representation

Catechols
~~~~~~~~~

* **catechol**: ``c1cc(O)c(O)cc1`` - Catechol (1,2-dihydroxybenzene)
* **catechol_A**: ``[OH]c1ccccc1[OH]`` - Alternative catechol pattern

Michael Acceptors
~~~~~~~~~~~~~~~~~

* **acyl_hydrazone**: ``C=CC(=O)N`` - Acyl hydrazone pattern
* **chalcone**: ``O=C(C=C)c1ccccc1`` - Chalcone scaffold

Rhodanines
~~~~~~~~~~

* **rhodanine**: ``O=C1CSC(=S)N1`` - Rhodanine core
* **thiazolidinone**: ``O=C1CSC(=O)N1`` - Thiazolidinone
* **ene_rhodanine**: ``S=C1SC(=O)N(C)C1=C`` - Ene-rhodanine

Other PAINS
~~~~~~~~~~~

* **barbiturate**: ``O=C1CC(=O)NC(=O)N1`` - Alkylidene barbiturate
* **phenol_hydrazone**: ``c1cc(O)ccc1N=N`` - Hydroxyphenyl hydrazone
* **mannich_phenol**: ``c1cc(O)c(CN)cc1`` - Phenolic Mannich base
* **curcumin**: ``O=C(C=C)C=C(C=C)C(=O)`` - Curcumin-like structure
* **beta_lactam**: ``C1C(=O)NC1`` - Beta-lactam ring
* **styrene**: ``C=Cc1ccccc1`` - Styrene pattern
* **nitro_aromatic**: ``c1ccc([N+](=O)[O-])cc1`` - Nitro aromatic
* **azo**: ``cN=Nc`` - Azo compound
* **thiourea**: ``NC(=S)N`` - Thiourea
* **isothiocyanate**: ``N=C=S`` - Isothiocyanate

Reactive Patterns
-----------------

These patterns identify chemically reactive functional groups that may cause issues in biological assays.

Electrophiles
~~~~~~~~~~~~~

* **acyl_halide**: ``[C;!$(C-[OH])][C;X3](=[OX1])[F,Cl,Br,I]`` - Acyl halide
* **aldehyde**: ``[CX3H1](=O)[#6]`` - Aldehyde group
* **alkyl_halide**: ``[CX4][F,Cl,Br,I]`` - Alkyl halide
* **anhydride**: ``[CX3](=[OX1])[OX2][CX3](=[OX1])`` - Anhydride
* **epoxide**: ``C1OC1`` - Epoxide (oxirane)
* **isocyanate**: ``N=C=O`` - Isocyanate
* **sulfonyl_halide**: ``[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]`` - Sulfonyl halide

Michael Acceptors
~~~~~~~~~~~~~~~~~

* **alpha_beta_unsaturated_carbonyl**: ``C=CC(=O)`` - α,β-unsaturated carbonyl
* **vinyl_sulfone**: ``C=CS(=O)(=O)`` - Vinyl sulfone
* **vinyl_sulfonamide**: ``C=CS(=O)(=O)N`` - Vinyl sulfonamide

Peroxides
~~~~~~~~~

* **peroxide**: ``[OX2][OX2]`` - Peroxide linkage
* **peroxy_acid**: ``C(=O)OO`` - Peroxyacid

Diazo Compounds
~~~~~~~~~~~~~~~

* **diazo**: ``[N-]=[N+]=C`` - Diazo group
* **diazonium**: ``[N+]#N`` - Diazonium salt

Other Reactive Groups
~~~~~~~~~~~~~~~~~~~~~

* **sulfonyl_chloride**: ``S(=O)(=O)Cl`` - Sulfonyl chloride
* **phosphoryl_halide**: ``P(=O)[F,Cl,Br,I]`` - Phosphoryl halide
* **phosphonate**: ``P(=O)(O)(O)`` - Phosphonate
* **hydroxamic_acid**: ``C(=O)N(O)`` - Hydroxamic acid
* **thiol**: ``[SH]`` - Thiol (sulfhydryl)
* **azide**: ``N=[N+]=[N-]`` - Azide
* **crown_ether**: ``C1COCCOCCOCCOC1`` - Crown ether pattern
* **hydrazine**: ``N-N`` - Hydrazine
* **peroxy_acid**: ``C(=O)OO`` - Peroxyacid

Heterocycle Patterns
--------------------

Five-Membered Rings (One Heteroatom)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **pyrrole**: ``c1cc[nH]c1`` - Pyrrole
* **furan**: ``c1ccoc1`` - Furan
* **thiophene**: ``c1ccsc1`` - Thiophene

Five-Membered Rings (Two Heteroatoms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **imidazole**: ``c1cnc[nH]1`` - Imidazole
* **pyrazole**: ``c1cn[nH]c1`` - Pyrazole
* **oxazole**: ``c1cnco1`` - Oxazole
* **isoxazole**: ``c1cnoc1`` - Isoxazole
* **thiazole**: ``c1cncs1`` - Thiazole

Six-Membered Rings (One Heteroatom)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **pyridine**: ``c1ccncc1`` - Pyridine
* **pyran**: ``C1CCOCC1`` - Pyran

Six-Membered Rings (Two Heteroatoms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **pyrimidine**: ``c1cncnc1`` - Pyrimidine
* **pyrazine**: ``c1cnccn1`` - Pyrazine
* **pyridazine**: ``c1cnncc1`` - Pyridazine
* **oxazine**: ``C1COCNC1`` - Oxazine

Six-Membered Rings (Three Heteroatoms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **triazine**: ``c1ncncn1`` - Triazine

Fused Bicyclic Systems
~~~~~~~~~~~~~~~~~~~~~~

* **indole**: ``c1ccc2[nH]ccc2c1`` - Indole
* **benzofuran**: ``c1ccc2occc2c1`` - Benzofuran
* **benzothiophene**: ``c1ccc2sccc2c1`` - Benzothiophene
* **benzimidazole**: ``c1ccc2[nH]cnc2c1`` - Benzimidazole
* **quinoline**: ``c1ccc2ncccc2c1`` - Quinoline
* **isoquinoline**: ``c1ccc2cnccc2c1`` - Isoquinoline
* **quinazoline**: ``c1ccc2ncncc2c1`` - Quinazoline
* **quinoxaline**: ``c1ccc2nccnc2c1`` - Quinoxaline
* **purine**: ``c1nc2[nH]cnc2n1`` - Purine

Saturated Heterocycles
~~~~~~~~~~~~~~~~~~~~~~

* **piperidine**: ``C1CCNCC1`` - Piperidine
* **piperazine**: ``C1CNCCN1`` - Piperazine
* **morpholine**: ``C1COCCN1`` - Morpholine
* **thiomorpholine**: ``C1CSCCN1`` - Thiomorpholine
* **pyrrolidine**: ``C1CCNC1`` - Pyrrolidine

Seven-Membered Rings
~~~~~~~~~~~~~~~~~~~~

* **azepine**: ``C1=CC=CNCC1`` - Azepine
* **oxepine**: ``C1=CC=COCC1`` - Oxepine

Property Ranges
---------------

Common property ranges used in drug discovery:

Lipinski's Rule of Five
~~~~~~~~~~~~~~~~~~~~~~~

* **Molecular Weight**: ≤ 500 Da
* **LogP**: ≤ 5
* **HBD**: ≤ 5
* **HBA**: ≤ 10

Veber's Rules
~~~~~~~~~~~~~

* **Rotatable Bonds**: ≤ 10
* **TPSA**: ≤ 140 Ų

Lead-like Properties
~~~~~~~~~~~~~~~~~~~~

* **Molecular Weight**: 200-350 Da
* **LogP**: 1-3
* **HBD**: ≤ 3
* **HBA**: ≤ 6

Drug-like Properties
~~~~~~~~~~~~~~~~~~~~

* **Molecular Weight**: 200-500 Da
* **LogP**: -1 to 5
* **HBD**: 0-5
* **HBA**: 0-10
* **TPSA**: 20-140 Ų
* **Rotatable Bonds**: ≤ 10
