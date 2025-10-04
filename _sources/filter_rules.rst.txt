Filter Rules Reference
======================

This page documents all the filter patterns used in the library, including SMARTS strings and visual representations of each pattern.

.. note::
   Each filter pattern includes:
   
   * **Pattern name**: The identifier used in the code
   * **SMARTS string**: The chemical pattern definition
   * **Description**: What the pattern identifies
   * **Visual representation**: RDKit-generated structure drawing

   The images shown are generated directly from the SMARTS patterns using RDKit, providing an accurate visual representation of the substructures being matched.

Overview
--------

The library includes four main filter types:

.. list-table:: Filter Types Summary
   :widths: 20 15 65
   :header-rows: 1

   * - Filter Type
     - Pattern Count
     - Description
   * - **PAINS Filter**
     - 19 patterns
     - Identifies Pan-Assay Interference Compounds that frequently show up as false positives in high-throughput screens
   * - **Reactive Filter**
     - 22 patterns
     - Detects chemically reactive functional groups that may cause issues in biological assays or drug development
   * - **Heterocycle Filter**
     - 31 patterns
     - Identifies common heterocyclic scaffolds in medicinal chemistry
   * - **Property Filter**
     - N/A
     - Filters molecules based on calculated properties (MW, LogP, TPSA, HBD, HBA, rotatable bonds, aromatic rings)

Total SMARTS Patterns: **72 patterns** with visual representations

PAINS Patterns
--------------

Pan-Assay Interference Compounds (PAINS) are molecules that frequently show up as hits in high-throughput screens but are likely false positives. These patterns identify problematic structures.

Quinones and Derivatives
~~~~~~~~~~~~~~~~~~~~~~~~

**quinone_A**: ``O=C1C=CC(=O)C=C1``
   p-Benzoquinone pattern

   .. image:: _static/pattern_images/pains/quinone_A.png
      :width: 250px
      :align: center

**quinone_B**: ``O=c1ccc(=O)cc1``
   Alternative quinone representation

   .. image:: _static/pattern_images/pains/quinone_B.png
      :width: 250px
      :align: center

Catechols
~~~~~~~~~

**catechol**: ``c1cc(O)c(O)cc1``
   Catechol (1,2-dihydroxybenzene)

   .. image:: _static/pattern_images/pains/catechol.png
      :width: 250px
      :align: center

**catechol_A**: ``[OH]c1ccccc1[OH]``
   Alternative catechol pattern

   .. image:: _static/pattern_images/pains/catechol_A.png
      :width: 250px
      :align: center

Michael Acceptors
~~~~~~~~~~~~~~~~~

**acyl_hydrazone**: ``C=CC(=O)N``
   Acyl hydrazone pattern

   .. image:: _static/pattern_images/pains/acyl_hydrazone.png
      :width: 250px
      :align: center

**chalcone**: ``O=C(C=C)c1ccccc1``
   Chalcone scaffold

   .. image:: _static/pattern_images/pains/chalcone.png
      :width: 250px
      :align: center

Rhodanines
~~~~~~~~~~

**rhodanine**: ``O=C1CSC(=S)N1``
   Rhodanine core

   .. image:: _static/pattern_images/pains/rhodanine.png
      :width: 250px
      :align: center

**thiazolidinone**: ``O=C1CSC(=O)N1``
   Thiazolidinone

   .. image:: _static/pattern_images/pains/thiazolidinone.png
      :width: 250px
      :align: center

**ene_rhodanine**: ``S=C1SC(=O)N(C)C1=C``
   Ene-rhodanine

   .. image:: _static/pattern_images/pains/ene_rhodanine.png
      :width: 250px
      :align: center

Other PAINS
~~~~~~~~~~~

**barbiturate**: ``O=C1CC(=O)NC(=O)N1``
   Alkylidene barbiturate

   .. image:: _static/pattern_images/pains/barbiturate.png
      :width: 250px
      :align: center

**phenol_hydrazone**: ``c1cc(O)ccc1N=N``
   Hydroxyphenyl hydrazone

   .. image:: _static/pattern_images/pains/phenol_hydrazone.png
      :width: 250px
      :align: center

**mannich_phenol**: ``c1cc(O)c(CN)cc1``
   Phenolic Mannich base

   .. image:: _static/pattern_images/pains/mannich_phenol.png
      :width: 250px
      :align: center

**curcumin**: ``O=C(C=C)C=C(C=C)C(=O)``
   Curcumin-like structure

   .. image:: _static/pattern_images/pains/curcumin.png
      :width: 250px
      :align: center

**beta_lactam**: ``C1C(=O)NC1``
   Beta-lactam ring

   .. image:: _static/pattern_images/pains/beta_lactam.png
      :width: 250px
      :align: center

**styrene**: ``C=Cc1ccccc1``
   Styrene pattern

   .. image:: _static/pattern_images/pains/styrene.png
      :width: 250px
      :align: center

**nitro_aromatic**: ``c1ccc([N+](=O)[O-])cc1``
   Nitro aromatic

   .. image:: _static/pattern_images/pains/nitro_aromatic.png
      :width: 250px
      :align: center

**azo**: ``cN=Nc``
   Azo compound

   .. image:: _static/pattern_images/pains/azo.png
      :width: 250px
      :align: center

**thiourea**: ``NC(=S)N``
   Thiourea

   .. image:: _static/pattern_images/pains/thiourea.png
      :width: 250px
      :align: center

**isothiocyanate**: ``N=C=S``
   Isothiocyanate

   .. image:: _static/pattern_images/pains/isothiocyanate.png
      :width: 250px
      :align: center

Reactive Patterns
-----------------

These patterns identify chemically reactive functional groups that may cause issues in biological assays.

Electrophiles
~~~~~~~~~~~~~

**acyl_halide**: ``[C;!$(C-[OH])][C;X3](=[OX1])[F,Cl,Br,I]``
   Acyl halide

   .. image:: _static/pattern_images/reactive/acyl_halide.png
      :width: 250px
      :align: center

**aldehyde**: ``[CX3H1](=O)[#6]``
   Aldehyde group

   .. image:: _static/pattern_images/reactive/aldehyde.png
      :width: 250px
      :align: center

**alkyl_halide**: ``[CX4][F,Cl,Br,I]``
   Alkyl halide

   .. image:: _static/pattern_images/reactive/alkyl_halide.png
      :width: 250px
      :align: center

**anhydride**: ``[CX3](=[OX1])[OX2][CX3](=[OX1])``
   Anhydride

   .. image:: _static/pattern_images/reactive/anhydride.png
      :width: 250px
      :align: center

**epoxide**: ``C1OC1``
   Epoxide (oxirane)

   .. image:: _static/pattern_images/reactive/epoxide.png
      :width: 250px
      :align: center

**isocyanate**: ``N=C=O``
   Isocyanate

   .. image:: _static/pattern_images/reactive/isocyanate.png
      :width: 250px
      :align: center

**sulfonyl_halide**: ``[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]``
   Sulfonyl halide

   .. image:: _static/pattern_images/reactive/sulfonyl_halide.png
      :width: 250px
      :align: center

Michael Acceptors
~~~~~~~~~~~~~~~~~

**alpha_beta_unsaturated_carbonyl**: ``C=CC(=O)``
   α,β-unsaturated carbonyl

   .. image:: _static/pattern_images/reactive/alpha_beta_unsaturated_carbonyl.png
      :width: 250px
      :align: center

**vinyl_sulfone**: ``C=CS(=O)(=O)``
   Vinyl sulfone

   .. image:: _static/pattern_images/reactive/vinyl_sulfone.png
      :width: 250px
      :align: center

**vinyl_sulfonamide**: ``C=CS(=O)(=O)N``
   Vinyl sulfonamide

   .. image:: _static/pattern_images/reactive/vinyl_sulfonamide.png
      :width: 250px
      :align: center

Peroxides
~~~~~~~~~

**peroxide**: ``[OX2][OX2]``
   Peroxide linkage

   .. image:: _static/pattern_images/reactive/peroxide.png
      :width: 250px
      :align: center

**peroxy_acid**: ``C(=O)OO``
   Peroxyacid

   .. image:: _static/pattern_images/reactive/peroxy_acid.png
      :width: 250px
      :align: center

Diazo Compounds
~~~~~~~~~~~~~~~

**diazo**: ``[N-]=[N+]=C``
   Diazo group

   .. image:: _static/pattern_images/reactive/diazo.png
      :width: 250px
      :align: center

**diazonium**: ``[N+]#N``
   Diazonium salt

   .. image:: _static/pattern_images/reactive/diazonium.png
      :width: 250px
      :align: center

Other Reactive Groups
~~~~~~~~~~~~~~~~~~~~~

**sulfonyl_chloride**: ``S(=O)(=O)Cl``
   Sulfonyl chloride

   .. image:: _static/pattern_images/reactive/sulfonyl_chloride.png
      :width: 250px
      :align: center

**phosphoryl_halide**: ``P(=O)[F,Cl,Br,I]``
   Phosphoryl halide

   .. image:: _static/pattern_images/reactive/phosphoryl_halide.png
      :width: 250px
      :align: center

**phosphonate**: ``P(=O)(O)(O)``
   Phosphonate

   .. image:: _static/pattern_images/reactive/phosphonate.png
      :width: 250px
      :align: center

**hydroxamic_acid**: ``C(=O)N(O)``
   Hydroxamic acid

   .. image:: _static/pattern_images/reactive/hydroxamic_acid.png
      :width: 250px
      :align: center

**thiol**: ``[SH]``
   Thiol (sulfhydryl)

   .. image:: _static/pattern_images/reactive/thiol.png
      :width: 250px
      :align: center

**azide**: ``N=[N+]=[N-]``
   Azide

   .. image:: _static/pattern_images/reactive/azide.png
      :width: 250px
      :align: center

**crown_ether**: ``C1COCCOCCOCCOC1``
   Crown ether pattern

   .. image:: _static/pattern_images/reactive/crown_ether.png
      :width: 250px
      :align: center

**hydrazine**: ``N-N``
   Hydrazine

   .. image:: _static/pattern_images/reactive/hydrazine.png
      :width: 250px
      :align: center

Heterocycle Patterns
--------------------

These patterns identify common heterocyclic scaffolds in medicinal chemistry.

Five-Membered Rings (One Heteroatom)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**pyrrole**: ``c1cc[nH]c1``
   Pyrrole

   .. image:: _static/pattern_images/heterocycle/pyrrole.png
      :width: 250px
      :align: center

**furan**: ``c1ccoc1``
   Furan

   .. image:: _static/pattern_images/heterocycle/furan.png
      :width: 250px
      :align: center

**thiophene**: ``c1ccsc1``
   Thiophene

   .. image:: _static/pattern_images/heterocycle/thiophene.png
      :width: 250px
      :align: center

Five-Membered Rings (Two Heteroatoms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**imidazole**: ``c1cnc[nH]1``
   Imidazole

   .. image:: _static/pattern_images/heterocycle/imidazole.png
      :width: 250px
      :align: center

**pyrazole**: ``c1cn[nH]c1``
   Pyrazole

   .. image:: _static/pattern_images/heterocycle/pyrazole.png
      :width: 250px
      :align: center

**oxazole**: ``c1cnco1``
   Oxazole

   .. image:: _static/pattern_images/heterocycle/oxazole.png
      :width: 250px
      :align: center

**isoxazole**: ``c1cnoc1``
   Isoxazole

   .. image:: _static/pattern_images/heterocycle/isoxazole.png
      :width: 250px
      :align: center

**thiazole**: ``c1cncs1``
   Thiazole

   .. image:: _static/pattern_images/heterocycle/thiazole.png
      :width: 250px
      :align: center

Six-Membered Rings (One Heteroatom)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**pyridine**: ``c1ccncc1``
   Pyridine

   .. image:: _static/pattern_images/heterocycle/pyridine.png
      :width: 250px
      :align: center

**pyran**: ``C1CCOCC1``
   Pyran

   .. image:: _static/pattern_images/heterocycle/pyran.png
      :width: 250px
      :align: center

Six-Membered Rings (Two Heteroatoms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**pyrimidine**: ``c1cncnc1``
   Pyrimidine

   .. image:: _static/pattern_images/heterocycle/pyrimidine.png
      :width: 250px
      :align: center

**pyrazine**: ``c1cnccn1``
   Pyrazine

   .. image:: _static/pattern_images/heterocycle/pyrazine.png
      :width: 250px
      :align: center

**pyridazine**: ``c1cnncc1``
   Pyridazine

   .. image:: _static/pattern_images/heterocycle/pyridazine.png
      :width: 250px
      :align: center

**oxazine**: ``C1COCNC1``
   Oxazine

   .. image:: _static/pattern_images/heterocycle/oxazine.png
      :width: 250px
      :align: center

Six-Membered Rings (Three Heteroatoms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**triazine**: ``c1ncncn1``
   Triazine

   .. image:: _static/pattern_images/heterocycle/triazine.png
      :width: 250px
      :align: center

Fused Bicyclic Systems
~~~~~~~~~~~~~~~~~~~~~~

**indole**: ``c1ccc2[nH]ccc2c1``
   Indole

   .. image:: _static/pattern_images/heterocycle/indole.png
      :width: 250px
      :align: center

**benzofuran**: ``c1ccc2occc2c1``
   Benzofuran

   .. image:: _static/pattern_images/heterocycle/benzofuran.png
      :width: 250px
      :align: center

**benzothiophene**: ``c1ccc2sccc2c1``
   Benzothiophene

   .. image:: _static/pattern_images/heterocycle/benzothiophene.png
      :width: 250px
      :align: center

**benzimidazole**: ``c1ccc2[nH]cnc2c1``
   Benzimidazole

   .. image:: _static/pattern_images/heterocycle/benzimidazole.png
      :width: 250px
      :align: center

**quinoline**: ``c1ccc2ncccc2c1``
   Quinoline

   .. image:: _static/pattern_images/heterocycle/quinoline.png
      :width: 250px
      :align: center

**isoquinoline**: ``c1ccc2cnccc2c1``
   Isoquinoline

   .. image:: _static/pattern_images/heterocycle/isoquinoline.png
      :width: 250px
      :align: center

**quinazoline**: ``c1ccc2ncncc2c1``
   Quinazoline

   .. image:: _static/pattern_images/heterocycle/quinazoline.png
      :width: 250px
      :align: center

**quinoxaline**: ``c1ccc2nccnc2c1``
   Quinoxaline

   .. image:: _static/pattern_images/heterocycle/quinoxaline.png
      :width: 250px
      :align: center

**purine**: ``c1nc2[nH]cnc2n1``
   Purine

   .. image:: _static/pattern_images/heterocycle/purine.png
      :width: 250px
      :align: center

Saturated Heterocycles
~~~~~~~~~~~~~~~~~~~~~~

**piperidine**: ``C1CCNCC1``
   Piperidine

   .. image:: _static/pattern_images/heterocycle/piperidine.png
      :width: 250px
      :align: center

**piperazine**: ``C1CNCCN1``
   Piperazine

   .. image:: _static/pattern_images/heterocycle/piperazine.png
      :width: 250px
      :align: center

**morpholine**: ``C1COCCN1``
   Morpholine

   .. image:: _static/pattern_images/heterocycle/morpholine.png
      :width: 250px
      :align: center

**thiomorpholine**: ``C1CSCCN1``
   Thiomorpholine

   .. image:: _static/pattern_images/heterocycle/thiomorpholine.png
      :width: 250px
      :align: center

**pyrrolidine**: ``C1CCNC1``
   Pyrrolidine

   .. image:: _static/pattern_images/heterocycle/pyrrolidine.png
      :width: 250px
      :align: center

Seven-Membered Rings
~~~~~~~~~~~~~~~~~~~~

**azepine**: ``C1=CC=CNCC1``
   Azepine

   .. image:: _static/pattern_images/heterocycle/azepine.png
      :width: 250px
      :align: center

**oxepine**: ``C1=CC=COCC1``
   Oxepine

   .. image:: _static/pattern_images/heterocycle/oxepine.png
      :width: 250px
      :align: center

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
