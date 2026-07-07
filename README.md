# Computational Post-Processing Pipelines and Functional Annotation Compilers for *Hypania invalida*

This repository contains the version-controlled custom computational source code, heuristic quality control filters, and multi-tool dataset integration engines developed to resolve the nuclear genome repeatome space and transcriptomic coding landscape of the invasive Ponto-Caspian freshwater ampharetid polychaete *Hypania invalida* (Grube, 1860).

## Institutional Affiliation
* **Developer/Author:** Shrikant Vichare
* **Institution:** Bio-Med-Chem Doctoral School, University of Lodz and Lodz Institutes of the Polish Academy of Sciences, Banacha 12/16, 90-237 Lodz, Poland
* **Department:** Department of Invertebrate Zoology and Hydrobiology, University of Lodz

---

## Repository Directory Structure

To guarantee computational reproducibility, the codebase is strictly segregated into structural modules under the `src/` directory tree:

```text
Hypania_invalida_genomic_resources/
├── LICENSE                         <-- Open-source license (MIT)
├── README.md                        <-- Primary documentation homepage (this file)
└── src/
    ├── repeat_curation/
    │   ├── TE_pre_curation_pipeline.py <-- Stage 1: Pre-Curation Translation & Consolidation Matrix
    │   ├── TE_library_curation.py   <-- Stage 2: Housekeeping & Structural Contamination Triage
    │   └── TE_nomenclature.py       <-- Stage 3: Automated Wicker Taxonomy & Suffix Assignation
    └── functional_annotation/
        └── functional_annotation_compiler.py <-- Final Multi-Tool Dataset Synthesis Engine
