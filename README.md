# ChemFlow-ED

**Chem**istry **Flow** for **E**ducation and **D**rug **D**iscovery

ChemFlow-ED is a modular, open-source computational workflow that lets
advanced undergraduate and graduate students carry out authentic
structure-based drug discovery beginning with nothing more than an
experimentally determined protein structure — no known ligand or inhibitor
required. Students identify a candidate binding pocket, build a
receptor-derived pharmacophore, run pharmacophore-guided virtual screening,
execute automated molecular docking, analyze protein–ligand interactions,
and prioritize candidates for future experimental evaluation. A companion
suite of Python scripts automates the repetitive parts of the process —
receptor/ligand preparation, docking, affinity extraction, and interaction
analysis — while leaving every intermediate result open for inspection.

This repository accompanies the manuscript:

> Perygin, D. H. ChemFlow-ED: A Generalized, Open-Source Computational
> Workflow for Teaching Structure-Based Drug Discovery in Advanced
> Medicinal Chemistry. *J. Chem. Educ.* (submitted).

The workflow is illustrated throughout using the canine Cu/Zn superoxide
dismutase E40K mutant (PDB ID: 7WX0) as a representative target, but it is
designed to be transferable to any experimentally determined protein
structure, and is especially well suited to targets that lack a known
small-molecule ligand.

## Repository structure

```
ChemFlow-ED/
├── docs/
│   └── WORKFLOW_GUIDE.md          Full step-by-step workflow (Steps 1-9)
├── scripts/
│   ├── verify_installation.py             Checks required software is installed
│   ├── pymol_vina_docking.py               Single-ligand docking
│   ├── pymol_vina_docking_batchv4.py       Batch docking (SDF hit list -> per-ligand results)
│   ├── docking_batch_with_create_results_file.py   Batch docking, alternate output layout
│   ├── extract_docking_affinities.py       Consolidates Vina affinities into one ranked CSV
│   ├── interactions.py                     Computes H-bond / hydrophobic contacts for a pose
│   ├── displayinteractions.py              Visualizes contacts and saves a PyMOL session
│   ├── pymol_docking_analysis.py           Multi-panel summary plots across a ligand set
│   ├── INSTALLATION_GUIDE.md               Software setup (Vina, PyMOL, MGLTools)
│   └── QUICK_START_GUIDE_batch.md          Condensed usage examples
├── example_data/
│   ├── protein/            Example receptor and pocket structures (7WX0)
│   ├── pharmacophore/      Example pharmacophore workflow notes
│   ├── ligands/            Example virtual-screening hit list (query_results.sdf)
│   └── docking/            Example docking results, organized per ligand
├── paper/                  Manuscript source files and figures
├── supporting_information/ Grading rubric, student assignment, project timeline,
│                           and a redacted example student manuscript
├── requirements.txt
├── LICENSE
├── CITATION.cff
└── RELEASE_INSTRUCTIONS.md
```

## Quick start

1. Install dependencies — see `scripts/INSTALLATION_GUIDE.md` (PyMOL, AutoDock
   Vina, MGLTools) and `requirements.txt` (Python packages).
2. Verify your setup: `python scripts/verify_installation.py`
3. Follow `docs/WORKFLOW_GUIDE.md` end to end, using `example_data/` to see
   expected inputs/outputs at each step, or substitute your own target.

```bash
pip install -r requirements.txt --break-system-packages
python scripts/verify_installation.py
```

## Citing ChemFlow-ED

If you use this workflow or its scripts, please cite both the manuscript and
the software release (see `CITATION.cff`, or the Data and Software
Availability section of the manuscript).

## License

This repository — including scripts, example data, and documentation — is
released under the MIT License. See `LICENSE`.

## Data and Software Availability

All external computational resources used in ChemFlow-ED (Protein Data Bank,
ProteinsPlus/DoGSiteScorer, PharmacoForge, Pharmit, PyMOL Open Source,
AutoDock Vina, MGLTools, Open Babel) are freely available from their
respective developers. The materials in this repository are archived via
Zenodo for long-term accessibility; see `CITATION.cff` for the current
version and DOI.

## Contact

Donna H. Perygin — Department of Chemistry and Geosciences, Jacksonville
State University — dperygin@jsu.edu — [ORCID: 0009-0000-6979-6033](https://orcid.org/0009-0000-6979-6033)
