# ChemFlow-ED: Complete Workflow Guide

This guide walks through a full ChemFlow-ED investigation end to end, using the
representative example (PDB 7WX0) as a worked illustration. Replace the PDB ID,
pocket selection, and file names with your own target throughout.

Before starting, run `python scripts/verify_installation.py` to confirm PyMOL,
AutoDock Vina, MGLTools, and required Python packages are installed correctly
(see `scripts/INSTALLATION_GUIDE.md` if anything is missing).

## Step 1 — Select and retrieve a protein structure

Choose a therapeutically relevant target from the Protein Data Bank
(https://www.rcsb.org). Download the structure file (e.g., `7WX0.pdb`).

## Step 2 — Identify candidate binding pockets (DoGSiteScorer)

1. Go to https://proteins.plus.
2. Enter your PDB ID (or upload your own structure file) and click **Go**.
3. Open the **DoGSiteScorer** panel.
4. Under analysis detail, choose **pocket and subpocket**; leave other options
   at their defaults, then click **Calculate**.
5. Review the candidate pockets returned. Each is scored on volume, surface
   area, enclosure, hydrophobicity, and predicted druggability — but scores
   alone should not decide the choice. Inspect each candidate visually in
   PyMOL for accessibility and biological plausibility before selecting one.
6. Click **Download** on your chosen pocket. The `residues/` subfolder of the
   downloaded archive contains the pocket PDB file (e.g., `7wx0_P_0_5_res.pdb`)
   used in the next step.

## Step 3 — Generate a receptor-derived pharmacophore (PharmacoForge)

1. Open the pocket PDB file in PyMOL and run:
   ```
   h_add
   select all
   ```
2. Generate the residue list for PharmacoForge with:
   ```python
   stored.res = set()
   iterate sele, stored.res.add(f"{chain}:{resi}")
   for r in sorted(stored.res): print(r)
   ```
3. Go to https://pharmit.csb.pitt.edu/forge/.
4. Upload the pocket PDB and paste in the residue list (or upload a ligand
   instead, if one is available for your target).
5. Click **Generate Pharmacophore**.
6. Review the returned interaction features. It is common for a small number
   of features to fall outside the pocket or to duplicate one another —
   inspect each and turn off any that are not chemically justified.

## Step 4 — Virtual screening (Pharmit)

1. From PharmacoForge, choose **Send to Pharmit** to search a compound
   database (e.g., MolPort).
2. Click **Search**.
3. If the result set is too large or too small, return to Step 3 and adjust
   individual pharmacophore features (tighten/relax distance tolerances, or
   add/remove a feature) and re-search. This iterative refinement is expected
   and is part of the scientific process, not a sign of error.
4. Once you have a workable hit list, click **Save** to download the results
   as an SDF file (e.g., `query_results.sdf`).

## Step 5 — Automated molecular docking

Place your receptor PDB and the downloaded ligand SDF file in a shared working
directory, then run:

```bash
python3 pymol_vina_docking_batchv4.py --receptor 7WX0/7wx0_P_0_5_res.pdb --sdf 7WX0/query_results.sdf
```

This automatically prepares the receptor and every ligand, runs AutoDock Vina
for each compound, and writes a standardized results directory per ligand
containing docking poses, log files, and a PyMOL session.

> **Note:** AutoDock/PyMOL cannot process some heavy or unusual atoms (e.g.,
> certain metal centers) automatically. If you see an error such as
> `ValueError: Could not find atomic number for ...`, inspect the receptor PDB
> for the offending HETATM record and remove or handle it manually before
> re-running.

## Step 6 — Extract and rank binding affinities

From the directory containing your docking results:

```bash
python3 extract_docking_affinities.py
```

This produces `docking_affinities_summary.csv`, ranking every ligand pose by
predicted binding affinity. AutoDock Vina reports up to nine poses per ligand;
mode 1 is the lowest-energy pose (RMSD = 0 by definition). The same compound
may appear more than once if it was present multiple times in the screening
hit list — each instance represents an independent docking run and may adopt a
different pose, so these can be evaluated separately.

## Step 7 — Interaction analysis

For each candidate pose you want to examine in detail:

1. Open the pose in PyMOL (`docking_session.pse` in that ligand's results
   folder) and disable all other poses so only the receptor and the pose of
   interest remain visible.
2. Copy `interactions.py` and `displayinteractions.py` into the working
   directory, then in the PyMOL command line:
   ```
   h_add
   select all
   run interactions.py
   ```
   This produces `interactions.csv`, summarizing hydrogen bonds, hydrophobic
   contacts, and other neighboring interactions for that pose.
3. Run `displayinteractions.py` to generate a visualization of the interaction
   pattern, and save the image using a name that identifies the ligand.

### Optional publication-quality rendering

```
util.chainbow receptor
select pocket, receptor within 5 of docked_ligand_0001
show surface, pocket
set transparency, 0.6, pocket
color gray80, pocket
util.cbag docked_ligand_0001
set depth_cue, 0
set ray_opaque_background, off
png pymol_pic.png
```

## Step 8 — Compare candidates and prioritize

Consolidate the `interactions.csv` files and the affinity summary across all
evaluated candidates. Prioritize compounds using multiple lines of evidence —
predicted affinity, interaction pattern, pharmacophore satisfaction, and visual
plausibility of the pose — rather than relying on any single score.

## Step 9 — Write up the investigation

Prepare a publication-style report describing the biological rationale,
computational methodology, results, and interpretation, following the
structure of the accompanying manuscript.
