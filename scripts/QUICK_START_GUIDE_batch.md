# Quick Start Guide: PyMOL + Vina Docking

## Basic Usage

### 1. Simplest Usage (Auto-detect binding site)


```bash

cd ~/Documents/pymol
mkdir docking
cd docking
or if exists
cd ~/Documents/pymol/docking
example usage:
python pymol_vina_docking.py --receptor 5NN3/protein.pdb --ligand 5NN3/J059-0151_mol.sdf
python pymol_vina_docking.py --receptor protein.pdb --ligand ligand.mol2
python pymol_vina_docking.py --receptor 3HS4.pdb --ligand benzamidine.pdb
To dock an entire sdf of ligands (3 input files: full receptor, pocket, sdf):
python pymol_vina_docking_batchv4.py --receptor 1KW0/1kw0.pdb --pocket 1KW0/1kw0_P_0_res.pdb --sdf 1KW0/query_results.sdf



```
All results saved in: docking_results

Key files:
  - PyMOL session: docking_results/docking_session.pse
  - Best pose: docking_results/best_pose.pdb
  - All poses: docking_results/docked_poses.pdbqt
  - Summary: docking_results/DOCKING_SUMMARY.txt
  - To find the best overall candidates: extract_docking_affinities.py
  - Analyze best candidates: docking_affinities_summary.csv
  - Output a csv file of interactions: interactions.py
  - Display key interactions: displayinteractions.py
  
  

To view results individually:
  pymol docking_results/docking_session.pse

This will:
- Clean and prepare your receptor
- Prepare your ligand
- Auto-detect the binding site (single-ligand script: from the receptor
  structure; batch script: from the pocket file — see below)
- Run docking with default parameters
- Generate visualizations and results

**Batch script inputs (`pymol_vina_docking_batchv4.py`):**
- `--receptor` — the **full** protein PDB. Every ligand is docked and scored
  against this structure.
- `--pocket` — the binding-pocket PDB (e.g. the DoGSiteScorer `*_res.pdb`
  file). It is used **only** to set the docking box, which encloses the
  pocket's heavy atoms plus 5 Å of padding on each side (change with
  `--padding`).
- `--sdf` — the hit list to dock (or `--ligand` for a single compound).

The pocket file and the full receptor must come from the same PDB entry so
they share coordinates. The script checks this and stops with an error if the
box contains no receptor atoms. Do **not** pass the pocket file as
`--receptor`: Vina would then dock against only the pocket fragment, with the
rest of the protein missing.

### 2. Specify Binding Site Manually
```bash
python pymol_vina_docking.py \
    --receptor protein.pdb \
    --ligand ligand.mol2 \
    --center 25.5 30.2 15.8 \
    --size 20 20 20 \
    --output my_docking_results
```

### 3. High-Thoroughness Docking
```bash
python pymol_vina_docking.py \
    --receptor protein.pdb \
    --ligand ligand.mol2 \
    --exhaustiveness 32 \
    --num-modes 20
```

### 4. Keep Waters and Cofactors
```bash
python pymol_vina_docking.py \
    --receptor protein.pdb \
    --ligand ligand.mol2 \
    --keep-waters \
    --keep-hetero
```

## Finding Your Binding Site Coordinates

If you don't know the binding site coordinates, you have several options:

### Option 1: Let the script auto-detect (easiest)
Just don't specify `--center` and `--size`. The single-ligand script
(`pymol_vina_docking.py`) will use the entire receptor; the batch script
(`pymol_vina_docking_batchv4.py`) will build the box around the `--pocket` file.
With the batch script, giving `--center` and `--size` together overrides the
pocket file (both are required; giving only one is an error).

### Option 2: Use PyMOL to find coordinates manually

1. Open your receptor in PyMOL:
```bash
pymol protein.pdb
```

2. Select the binding site region (e.g., residues around a known ligand):
```
select binding_site, resi 100-120
```

3. Get the center:
```
get_extent binding_site
```

4. PyMOL will show you the min/max coordinates. Calculate center:
- center_x = (min_x + max_x) / 2
- center_y = (min_y + max_y) / 2
- center_z = (min_z + max_z) / 2

5. Use a box size that covers the region with some padding (typically 20-30 Å per side)

### Option 3: If you have a reference ligand in your PDB

The script can be modified to use a reference ligand. Let me know if you need this!

## Understanding the Output

After running, you'll find in your output directory:

### Files Generated:
```
docking_results/
├── receptor_clean.pdb       # Cleaned receptor structure
├── receptor.pdbqt           # Receptor in PDBQT format
├── ligand.pdb               # Ligand in PDB format
├── ligand.pdbqt             # Ligand in PDBQT format
├── docked_poses.pdbqt       # All docked poses
├── best_pose.pdb            # Best scoring pose (PDB)
├── best_pose.mol2           # Best scoring pose (MOL2)
├── docking_session.pse      # PyMOL session file
├── docking_result.png       # Visualization image
├── vina_log.txt             # Vina output log
└── DOCKING_SUMMARY.txt      # Summary report
```

### Interpreting Results:

**Binding Affinity Values:**
- More negative = stronger binding
- Typical ranges:
  - < -10 kcal/mol: Very strong binding
  - -8 to -10 kcal/mol: Strong binding
  - -6 to -8 kcal/mol: Moderate binding
  - > -6 kcal/mol: Weak binding

**Multiple Poses:**
- Pose 1 is always the best scoring
- Other poses represent alternative binding modes
- Similar affinities suggest flexible binding

## Viewing Results

### Open in PyMOL:
```bash
pymol docking_results/docking_session.pse
```

In the session:
- **Receptor**: cyan, cartoon + lines
- **Best pose**: green, sticks (docked_ligand_0001)
- **Alternative poses**: yellow, lines (docked_ligand_0002, etc.)

### PyMOL commands for analysis:
```python
# Show only best pose
hide everything, docked_ligand_*
show sticks, docked_ligand_0001

# Measure distances
distance myDist, receptor and resi 100 and name CA, docked_ligand_0001

# Show hydrogen bonds (if receptor has H atoms)
distance HBonds, receptor, docked_ligand_0001, mode=2

# Color by element
util.cbag docked_ligand_0001

# Show surface
show surface, receptor
set transparency, 0.5
```

## Common Issues & Solutions

### Issue: "vina: command not found"
**Solution:** Make sure Vina is installed and in PATH. Try:
```bash
which vina
# or
conda install -c conda-forge vina
```

### Issue: "prepare_receptor4.py: not found"
**Solution:** Install MGLTools:
```bash
conda install -c conda-forge mgltools
```

### Issue: Poor docking results
**Solutions:**
1. Increase exhaustiveness (try 16 or 32)
2. Verify binding site is correct
3. Check if your ligand preparation is correct
4. Try different box sizes

### Issue: Script too slow
**Solutions:**
1. Decrease exhaustiveness (try 4 for quick testing)
2. Reduce num_modes
3. Make search box smaller

## Advanced: Batch Docking

The easiest way to dock many ligands is `pymol_vina_docking_batchv4.py`,
which takes a multi-compound SDF directly:

```bash
python pymol_vina_docking_batchv4.py \
    --receptor protein.pdb \
    --pocket pocket_res.pdb \
    --sdf query_results.sdf
```

Alternatively, to loop the single-ligand script over separate files, create a
simple bash script:

```bash
#!/bin/bash
# batch_dock.sh

for ligand in ligands/*.mol2; do
    name=$(basename "$ligand" .mol2)
    echo "Docking $name..."
    python pymol_vina_docking.py \
        --receptor protein.pdb \
        --ligand "$ligand" \
        --center 25.5 30.2 15.8 \
        --size 20 20 20 \
        --output "results_$name"
done
```

Run with:
```bash
chmod +x batch_dock.sh
./batch_dock.sh
```

## Customization

The script is designed to be easily modified. Common customizations:

### 1. Change Vina parameters:
Edit line ~320 in the script to add more Vina options:
```python
vina_cmd = [
    VINA_EXECUTABLE,
    # ... existing options ...
    "--cpu", "4",  # Use 4 CPU cores
    "--seed", "12345",  # For reproducibility
]
```

### 2. Different receptor preparation:
Modify the `prepare_receptor()` method to customize cleaning.

### 3. Add more visualization:
Modify `visualize_results()` to customize PyMOL rendering.

## Getting Help

Run with `--help` to see all options:
```bash
python pymol_vina_docking.py --help
```

## Example with Real Data

Here's a complete example using a test system:

```bash
# 1. Download example files (HIV protease + ligand)
wget https://files.rcsb.org/download/1HVR.pdb -O protein.pdb
wget https://files.rcsb.org/ligands/view/XK2.sdf -O ligand.sdf

# 2. Convert SDF to MOL2 using Open Babel (if needed)
obabel ligand.sdf -O ligand.mol2

# 3. Run docking
python pymol_vina_docking.py \
    --receptor protein.pdb \
    --ligand ligand.mol2 \
    --center 15.0 18.0 24.0 \
    --size 25 25 25 \
    --output hiv_protease_docking

# 4. View results
pymol hiv_protease_docking/docking_session.pse
```

## Tips for Best Results

1. **Always check your binding site** - incorrect box position is the #1 cause of poor results
2. **Start with default parameters** - only increase exhaustiveness if needed
3. **Prepare your structures carefully** - remove unnecessary molecules
4. **Validate results** - check if poses make chemical sense
5. **Use multiple starting conformations** - if ligand is very flexible

## Next Steps

Once you're comfortable with basic docking:
- Try rescoring with different force fields
- Perform MD simulations on best poses
- Calculate binding free energies with MM-PBSA
- Dock libraries of compounds for virtual screening
