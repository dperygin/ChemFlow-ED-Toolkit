# Quick Start Guide: PyMOL + Vina Docking

## Basic Usage

Work from your docking folder (create it the first time):

```bash
mkdir -p ~/Documents/pymol/docking
cd ~/Documents/pymol/docking
```

### 1. Batch Docking a Hit List (3 input files)

`pymol_vina_docking_batchv4.py` needs three files:

| Option | File | Used for |
|---|---|---|
| `--receptor` | Full protein PDB (e.g. `7WX0.pdb` from the PDB) | The structure every ligand is docked and scored against |
| `--pocket` | Pocket PDB from DoGSiteScorer (e.g. `7wx0_P_0_5_res.pdb`) | Only to set the docking box |
| `--sdf` | Hit list from Pharmit (e.g. `query_results.sdf`) | The compounds to dock |

```bash
python3 pymol_vina_docking_batchv4.py \
    --receptor 7WX0/7WX0.pdb \
    --pocket 7WX0/7wx0_P_0_5_res.pdb \
    --sdf 7WX0/query_results.sdf
```

The docking box is computed automatically from the pocket file: it encloses
the pocket's heavy atoms plus 5 Å of padding on each side (change with
`--padding`). You do not need to enter any coordinates.

- The pocket file and the full receptor must come from the same PDB entry so
  they share coordinates. The script checks this and stops with an error if
  the box contains no receptor atoms.
- Do **not** pass the pocket file as `--receptor`. Vina would then dock
  against only the pocket fragment, with the rest of the protein missing.
- To dock a single compound with the same box, use `--ligand ligand.mol2`
  in place of `--sdf`.

The script will:
- Clean and prepare the full receptor (once for the whole batch)
- Compute the docking box from the pocket file
- Split the SDF and prepare each ligand
- Run AutoDock Vina for every compound
- Write a results folder per compound plus `BATCH_SUMMARY.txt`

Results are saved in `docking_results/`:

```
docking_results/
├── BATCH_SUMMARY.txt        # Top 2 poses for every compound
├── receptor/                # Prepared receptor (receptor_clean.pdb, receptor.pdbqt)
├── _sdf_split/              # One SDF per compound from the hit list
└── <compound name>/         # One folder per compound
    ├── docked_poses.pdbqt   # All docked poses
    ├── best_pose.pdb        # Best scoring pose (PDB)
    ├── best_pose.mol2       # Best scoring pose (MOL2)
    ├── docking_session.pse  # PyMOL session file
    ├── docking_result.png   # Visualization image
    ├── vina_log.txt         # Vina output log
    └── DOCKING_SUMMARY.txt  # Summary report for this compound
```

Next steps after docking:
- Rank all compounds: `extract_docking_affinities.py` → `docking_affinities_summary.csv`
- Output a csv file of interactions for a pose: `interactions.py`
- Display key interactions: `displayinteractions.py`
- View one result: `pymol docking_results/<compound name>/docking_session.pse`

### 2. Single Ligand (`pymol_vina_docking.py`)

```bash
python3 pymol_vina_docking.py --receptor protein.pdb --ligand ligand.mol2
```

This single-ligand script builds the box around the entire receptor unless
you give `--center` and `--size` (see below). Its results are written
directly into `docking_results/` with the same per-compound files listed
above.

### 3. Specify Binding Site Manually (single-ligand script)
```bash
python pymol_vina_docking.py \
    --receptor protein.pdb \
    --ligand ligand.mol2 \
    --center 25.5 30.2 15.8 \
    --size 20 20 20 \
    --output my_docking_results
```

### 4. High-Thoroughness Docking
```bash
python pymol_vina_docking.py \
    --receptor protein.pdb \
    --ligand ligand.mol2 \
    --exhaustiveness 32 \
    --num-modes 20
```

### 5. Keep Waters and Cofactors
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
See the `docking_results/` layout under **1. Batch Docking a Hit List** above.

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
python3 pymol_vina_docking_batchv4.py --help
python3 pymol_vina_docking.py --help
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
