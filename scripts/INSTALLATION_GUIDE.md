# AutoDock Vina + PyMOL Installation Guide

## System Requirements
- Python 3.7 or higher
- PyMOL (open-source educational version)
- AutoDock Vina
- Required Python packages

## Step 1: Install AutoDock Vina

### For Linux:
```bash
# Download Vina
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64
sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina

# Or install via conda (recommended)
conda install -c conda-forge vina
```

### For macOS:
```bash
# Using Homebrew
brew install autodock-vina

# Or using conda
conda install -c conda-forge vina
```

### For Windows:
```bash
# Download from: https://github.com/ccsb-scripps/AutoDock-Vina/releases
# Extract vina.exe to a location in your PATH or note the full path

# Or use conda
conda install -c conda-forge vina
```

## Step 2: Install Required Python Packages
# make sure venv support is installed
sudo apt update
sudo apt install python3-venv python3-full

# create a virtual environment (pick a sensible location)
python3 -m venv ~/docking-env

# activate it
source ~/docking-env/bin/activate

# upgrade pip inside the venv
pip install --upgrade pip

# install your packages
pip install biopython numpy

```bash
or if you don't have PyMOL installed:
pip install biopython numpy pymol-open-source

or if you want it system-wide:
sudo apt install python3-numpy python3-biopython

## Step 3: Install MGLTools (for PDBQT conversion)

### Option A: Use conda (Easiest)

conda install -c conda-forge mgltools
if that doesn't work try: conda install mgltools
or use the installer from ccsb.scripps.edu
if you need conda please see external resources. Here is how I did it:
cd ~
wget https://repo.anaconda.com/archive/Anaconda3-2025.12-2-Linux-x86_64.sh
bash Anaconda3-2025.12-2-Linux-x86_64.sh
source ~/.bashrc (/home/username/anaconda3)
conda --version
✔️ (Recommended) add conda-forge
conda config --add channels conda-forge
conda config --set channel_priority strict
✔️ Create environment if you prefer to use conda and meeko
conda create -n docking python=3.11 numpy biopython
conda activate docking
pip install meeko
✔️ (Optional: MGLTools inside env)
conda create -n mgl python=2.7
conda activate mgl
conda install mgltools

```

### Option B: Manual installation
1. Download from: http://mgltools.scripps.edu/downloads
2. Install following platform-specific instructions
3. Note the installation path (you'll need paths to `prepare_receptor4.py` and `prepare_ligand4.py`)

## Step 4: Verify Installation

```bash
# Check Vina
vina --version

# Check Python packages
python -c "import pymol; import Bio; print('All packages installed successfully')"

# Check MGLTools scripts (adjust paths as needed)
pythonsh /path/to/MGLTools/bin/prepare_receptor4.py --help
pythonsh /path/to/MGLTools/bin/prepare_ligand4.py --help
```

## Step 5: Configure the Docking Script

After installation, you'll need to edit the `pymol_vina_docking.py` script and set:

```python
# If Vina is not in PATH, specify full path
VINA_EXECUTABLE = "vina"  # or "/full/path/to/vina"

# MGLTools script paths (if not using conda version)
PREPARE_RECEPTOR = "prepare_receptor4.py"  # or full path
PREPARE_LIGAND = "prepare_ligand4.py"      # or full path
```

## Troubleshooting

### "vina: command not found"
- Add Vina to your PATH, or use full path in script

### "prepare_receptor4.py: not found"
- Install MGLTools or use conda version
- Set full paths in the script configuration

### PyMOL import errors
- Make sure you're using the same Python environment where PyMOL is installed
- Try: `python -m pip list | grep pymol`

### Permission denied errors (Linux/Mac)
- Make sure Vina is executable: `chmod +x /path/to/vina`

## Quick Test

Once everything is installed, test with:
```bash
python pymol_vina_docking.py --help
```

You should see the help message with all available options.
