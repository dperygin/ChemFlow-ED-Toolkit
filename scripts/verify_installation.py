#!/usr/bin/env python3
"""
Installation verification script for PyMOL + Vina docking setup

This script checks if all required components are properly installed.
"""

import sys
import subprocess
from pathlib import Path

def check_python_version():
    """Check Python version"""
    print("Checking Python version...")
    version = sys.version_info
    if version.major >= 3 and version.minor >= 7:
        print(f"  ✓ Python {version.major}.{version.minor}.{version.micro}")
        return True
    else:
        print(f"  ✗ Python {version.major}.{version.minor}.{version.micro} (need 3.7+)")
        return False

def check_module(module_name, import_name=None):
    """Check if a Python module is installed"""
    if import_name is None:
        import_name = module_name
    
    print(f"Checking {module_name}...")
    try:
        __import__(import_name)
        print(f"  ✓ {module_name} is installed")
        return True
    except ImportError:
        print(f"  ✗ {module_name} is NOT installed")
        print(f"     Install with: pip install {module_name}")
        return False

def check_executable(name, cmd):
    """Check if an executable is available"""
    print(f"Checking {name}...")
    try:
        result = subprocess.run([cmd, "--version"], 
                              capture_output=True, 
                              text=True, 
                              timeout=5)
        if result.returncode == 0:
            version_info = result.stdout.strip() or result.stderr.strip()
            print(f"  ✓ {name} is installed")
            print(f"     {version_info.split(chr(10))[0]}")
            return True
        else:
            print(f"  ✗ {name} not found or error running")
            return False
    except FileNotFoundError:
        print(f"  ✗ {name} is NOT installed or not in PATH")
        return False
    except subprocess.TimeoutExpired:
        print(f"  ✗ {name} command timed out")
        return False

def check_mgltools_script(script_name):
    """Check if MGLTools script is available"""
    print(f"Checking {script_name}...")
    
    # Try different possible locations
    locations = [
        script_name,  # If in PATH
        f"prepare_{script_name}",
    ]
    
    for loc in locations:
        try:
            # Try to run with --help
            result = subprocess.run([loc, "--help"], 
                                  capture_output=True, 
                                  timeout=5)
            if "prepare" in result.stdout.decode().lower() or \
               "prepare" in result.stderr.decode().lower():
                print(f"  ✓ {script_name} found")
                return True
        except:
            continue
    
    # Try pythonsh with full script name
    try:
        result = subprocess.run(["pythonsh", script_name, "--help"],
                              capture_output=True,
                              timeout=5)
        if result.returncode == 0:
            print(f"  ✓ {script_name} found (via pythonsh)")
            return True
    except:
        pass
    
    print(f"  ✗ {script_name} NOT found")
    print(f"     Install MGLTools: conda install -c conda-forge mgltools")
    return False

def main():
    print("=" * 70)
    print("PyMOL + AutoDock Vina Installation Verification")
    print("=" * 70)
    print()
    
    results = {}
    
    # Check Python version
    results['python'] = check_python_version()
    print()
    
    # Check Python modules
    results['pymol'] = check_module('pymol', 'pymol')
    results['biopython'] = check_module('biopython', 'Bio')
    results['numpy'] = check_module('numpy')
    print()
    
    # Check executables
    results['vina'] = check_executable('AutoDock Vina', 'vina')
    print()
    
    # Check MGLTools scripts
    results['prepare_receptor'] = check_mgltools_script('prepare_receptor4.py')
    results['prepare_ligand'] = check_mgltools_script('prepare_ligand4.py')
    print()
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    all_passed = all(results.values())
    
    if all_passed:
        print("\n✓ All components are installed correctly!")
        print("\nYou're ready to run docking:")
        print("  python pymol_vina_docking.py --receptor protein.pdb --ligand ligand.mol2")
    else:
        print("\n✗ Some components are missing. Please install:")
        print()
        
        if not results['python']:
            print("  - Upgrade Python to 3.7 or higher")
        
        if not results['pymol']:
            print("  - PyMOL: pip install pymol-open-source")
        
        if not results['biopython']:
            print("  - BioPython: pip install biopython")
        
        if not results['numpy']:
            print("  - NumPy: pip install numpy")
        
        if not results['vina']:
            print("  - AutoDock Vina:")
            print("      Linux: wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/latest")
            print("      macOS: brew install autodock-vina")
            print("      conda: conda install -c conda-forge vina")
        
        if not results['prepare_receptor'] or not results['prepare_ligand']:
            print("  - MGLTools: conda install -c conda-forge mgltools")
        
        print("\nFor detailed instructions, see INSTALLATION_GUIDE.md")
    
    print()
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())
