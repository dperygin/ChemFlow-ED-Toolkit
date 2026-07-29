#!/usr/bin/env python3
"""
PyMOL + AutoDock Vina Flexible Docking Script

This script automates flexible ligand docking into a rigid receptor using:
- PyMOL for structure preparation and visualization
- AutoDock Vina for molecular docking

Usage:
    python pymol_vina_docking.py --receptor protein.pdb --ligand ligand.mol2 \
           --center 10.0 20.0 30.0 --size 20 20 20 --output results

Requirements:
    - PyMOL (open-source)
    - AutoDock Vina
    - MGLTools (for PDBQT conversion)
    - Python packages: pymol, biopython, numpy
"""

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path

# Try to import PyMOL
try:
    import pymol
    from pymol import cmd
except ImportError:
    print("ERROR: PyMOL not found. Please install: pip install pymol-open-source")
    sys.exit(1)

# Configuration - EDIT THESE IF NEEDED
VINA_EXECUTABLE = "/usr/local/bin/vina"  # or full path: "/usr/local/bin/vina"
PREPARE_RECEPTOR = "/home/username/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"  # or full path
PREPARE_LIGAND = "/home/username/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"      # or full path
PYTHONSH = "/home/username/MGLTools-1.5.7/bin/pythonsh"  # MGLTools python, or just "python" if conda version


class VinaDocker:
    """Main class for PyMOL + Vina docking workflow"""
    
    def __init__(self, receptor_file, ligand_file, output_dir="docking_results"):
        self.receptor_file = Path(receptor_file)
        self.ligand_file = Path(ligand_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize PyMOL in quiet mode
        pymol.finish_launching(['pymol', '-qc'])
        
    def prepare_receptor(self, remove_waters=True, remove_hetero=True):
        """Prepare receptor using PyMOL and convert to PDBQT"""
        print("=" * 60)
        print("STEP 1: Preparing Receptor")
        print("=" * 60)
        
        # Load receptor
        cmd.load(str(self.receptor_file), "receptor")
        print(f"✓ Loaded receptor: {self.receptor_file}")
        
        # Clean structure
        if remove_waters:
            cmd.remove("resn HOH")
            print("✓ Removed water molecules")
            
        if remove_hetero:
            cmd.remove("hetatm")
            print("✓ Removed heteroatoms")
        
        # Add hydrogens
        cmd.h_add("receptor")
        print("✓ Added hydrogens")
        
        # Save cleaned PDB
        receptor_pdb = self.output_dir / "receptor_clean.pdb"
        cmd.save(str(receptor_pdb), "receptor")
        print(f"✓ Saved cleaned receptor: {receptor_pdb}")
        
        # Convert to PDBQT using MGLTools
        receptor_pdbqt = self.output_dir / "receptor.pdbqt"
        try:
            cmd_list = [
                PYTHONSH, PREPARE_RECEPTOR,
                "-r", str(receptor_pdb),
                "-o", str(receptor_pdbqt),
                "-A", "hydrogens",  # Add hydrogens if not present
                "-U", "nphs"  # Unite non-polar hydrogens
            ]
            subprocess.run(cmd_list, check=True, capture_output=True)
            print(f"✓ Converted to PDBQT: {receptor_pdbqt}")
        except subprocess.CalledProcessError as e:
            print(f"ERROR converting receptor to PDBQT: {e}")
            print(f"STDERR: {e.stderr.decode()}")
            sys.exit(1)
        
        cmd.delete("receptor")
        return receptor_pdbqt
    
    def prepare_ligand(self):
        """Prepare ligand using PyMOL and convert to PDBQT"""
        print("\n" + "=" * 60)
        print("STEP 2: Preparing Ligand")
        print("=" * 60)
        
        # Load ligand
        cmd.load(str(self.ligand_file), "ligand")
        print(f"✓ Loaded ligand: {self.ligand_file}")
        
        # Add hydrogens if needed
        cmd.h_add("ligand")
        print("✓ Added hydrogens")
        
        # Save as PDB
        ligand_pdb = self.output_dir / "ligand.pdb"
        cmd.save(str(ligand_pdb), "ligand")
        print(f"✓ Saved ligand as PDB: {ligand_pdb}")
        
        # Convert to PDBQT using MGLTools
        ligand_pdbqt = self.output_dir / "ligand.pdbqt"
        try:
            cmd_list = [
                PYTHONSH, PREPARE_LIGAND,
                "-l", str(ligand_pdb),
                "-o", str(ligand_pdbqt),
                "-A", "hydrogens"
            ]
            subprocess.run(cmd_list, check=True, capture_output=True)
            print(f"✓ Converted to PDBQT: {ligand_pdbqt}")
        except subprocess.CalledProcessError as e:
            print(f"ERROR converting ligand to PDBQT: {e}")
            print(f"STDERR: {e.stderr.decode()}")
            sys.exit(1)
        
        cmd.delete("ligand")
        return ligand_pdbqt
    
    def calculate_binding_site(self, selection="all", padding=5.0):
        """Calculate binding site center and size from PyMOL selection"""
        print("\n" + "=" * 60)
        print("STEP 3: Calculating Binding Site (Auto-detect)")
        print("=" * 60)
        
        cmd.load(str(self.receptor_file), "receptor")
        
        # Get extent of selection
        extent = cmd.get_extent(selection)
        
        # Calculate center
        center = [
            (extent[0][0] + extent[1][0]) / 2,
            (extent[0][1] + extent[1][1]) / 2,
            (extent[0][2] + extent[1][2]) / 2
        ]
        
        # Calculate size with padding
        size = [
            extent[1][0] - extent[0][0] + padding * 2,
            extent[1][1] - extent[0][1] + padding * 2,
            extent[1][2] - extent[0][2] + padding * 2
        ]
        
        print(f"✓ Binding site center: {center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}")
        print(f"✓ Search box size: {size[0]:.2f} x {size[1]:.2f} x {size[2]:.2f} Å")
        
        cmd.delete("receptor")
        return center, size
    
    def run_vina(self, receptor_pdbqt, ligand_pdbqt, center, size, 
                 exhaustiveness=8, num_modes=9, energy_range=3):
        """Run AutoDock Vina docking"""
        print("\n" + "=" * 60)
        print("STEP 4: Running AutoDock Vina")
        print("=" * 60)
        
        output_pdbqt = self.output_dir / "docked_poses.pdbqt"
        log_file = self.output_dir / "vina_log.txt"
        
        # Construct Vina command
        vina_cmd = [
            VINA_EXECUTABLE,
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--out", str(output_pdbqt),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--energy_range", str(energy_range)
        ]
        
        print(f"Running: {' '.join(vina_cmd)}")
        print("\nDocking in progress... (this may take several minutes)")
        
        try:
            result = subprocess.run(vina_cmd, capture_output=True, text=True, check=True)
            
            # Save log
            with open(log_file, 'w') as f:
                f.write(result.stdout)
            
            print("✓ Docking completed successfully!")
            print(f"\nVina Output:\n{result.stdout}")
            
        except subprocess.CalledProcessError as e:
            print(f"ERROR running Vina: {e}")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            sys.exit(1)
        
        return output_pdbqt, log_file
    
    def visualize_results(self, receptor_pdbqt, docked_pdbqt):
        """Load and visualize docking results in PyMOL"""
        print("\n" + "=" * 60)
        print("STEP 5: Visualizing Results")
        print("=" * 60)
        
        # Clear everything
        cmd.delete("all")
        
        # Load receptor
        cmd.load(str(self.receptor_file), "receptor")
        cmd.color("cyan", "receptor")
        cmd.show("cartoon", "receptor")
        cmd.show("lines", "receptor")
        print("✓ Loaded receptor")
        
        # Load all docked poses
        cmd.load(str(docked_pdbqt), "docked_ligand")
        
        # Split into separate states/models
        cmd.split_states("docked_ligand")
        cmd.delete("docked_ligand")
        
        # Color poses
        num_poses = cmd.count_states("docked_ligand_0001")
        for i in range(1, num_poses + 1):
            obj_name = f"docked_ligand_{i:04d}"
            if i == 1:
                # Best pose in green
                cmd.color("green", obj_name)
                cmd.show("sticks", obj_name)
            else:
                # Other poses in yellow
                cmd.color("yellow", obj_name)
                cmd.show("lines", obj_name)
        
        print(f"✓ Loaded {num_poses} docked poses")
        print("  - Best pose (green): docked_ligand_0001")
        print("  - Other poses (yellow): docked_ligand_0002 to docked_ligand_00XX")
        
        # Zoom to ligand
        cmd.zoom("docked_ligand_0001", buffer=5)
        
        # Save PyMOL session
        session_file = self.output_dir / "docking_session.pse"
        cmd.save(str(session_file))
        print(f"\n✓ Saved PyMOL session: {session_file}")
        
        # Save image
        image_file = self.output_dir / "docking_result.png"
        cmd.png(str(image_file), width=1200, height=1200, dpi=300, ray=1)
        print(f"✓ Saved image: {image_file}")
        
        return session_file
    
    def extract_best_pose(self, docked_pdbqt):
        """Extract the best scoring pose as separate files"""
        print("\n" + "=" * 60)
        print("STEP 6: Extracting Best Pose")
        print("=" * 60)
        
        # Load docked poses
        cmd.delete("all")
        cmd.load(str(docked_pdbqt), "docked")
        cmd.split_states("docked")
        
        # Save best pose
        best_pose_pdb = self.output_dir / "best_pose.pdb"
        best_pose_mol2 = self.output_dir / "best_pose.mol2"
        
        cmd.save(str(best_pose_pdb), "docked_0001")
        cmd.save(str(best_pose_mol2), "docked_0001")
        
        print(f"✓ Saved best pose as PDB: {best_pose_pdb}")
        print(f"✓ Saved best pose as MOL2: {best_pose_mol2}")
        
        return best_pose_pdb
    
    def generate_summary_report(self, log_file):
        """Parse Vina log and create summary report"""
        print("\n" + "=" * 60)
        print("DOCKING SUMMARY")
        print("=" * 60)
        
        with open(log_file, 'r') as f:
            log_content = f.read()
        
        # Extract binding affinities
        if "mode |   affinity" in log_content:
            lines = log_content.split('\n')
            affinity_section = False
            results = []
            
            for line in lines:
                if "mode |   affinity" in line:
                    affinity_section = True
                    continue
                if affinity_section and line.strip().startswith(('-', '+')):
                    break
                if affinity_section and line.strip():
                    parts = line.split()
                    if len(parts) >= 2 and parts[0].isdigit():
                        mode = parts[0]
                        affinity = parts[1]
                        results.append((mode, affinity))
            
            if results:
                print("\nBinding Affinities (kcal/mol):")
                print("-" * 40)
                for mode, affinity in results[:5]:  # Show top 5
                    print(f"  Pose {mode}: {affinity} kcal/mol")
                print(f"\n✓ Best binding affinity: {results[0][1]} kcal/mol (Pose 1)")
        
        # Create summary file
        summary_file = self.output_dir / "DOCKING_SUMMARY.txt"
        with open(summary_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("AutoDock Vina Docking Summary\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Receptor: {self.receptor_file}\n")
            f.write(f"Ligand: {self.ligand_file}\n")
            f.write(f"Output directory: {self.output_dir}\n\n")
            f.write("Results:\n")
            f.write("-" * 60 + "\n")
            if results:
                for mode, affinity in results:
                    f.write(f"Pose {mode}: {affinity} kcal/mol\n")
            f.write("\n" + log_content)
        
        print(f"\n✓ Summary report saved: {summary_file}")
        
        return summary_file


def main():
    parser = argparse.ArgumentParser(
        description="Flexible ligand docking using PyMOL + AutoDock Vina",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect binding site from entire receptor:
  python pymol_vina_docking.py --receptor protein.pdb --ligand ligand.mol2
  
  # Specify binding site manually:
  python pymol_vina_docking.py --receptor protein.pdb --ligand ligand.mol2 \\
         --center 10.0 20.0 30.0 --size 20 20 20
  
  # Keep waters and heteroatoms:
  python pymol_vina_docking.py --receptor protein.pdb --ligand ligand.mol2 \\
         --keep-waters --keep-hetero
        """
    )
    
    parser.add_argument("--receptor", required=True, help="Receptor PDB file")
    parser.add_argument("--ligand", required=True, help="Ligand file (PDB, MOL2, SDF)")
    parser.add_argument("--output", default="docking_results", help="Output directory")
    
    # Binding site options
    parser.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"),
                       help="Center of binding site (x y z). If not provided, auto-detected.")
    parser.add_argument("--size", nargs=3, type=float, metavar=("X", "Y", "Z"),
                       help="Size of search box (x y z). If not provided, auto-detected.")
    parser.add_argument("--padding", type=float, default=5.0,
                       help="Padding around binding site for auto-detection (default: 5.0 Å)")
    
    # Receptor preparation options
    parser.add_argument("--keep-waters", action="store_true",
                       help="Keep water molecules in receptor")
    parser.add_argument("--keep-hetero", action="store_true",
                       help="Keep heteroatoms in receptor")
    
    # Vina options
    parser.add_argument("--exhaustiveness", type=int, default=8,
                       help="Vina exhaustiveness (default: 8, higher = more thorough)")
    parser.add_argument("--num-modes", type=int, default=9,
                       help="Number of binding modes to generate (default: 9)")
    parser.add_argument("--energy-range", type=float, default=3,
                       help="Energy range for poses in kcal/mol (default: 3)")
    
    args = parser.parse_args()
    
    # Check if files exist
    if not Path(args.receptor).exists():
        print(f"ERROR: Receptor file not found: {args.receptor}")
        sys.exit(1)
    if not Path(args.ligand).exists():
        print(f"ERROR: Ligand file not found: {args.ligand}")
        sys.exit(1)
    
    # Check if Vina is available
    try:
        result = subprocess.run([VINA_EXECUTABLE, "--version"], 
                              capture_output=True, text=True, check=True)
        print(f"Found AutoDock Vina: {result.stdout.strip()}\n")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: AutoDock Vina not found!")
        print("Please install Vina and make sure it's in your PATH.")
        print("See INSTALLATION_GUIDE.md for instructions.")
        sys.exit(1)
    
    # Initialize docker
    docker = VinaDocker(args.receptor, args.ligand, args.output)
    
    try:
        # Prepare receptor
        receptor_pdbqt = docker.prepare_receptor(
            remove_waters=not args.keep_waters,
            remove_hetero=not args.keep_hetero
        )
        
        # Prepare ligand
        ligand_pdbqt = docker.prepare_ligand()
        
        # Determine binding site
        if args.center and args.size:
            center = args.center
            size = args.size
            print(f"\n✓ Using user-specified binding site")
            print(f"  Center: {center}")
            print(f"  Size: {size}")
        else:
            center, size = docker.calculate_binding_site(padding=args.padding)
        
        # Run docking
        docked_pdbqt, log_file = docker.run_vina(
            receptor_pdbqt, ligand_pdbqt, center, size,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            energy_range=args.energy_range
        )
        
        # Visualize results
        session_file = docker.visualize_results(receptor_pdbqt, docked_pdbqt)
        
        # Extract best pose
        best_pose = docker.extract_best_pose(docked_pdbqt)
        
        # Generate summary
        summary_file = docker.generate_summary_report(log_file)
        
        print("\n" + "=" * 60)
        print("DOCKING COMPLETE!")
        print("=" * 60)
        print(f"\nAll results saved in: {docker.output_dir}")
        print(f"\nKey files:")
        print(f"  - PyMOL session: {session_file}")
        print(f"  - Best pose: {best_pose}")
        print(f"  - All poses: {docked_pdbqt}")
        print(f"  - Summary: {summary_file}")
        print(f"\nTo view results:")
        print(f"  pymol {session_file}")
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        # Clean up PyMOL
        pymol.cmd.quit()


if __name__ == "__main__":
    main()
