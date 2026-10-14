#!/usr/bin/env python3
"""
PyMOL + AutoDock Vina Flexible Docking Script
Supports single-ligand and multi-ligand SDF batch docking.

Usage (batch SDF):
    python pymol_vina_docking_batchv4.py --receptor protein.pdb \
           --pocket pocket.pdb --sdf compounds.sdf

Usage (single ligand):
    python pymol_vina_docking_batchv4.py --receptor protein.pdb \
           --pocket pocket.pdb --ligand ligand.mol2

Three inputs:
    --receptor  the FULL protein structure. Ligands are docked and scored
                against this, so every atom around the pocket is present.
    --pocket    a binding-pocket PDB (e.g. a DoGSiteScorer *_res.pdb file)
                in the same coordinate frame as the receptor. It is used
                ONLY to compute the docking box: the box encloses every
                pocket atom plus --padding on each side.
    --sdf       the hit list (or --ligand for a single compound).

--center and --size can still be given together to set the box by hand,
in which case --pocket is not needed.

Requirements:
    - PyMOL (open-source)
    - AutoDock Vina
    - MGLTools (for PDBQT conversion)
    - Python packages: pymol, numpy
"""

import os
import re
import sys
import argparse
import subprocess
import shutil
from pathlib import Path

# ── Try to import PyMOL ──────────────────────────────────────────────────────
try:
    import pymol
    from pymol import cmd
except ImportError:
    print("ERROR: PyMOL not found.  pip install pymol-open-source")
    sys.exit(1)

# ── Configuration – edit paths to match your installation ───────────────────
VINA_EXECUTABLE  = "/usr/local/bin/vina"
PREPARE_RECEPTOR = (
    "/home/user/MGLTools-1.5.7/MGLToolsPckgs/"
    "AutoDockTools/Utilities24/prepare_receptor4.py"
)
PREPARE_LIGAND   = (
    "/home/user/MGLTools-1.5.7/MGLToolsPckgs/"
    "AutoDockTools/Utilities24/prepare_ligand4.py"
)
PYTHONSH         = "/home/user/MGLTools-1.5.7/bin/pythonsh"

# Absolute working directory where MGLTools reads/writes protein.pdb and
# ligand.pdb.  MGLTools resolves paths relative to its own cwd, so we
# always use absolute paths rooted here.
MGLTOOLS_WORK_DIR = Path("/home/user/Documents/pymol/docking")

_PYMOL_LAUNCHED  = False   # module-level flag – launch PyMOL only once


# ── Utility: split a multi-molecule SDF into individual .sdf files ──────────

def split_sdf_file(sdf_path: Path, work_dir: Path):
    """
    Parse a multi-molecule SDF and write one .sdf file per entry.

    Returns:
        list of (sanitised_name, Path-to-sdf) tuples, in order.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    with open(sdf_path, "r") as fh:
        raw = fh.read()

    blocks = [b.strip() for b in raw.split("$$$$") if b.strip()]

    seen_names = {}
    ligands    = []

    for idx, block in enumerate(blocks, start=1):
        lines      = block.split("\n")
        first_line = lines[0].strip() if lines else ""

        if first_line:
            raw_name = re.sub(r"[^\w\-\.]", "_", first_line)
            raw_name = raw_name[:60]
        else:
            raw_name = f"ligand_{idx:04d}"

        if raw_name in seen_names:
            seen_names[raw_name] += 1
            name = f"{raw_name}_{seen_names[raw_name]}"
        else:
            seen_names[raw_name] = 1
            name = raw_name

        sdf_out = work_dir / f"{name}.sdf"
        with open(sdf_out, "w") as fh:
            fh.write(block + "\n$$$$\n")

        ligands.append((name, sdf_out))

    print(f"  Split SDF into {len(ligands)} individual ligand files -> {work_dir}")
    return ligands


# ── Main docking class ───────────────────────────────────────────────────────

class VinaDocker:
    """
    Manages the prepare -> dock -> visualise -> report pipeline for ONE ligand.

    receptor_file  - original receptor PDB (for visualisation)
    ligand_file    - ligand file (SDF, MOL2, PDB) for this run
    output_dir     - per-ligand result directory
    receptor_pdbqt - already-prepared PDBQT (shared across batch runs)
    """

    def __init__(self, receptor_file, ligand_file, output_dir,
                 receptor_pdbqt=None):
        self.receptor_file  = Path(receptor_file)
        self.ligand_file    = Path(ligand_file)
        self.output_dir     = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._shared_receptor_pdbqt = (
            Path(receptor_pdbqt) if receptor_pdbqt else None
        )

    # ── Receptor preparation ─────────────────────────────────────────────────

    def prepare_receptor(self, remove_waters=True, remove_hetero=True,
                         shared_out_dir=None):
        """
        Clean the receptor in PyMOL and convert to PDBQT.
        If shared_out_dir is given the files are written there so the
        same PDBQT can be reused across all batch ligands.
        """
        out = Path(shared_out_dir) if shared_out_dir else self.output_dir
        out.mkdir(parents=True, exist_ok=True)

        print("=" * 60)
        print("STEP 1: Preparing Receptor")
        print("=" * 60)

        cmd.load(str(self.receptor_file), "receptor")
        print(f"  Loaded receptor: {self.receptor_file}")

        if remove_waters:
            cmd.remove("resn HOH")
            print("  Removed water molecules")
        if remove_hetero:
            cmd.remove("hetatm")
            print("  Removed heteroatoms")

        cmd.h_add("receptor")
        print("  Added hydrogens")

        receptor_pdb   = out / "receptor_clean.pdb"
        receptor_pdbqt = out / "receptor.pdbqt"

        # Also write to the absolute MGLTools work dir so prepare_receptor4.py
        # can find it regardless of the script's cwd.
        MGLTOOLS_WORK_DIR.mkdir(parents=True, exist_ok=True)
        mgl_protein = MGLTOOLS_WORK_DIR / "protein.pdb"
        mgl_receptor_pdbqt = MGLTOOLS_WORK_DIR / "protein.pdbqt"

        cmd.save(str(receptor_pdb), "receptor")   # keep copy in results dir
        cmd.save(str(mgl_protein),  "receptor")   # absolute path for MGLTools
        print(f"  Saved cleaned receptor: {receptor_pdb}")
        print(f"  Saved MGLTools copy   : {mgl_protein}")
        cmd.delete("receptor")

        try:
            subprocess.run(
                [PYTHONSH, PREPARE_RECEPTOR,
                 "-r", str(mgl_protein),
                 "-o", str(mgl_receptor_pdbqt),
                 "-A", "hydrogens",
                 "-U", "nphs"],
                check=True, capture_output=True
            )
            # Copy PDBQT back to the results directory too
            import shutil as _shutil
            _shutil.copy2(str(mgl_receptor_pdbqt), str(receptor_pdbqt))
            print(f"  Converted to PDBQT: {mgl_receptor_pdbqt}")
        except subprocess.CalledProcessError as e:
            print(f"ERROR converting receptor: {e.stderr.decode()}")
            sys.exit(1)

        # Return the absolute path so Vina always gets a resolvable path
        return mgl_receptor_pdbqt

    # ── Ligand preparation ───────────────────────────────────────────────────

    def prepare_ligand(self, work_dir=None):
        """
        Convert ligand to PDBQT via MGLTools.

        work_dir - directory where ligand.pdb and ligand.pdbqt are written.
                   Defaults to self.output_dir.  Pass the main docking_results
                   directory so the files sit alongside the receptor PDBQT
                   and are overwritten each iteration rather than scattered
                   across every per-ligand subdirectory.
        """
        print("\n" + "=" * 60)
        print("STEP 2: Preparing Ligand")
        print("=" * 60)

        out = Path(work_dir) if work_dir else self.output_dir

        cmd.delete("ligand")
        cmd.load(str(self.ligand_file), "ligand")
        print(f"  Loaded ligand: {self.ligand_file}")

        cmd.h_add("ligand")
        print("  Added hydrogens")

        # Write to both the results dir (for record-keeping) and the
        # absolute MGLTools work dir (so prepare_ligand4.py can find it).
        ligand_pdb_results = out / "ligand.pdb"
        MGLTOOLS_WORK_DIR.mkdir(parents=True, exist_ok=True)
        mgl_ligand     = MGLTOOLS_WORK_DIR / "ligand.pdb"
        mgl_ligand_qt  = MGLTOOLS_WORK_DIR / "ligand.pdbqt"
        ligand_pdbqt   = out / "ligand.pdbqt"   # final copy in results dir

        cmd.save(str(ligand_pdb_results), "ligand")  # results dir copy
        cmd.save(str(mgl_ligand),         "ligand")  # absolute path for MGLTools
        print(f"  Saved ligand PDB: {ligand_pdb_results}")
        print(f"  Saved MGLTools copy: {mgl_ligand}")
        cmd.delete("ligand")

        try:
            subprocess.run(
                [PYTHONSH, PREPARE_LIGAND,
                 "-l", str(mgl_ligand),
                 "-o", str(mgl_ligand_qt),
                 "-A", "hydrogens"],
                check=True, capture_output=True
            )
            import shutil as _shutil
            _shutil.copy2(str(mgl_ligand_qt), str(ligand_pdbqt))
            print(f"  Converted to PDBQT: {mgl_ligand_qt}")
        except subprocess.CalledProcessError as e:
            print(f"ERROR converting ligand: {e.stderr.decode()}")
            raise

        # Return absolute path so Vina always gets a resolvable path
        return mgl_ligand_qt

    # ── Binding-site calculation ─────────────────────────────────────────────

    def calculate_binding_site(self, box_file, padding=5.0):
        """
        Compute center + size of the search box from the extent of
        box_file (the binding-pocket PDB), plus `padding` on every side.
        Waters and hydrogens in the pocket file are ignored, so the box is
        the same whether or not the pocket file has been protonated.
        """
        print("\n" + "=" * 60)
        print("STEP 3: Calculating Binding Box from Pocket File")
        print("=" * 60)
        print(f"  Pocket file : {box_file}")

        cmd.load(str(box_file), "pocket_box")
        n_atoms = cmd.count_atoms("pocket_box and not solvent and not hydro")
        if n_atoms == 0:
            cmd.delete("pocket_box")
            raise ValueError(f"No atoms found in pocket file: {box_file}")
        # Measure ONLY the pocket object, not anything else loaded in PyMOL
        extent = cmd.get_extent("pocket_box and not solvent and not hydro")
        cmd.delete("pocket_box")
        print(f"  Pocket atoms: {n_atoms}")

        center = [
            (extent[0][0] + extent[1][0]) / 2,
            (extent[0][1] + extent[1][1]) / 2,
            (extent[0][2] + extent[1][2]) / 2,
        ]
        size = [
            extent[1][0] - extent[0][0] + padding * 2,
            extent[1][1] - extent[0][1] + padding * 2,
            extent[1][2] - extent[0][2] + padding * 2,
        ]

        print(f"  Center : {center[0]:.2f}  {center[1]:.2f}  {center[2]:.2f}")
        print(f"  Box    : {size[0]:.2f} x {size[1]:.2f} x {size[2]:.2f} A")
        return center, size

    # ── Run Vina ─────────────────────────────────────────────────────────────

    def run_vina(self, receptor_pdbqt, ligand_pdbqt, center, size,
                 exhaustiveness=8, num_modes=9, energy_range=3):
        """Run AutoDock Vina for this ligand."""
        print("\n" + "=" * 60)
        print("STEP 4: Running AutoDock Vina")
        print("=" * 60)

        output_pdbqt = self.output_dir / "docked_poses.pdbqt"
        log_file     = self.output_dir / "vina_log.txt"

        vina_cmd = [
            VINA_EXECUTABLE,
            "--receptor",       str(receptor_pdbqt),
            "--ligand",         str(ligand_pdbqt),
            "--center_x",       str(round(center[0], 4)),
            "--center_y",       str(round(center[1], 4)),
            "--center_z",       str(round(center[2], 4)),
            "--size_x",         str(round(size[0], 4)),
            "--size_y",         str(round(size[1], 4)),
            "--size_z",         str(round(size[2], 4)),
            "--out",            str(output_pdbqt),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes",      str(num_modes),
            "--energy_range",   str(energy_range),
        ]

        print(f"CMD: {' '.join(vina_cmd)}")
        print("Docking in progress ...")

        try:
            result = subprocess.run(vina_cmd, capture_output=True, text=True,
                                    check=True)
            with open(log_file, "w") as fh:
                fh.write(result.stdout)
            print("  Docking completed!")
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"ERROR running Vina: {e.stderr}")
            raise

        return output_pdbqt, log_file

    # ── Visualise ────────────────────────────────────────────────────────────

    def visualize_results(self, receptor_pdbqt, docked_pdbqt):
        """Build a PyMOL session showing receptor + all docked poses."""
        print("\n" + "=" * 60)
        print("STEP 5: Visualising Results")
        print("=" * 60)

        cmd.delete("all")
        cmd.load(str(self.receptor_file), "receptor")
        cmd.color("cyan", "receptor")
        cmd.show("cartoon", "receptor")
        cmd.show("lines",   "receptor")
        print("  Loaded receptor")

        cmd.load(str(docked_pdbqt), "docked_ligand")
        cmd.split_states("docked_ligand")
        cmd.delete("docked_ligand")

        try:
            num_poses = cmd.count_states("docked_ligand_0001")
        except Exception:
            num_poses = 1

        for i in range(1, num_poses + 1):
            obj = f"docked_ligand_{i:04d}"
            if i == 1:
                cmd.color("green", obj)
                cmd.show("sticks", obj)
            else:
                cmd.color("yellow", obj)
                cmd.show("lines", obj)

        print(f"  Loaded {num_poses} docked pose(s)")
        cmd.zoom("docked_ligand_0001", buffer=5)

        session_file = self.output_dir / "docking_session.pse"
        image_file   = self.output_dir / "docking_result.png"
        cmd.save(str(session_file))
        cmd.png(str(image_file), width=1200, height=1200, dpi=300, ray=1)
        print(f"  Session : {session_file}")
        print(f"  Image   : {image_file}")

        return session_file

    # ── Extract best pose ────────────────────────────────────────────────────

    def extract_best_pose(self, docked_pdbqt):
        """Save pose 1 (best) as PDB and MOL2."""
        print("\n" + "=" * 60)
        print("STEP 6: Extracting Best Pose")
        print("=" * 60)

        cmd.delete("docked_best")
        cmd.load(str(docked_pdbqt), "docked_best")
        cmd.split_states("docked_best")
        cmd.delete("docked_best")

        best_pose_pdb  = self.output_dir / "best_pose.pdb"
        best_pose_mol2 = self.output_dir / "best_pose.mol2"

        cmd.save(str(best_pose_pdb),  "docked_best_0001")
        cmd.save(str(best_pose_mol2), "docked_best_0001")
        print(f"  Best pose PDB  : {best_pose_pdb}")
        print(f"  Best pose MOL2 : {best_pose_mol2}")

        return best_pose_pdb

    # ── Per-ligand summary report ────────────────────────────────────────────

    def generate_summary_report(self, log_file, ligand_name=""):
        """
        Parse Vina log -> DOCKING_SUMMARY.txt.

        Returns:
            (summary_file, top2_rows)
            top2_rows is a list of up to 2 tuples:
                (mode, affinity, rmsd_lb, rmsd_ub)
        """
        print("\n" + "=" * 60)
        print("DOCKING SUMMARY")
        print("=" * 60)

        with open(log_file, "r") as fh:
            log_content = fh.read()

        # Parse all 4 columns Vina reports: mode, affinity, rmsd_lb, rmsd_ub
        # The score rows sit BELOW the "-----+------------+" divider that
        # follows the column headers, so start collecting after the divider
        # and stop at the first non-score line once rows have been read.
        results = []
        if "mode |   affinity" in log_content:
            in_table = False
            for line in log_content.split("\n"):
                if line.strip().startswith("-----+"):
                    in_table = True
                    continue
                if not in_table:
                    continue
                parts = line.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    results.append((parts[0], parts[1], parts[2], parts[3]))
                elif len(parts) == 2 and parts[0].isdigit():
                    # Fallback if rmsd columns are absent
                    results.append((parts[0], parts[1], "N/A", "N/A"))
                elif results:
                    break

            if results:
                print("\nBinding Affinities (kcal/mol):")
                for mode, aff, *_ in results[:5]:
                    print(f"  Pose {mode}: {aff} kcal/mol")
                print(f"\n  Best affinity: {results[0][1]} kcal/mol (Pose 1)")

        summary_file = self.output_dir / "DOCKING_SUMMARY.txt"
        with open(summary_file, "w") as fh:
            fh.write("=" * 60 + "\n")
            fh.write("AutoDock Vina Docking Summary\n")
            fh.write("=" * 60 + "\n\n")
            fh.write(
                f"Ligand name      : {ligand_name or self.ligand_file.name}\n"
            )
            fh.write(f"Ligand file      : {self.ligand_file}\n")
            fh.write(f"Receptor file    : {self.receptor_file}\n")
            fh.write(f"Output directory : {self.output_dir}\n\n")
            fh.write("Results:\n")
            fh.write("-" * 60 + "\n")
            for mode, aff, rmsd_lb, rmsd_ub in results:
                fh.write(
                    f"Pose {mode}: {aff} kcal/mol  "
                    f"rmsd_lb={rmsd_lb}  rmsd_ub={rmsd_ub}\n"
                )
            fh.write("\n" + log_content)

        print(f"  Summary saved : {summary_file}")
        # Return only the top 2 poses for the batch summary table
        return summary_file, results[:2]


# ── Shared receptor preparation helper ──────────────────────────────────────

def prepare_shared_receptor(receptor_file, base_output_dir,
                             remove_waters=True, remove_hetero=True):
    """
    Prepare the receptor PDBQT once and store it in
    <base_output_dir>/receptor/.
    Returns Path to receptor.pdbqt.
    """
    shared_dir = Path(base_output_dir) / "receptor"
    shared_dir.mkdir(parents=True, exist_ok=True)

    dummy = VinaDocker(receptor_file, receptor_file, shared_dir)
    return dummy.prepare_receptor(
        remove_waters=remove_waters,
        remove_hetero=remove_hetero,
        shared_out_dir=shared_dir,
    )


def check_box_against_receptor(receptor_file, center, size):
    """
    Count full-receptor atoms inside the docking box. Zero means the pocket
    file and the receptor are not in the same coordinate frame (e.g. the
    pocket came from a different PDB entry or a re-centred structure).
    """
    cmd.load(str(receptor_file), "receptor_check")
    coords = cmd.get_coords("receptor_check and not solvent")
    cmd.delete("receptor_check")
    lo = [center[i] - size[i] / 2 for i in range(3)]
    hi = [center[i] + size[i] / 2 for i in range(3)]
    inside = 0
    if coords is not None:
        for x, y, z in coords:
            if lo[0] <= x <= hi[0] and lo[1] <= y <= hi[1] and lo[2] <= z <= hi[2]:
                inside += 1
    print(f"  Receptor atoms inside box: {inside}")
    if inside == 0:
        print("ERROR: The docking box contains no receptor atoms. Check that the "
              "pocket file and the receptor come from the same structure and "
              "coordinate frame.")
        sys.exit(1)


# ── Batch summary table ──────────────────────────────────────────────────────

def write_batch_summary(batch_pose_rows, base_output_dir):
    """
    Write docking_results/BATCH_SUMMARY.txt

    Format mirrors Vina's own results table with a molecule-name column
    prepended.  Only the top 2 poses per ligand are included.

    Example:
        molecule name                          | mode |   affinity | dist from best mode
                                               |      | (kcal/mol) | rmsd l.b.| rmsd u.b.
        ---------------------------------------+------+------------+----------+----------
        MolPort-002-793-101                    |    1 |     -8.300 |    0.000 |    0.000
                                               |    2 |     -8.100 |    1.256 |    2.341

    batch_pose_rows: list of (lig_name, top2_poses, status)
        top2_poses: list of up to 2 tuples (mode, affinity, rmsd_lb, rmsd_ub)
    """
    NAME_W = 45   # width of the molecule-name column

    header1 = (
        f"{'molecule name':<{NAME_W}} | mode |   affinity | dist from best mode"
    )
    header2 = (
        f"{'':<{NAME_W}} |      | (kcal/mol) | rmsd l.b.| rmsd u.b."
    )
    divider = "-" * NAME_W + "-+------+------------+----------+----------"

    out_path = Path(base_output_dir) / "BATCH_SUMMARY.txt"
    with open(out_path, "w") as fh:
        fh.write(header1 + "\n")
        fh.write(header2 + "\n")
        fh.write(divider + "\n")

        for lig_name, top2_poses, status in batch_pose_rows:
            if status != "OK" or not top2_poses:
                fh.write(
                    f"{lig_name:<{NAME_W}} | {'':>4} | "
                    f"{'FAILED / NO RESULT':>10} | {'':>8} | {'':>8}\n"
                )
                continue

            for row_idx, (mode, aff, rmsd_lb, rmsd_ub) in enumerate(top2_poses):
                # Only the first row shows the molecule name; second row is blank
                name_col = lig_name if row_idx == 0 else ""
                fh.write(
                    f"{name_col:<{NAME_W}} | {mode:>4} | {aff:>10} "
                    f"| {rmsd_lb:>8} | {rmsd_ub:>8}\n"
                )

    print(f"\n  Batch summary : {out_path}")
    return out_path


# ── Argument parsing ─────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Flexible ligand docking - PyMOL + AutoDock Vina",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Batch SDF, box computed from the pocket file:
  python pymol_vina_docking_batchv4.py --receptor 7wx0.pdb \\
         --pocket 7wx0_P_1_res.pdb --sdf query_results.sdf

  # Single ligand:
  python pymol_vina_docking_batchv4.py --receptor 7wx0.pdb \\
         --pocket 7wx0_P_1_res.pdb --ligand ligand.mol2

  # Set the box by hand instead (no pocket file needed):
  python pymol_vina_docking_batchv4.py --receptor 7wx0.pdb --sdf compounds.sdf \\
         --center 14.5 -22.0 104.0 --size 30 30 30
        """
    )

    parser.add_argument("--receptor", required=True,
                        help="FULL receptor PDB file (docked against)")
    parser.add_argument("--pocket",
                        help="Binding-pocket PDB, same coordinate frame as the "
                             "receptor; used only to compute the docking box")
    lig_group = parser.add_mutually_exclusive_group(required=True)
    lig_group.add_argument("--ligand",
                           help="Single ligand file (PDB, MOL2, SDF)")
    lig_group.add_argument("--sdf",
                           help="Multi-molecule SDF file for batch docking")

    parser.add_argument("--output", default="docking_results",
                        help="Base output directory (default: docking_results)")

    parser.add_argument("--center", nargs=3, type=float,
                        metavar=("X", "Y", "Z"),
                        help="Manual box center (use with --size instead of --pocket)")
    parser.add_argument("--size", nargs=3, type=float,
                        metavar=("X", "Y", "Z"),
                        help="Manual box size in Angstroms (use with --center)")
    parser.add_argument("--padding", type=float, default=5.0,
                        help="Padding around the pocket for the auto box (default: 5.0 A)")

    parser.add_argument("--keep-waters",  action="store_true",
                        help="Keep water molecules")
    parser.add_argument("--keep-hetero",  action="store_true",
                        help="Keep heteroatoms")

    parser.add_argument("--exhaustiveness", type=int, default=8,
                        help="Vina exhaustiveness (default: 8)")
    parser.add_argument("--num-modes",  type=int, default=9,
                        help="Number of binding modes (default: 9)")
    parser.add_argument("--energy-range", type=float, default=3,
                        help="Energy range in kcal/mol (default: 3)")

    parser.add_argument("--skip-viz", action="store_true",
                        help="Skip PyMOL visualisation (faster for large batches)")

    return parser.parse_args()


# ── Entry point ──────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if not Path(args.receptor).exists():
        print(f"ERROR: Receptor not found: {args.receptor}")
        sys.exit(1)
    manual_box = args.center is not None or args.size is not None
    if manual_box and (args.center is None or args.size is None):
        print("ERROR: --center and --size must be given together.")
        sys.exit(1)
    if not manual_box and not args.pocket:
        print("ERROR: Provide --pocket (pocket PDB used to set the docking box), "
              "or set the box by hand with --center and --size.")
        sys.exit(1)
    if args.pocket and not Path(args.pocket).exists():
        print(f"ERROR: Pocket file not found: {args.pocket}")
        sys.exit(1)
    if args.ligand and not Path(args.ligand).exists():
        print(f"ERROR: Ligand not found: {args.ligand}")
        sys.exit(1)
    if args.sdf and not Path(args.sdf).exists():
        print(f"ERROR: SDF not found: {args.sdf}")
        sys.exit(1)

    try:
        result = subprocess.run([VINA_EXECUTABLE, "--version"],
                                capture_output=True, text=True, check=True)
        print(f"Found AutoDock Vina: {result.stdout.strip()}\n")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: AutoDock Vina not found! Check VINA_EXECUTABLE path.")
        sys.exit(1)

    global _PYMOL_LAUNCHED
    if not _PYMOL_LAUNCHED:
        pymol.finish_launching(["pymol", "-qc"])
        _PYMOL_LAUNCHED = True

    base_output = Path(args.output)
    base_output.mkdir(parents=True, exist_ok=True)

    # Prepare receptor ONCE
    print("\n" + "=" * 60)
    print("RECEPTOR PREPARATION  (done once for entire batch)")
    print("=" * 60)
    receptor_pdbqt = prepare_shared_receptor(
        args.receptor,
        base_output,
        remove_waters=not args.keep_waters,
        remove_hetero=not args.keep_hetero,
    )

    # Determine binding box ONCE
    if manual_box:
        center = args.center
        size   = args.size
        print(f"\n  Using manual binding box")
        print(f"  Center: {center}")
        print(f"  Size  : {size}")
        if args.pocket:
            print("  (--pocket ignored because --center/--size were given)")
    else:
        _tmp = VinaDocker(args.receptor, args.receptor, base_output,
                          receptor_pdbqt=receptor_pdbqt)
        center, size = _tmp.calculate_binding_site(args.pocket,
                                                   padding=args.padding)
    check_box_against_receptor(args.receptor, center, size)

    # Build list of (name, path) for all ligands
    if args.sdf:
        sdf_work_dir = base_output / "_sdf_split"
        ligand_list  = split_sdf_file(Path(args.sdf), sdf_work_dir)
    else:
        stem = Path(args.ligand).stem
        ligand_list = [(stem, Path(args.ligand))]

    total = len(ligand_list)
    print(f"\n{'=' * 60}")
    print(f"BATCH DOCKING: {total} ligand(s)")
    print(f"{'=' * 60}\n")

    # batch_pose_rows: list of (lig_name, top2_poses, status)
    batch_pose_rows = []

    for i, (lig_name, lig_path) in enumerate(ligand_list, start=1):
        print(f"\n{'#' * 60}")
        print(f"  Ligand {i}/{total}: {lig_name}")
        print(f"{'#' * 60}")

        lig_out = base_output / lig_name

        docker = VinaDocker(
            receptor_file  = args.receptor,
            ligand_file    = lig_path,
            output_dir     = lig_out,
            receptor_pdbqt = receptor_pdbqt,
        )

        try:
            ligand_pdbqt = docker.prepare_ligand(work_dir=base_output)

            docked_pdbqt, log_file = docker.run_vina(
                receptor_pdbqt, ligand_pdbqt, center, size,
                exhaustiveness = args.exhaustiveness,
                num_modes      = args.num_modes,
                energy_range   = args.energy_range,
            )

            if not args.skip_viz:
                docker.visualize_results(receptor_pdbqt, docked_pdbqt)
            else:
                print("(Visualisation skipped)")

            docker.extract_best_pose(docked_pdbqt)

            _, top2_poses = docker.generate_summary_report(
                log_file, ligand_name=lig_name
            )

            batch_pose_rows.append((lig_name, top2_poses, "OK"))

        except Exception as exc:
            print(f"\n  WARNING: FAILED for {lig_name}: {exc}")
            batch_pose_rows.append((lig_name, [], f"FAILED: {exc}"))

    # Write BATCH_SUMMARY.txt
    write_batch_summary(batch_pose_rows, base_output)

    # Console final summary
    print("\n" + "=" * 60)
    print("ALL DOCKING RUNS COMPLETE")
    print("=" * 60)

    ok_runs   = [r for r in batch_pose_rows if r[2] == "OK"]
    fail_runs = [r for r in batch_pose_rows if r[2] != "OK"]

    print(f"\n  Successful : {len(ok_runs)}/{total}")
    if fail_runs:
        print(f"  Failed     : {len(fail_runs)}/{total}")
        for name, _, status in fail_runs:
            print(f"    * {name} -> {status}")

    if ok_runs:
        try:
            sorted_ok = sorted(
                ok_runs,
                key=lambda r: float(r[1][0][1]) if r[1] else 0
            )
            print(f"\nTop 5 hits by binding affinity:")
            for name, poses, _ in sorted_ok[:5]:
                aff = poses[0][1] if poses else "N/A"
                print(f"  {name:<50s}  {aff} kcal/mol")
        except (ValueError, IndexError):
            pass

    print(f"\nAll results in  : {base_output}/")
    print(f"Batch summary   : {base_output}/BATCH_SUMMARY.txt")

    pymol.cmd.quit()


if __name__ == "__main__":
    main()
