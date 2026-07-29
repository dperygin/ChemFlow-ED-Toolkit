#!/usr/bin/env python3
"""
extract_docking_affinities.py

Recursively finds all DOCKING_SUMMARY.txt files under a root directory,
extracts ligand names and binding affinities (kcal/mol) for every docking
mode, and writes a combined CSV summary.

Usage:
    python extract_docking_affinities.py [root_dir] [output_csv]

Defaults:
    root_dir   = current working directory (.)
    output_csv = docking_affinities_summary.csv
"""

import os
import re
import csv
import sys

def parse_docking_summary(filepath):
    """
    Parse a single DOCKING_SUMMARY.txt file.
    Returns a dict with:
        ligand_name  : str
        receptor     : str
        modes        : list of dicts {mode, affinity_kcal_mol, rmsd_lb, rmsd_ub}
    """
    result = {
        "ligand_name": "",
        "receptor": "",
        "modes": []
    }

    with open(filepath, "r") as fh:
        lines = fh.readlines()

    # --- extract header fields ---
    for line in lines:
        m = re.match(r"Ligand name\s*:\s*(.+)", line)
        if m:
            result["ligand_name"] = m.group(1).strip()

        m = re.match(r"Receptor file\s*:\s*(.+)", line)
        if m:
            result["receptor"] = m.group(1).strip()

    # --- extract mode table ---
    # The table starts after the separator line "-----+----..."
    in_table = False
    for line in lines:
        if re.match(r"\s*-----\+", line):
            in_table = True
            continue
        if in_table:
            # Each data row: "   1       -7.512          0          0"
            m = re.match(
                r"\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", line
            )
            if m:
                result["modes"].append({
                    "mode":              int(m.group(1)),
                    "affinity_kcal_mol": float(m.group(2)),
                    "rmsd_lb":           float(m.group(3)),
                    "rmsd_ub":           float(m.group(4)),
                })
            else:
                # Stop at first non-matching line after the table begins
                if result["modes"]:
                    break

    return result


def find_summary_files(root_dir):
    matches = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname == "DOCKING_SUMMARY.txt":
                matches.append(os.path.join(dirpath, fname))
    return sorted(matches)


def main():
    root_dir   = sys.argv[1] if len(sys.argv) > 1 else "."
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "docking_affinities_summary.csv"

    summary_files = find_summary_files(root_dir)
    if not summary_files:
        print(f"No DOCKING_SUMMARY.txt files found under: {root_dir}")
        sys.exit(1)

    print(f"Found {len(summary_files)} DOCKING_SUMMARY.txt file(s).")

    rows = []
    for fpath in summary_files:
        subfolder = os.path.relpath(os.path.dirname(fpath), root_dir)
        data = parse_docking_summary(fpath)

        if not data["modes"]:
            print(f"  WARNING: no affinity data found in {fpath}")
            continue

        for mode_info in data["modes"]:
            rows.append({
                "subfolder":         subfolder,
                "ligand_name":       data["ligand_name"] or subfolder,
                "receptor":          data["receptor"],
                "mode":              mode_info["mode"],
                "affinity_kcal_mol": mode_info["affinity_kcal_mol"],
                "rmsd_lb":           mode_info["rmsd_lb"],
                "rmsd_ub":           mode_info["rmsd_ub"],
            })

    # Write CSV
    fieldnames = [
        "subfolder", "ligand_name", "receptor",
        "mode", "affinity_kcal_mol", "rmsd_lb", "rmsd_ub"
    ]
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Results written to: {output_csv}")
    print(f"Total rows (ligand × mode): {len(rows)}")

    # Print a quick best-mode summary
    print("\n--- Best Mode (mode 1) Summary ---")
    print(f"{'Ligand':<50}  {'Affinity (kcal/mol)':>20}")
    print("-" * 74)
    best = [r for r in rows if r["mode"] == 1]
    best.sort(key=lambda r: r["affinity_kcal_mol"])
    for r in best:
        name = r["ligand_name"] if len(r["ligand_name"]) <= 50 else r["ligand_name"][:47] + "..."
        print(f"{name:<50}  {r['affinity_kcal_mol']:>20.3f}")


if __name__ == "__main__":
    main()
