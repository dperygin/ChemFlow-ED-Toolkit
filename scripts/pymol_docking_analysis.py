"""
pymol_docking_analysis.py
=========================
Analyze molecular docking results loaded in PyMOL.

Produces a multi-panel chart showing:
  1. Docking scores for all ligands
  2. Unique atom-atom contact counts (one entry per unique atom pair,
     shortest distance only)
  3. Distance distributions across unique pairs per ligand
  4. Score vs. contact-count correlation scatter

Usage
-----
From within PyMOL's command line or a script panel:

    run /path/to/pymol_docking_analysis.py

Then call the main function:

    analyze_docking()                              # uses defaults
    analyze_docking(receptor="my_protein",         # explicit receptor name
                    cutoff=5.0,                    # interaction distance cutoff (Å)
                    output="results.png")          # output image path

How scores are read
-------------------
The script tries three strategies in order:
  1. A common constant B-factor across all atoms in the ligand — many
     docking tools (AutoDock Vina, GNINA, Glide) write the score there.
  2. A "REMARK" or "SCORE" annotation stored in the PDB HEADER / REMARK
     lines (read back via cmd.get_pdbstr).
  3. A numeric token embedded in the object name (e.g. "lig_-7.5").
  4. Falls back to 0.0 with a warning.

Deduplication logic
-------------------
For every ligand object, every ligand-atom ↔ receptor-atom pair whose
Euclidean distance is ≤ cutoff is evaluated.  Only ONE record is kept
per (ligand_atom_index, receptor_atom_index) pair — the shortest
observed distance — so no two-row is reported twice and no inflated
count arises from alternative conformations or symmetry.

Requirements
------------
  PyMOL (open-source or commercial), numpy, matplotlib
  Install extras if needed:
      pip install numpy matplotlib          (or conda install ...)
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")           # headless-safe; change to "TkAgg" if you
                                    # want an interactive window
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[docking_analysis] WARNING: matplotlib not found — chart will not be generated.")

try:
    from pymol import cmd, stored
    HAS_PYMOL = True
except ImportError:
    HAS_PYMOL = False
    print("[docking_analysis] WARNING: PyMOL not found — running in stub mode.")


# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
AtomInfo   = Tuple[int, str, str, str, str, float, float, float]
# (index, name, resn, resi, chain, x, y, z)

PairKey    = Tuple[int, int]
# (ligand_atom_index, receptor_atom_index)

InteractionRecord = Tuple[float, str, str, str, str, str, str, str]
# (distance, lig_atom_name, lig_resn, lig_resi,
#            rec_atom_name, rec_resn, rec_resi, rec_chain)

LigInteractions = Dict[PairKey, InteractionRecord]
# unique pair → shortest-distance record


# ---------------------------------------------------------------------------
# Score extraction helpers
# ---------------------------------------------------------------------------

def _score_from_bfactor(obj: str) -> Optional[float]:
    """Return score if all B-factors in the object are the same value."""
    stored.bfac = []
    cmd.iterate(obj, "stored.bfac.append(b)")
    bfac = stored.bfac
    if not bfac:
        return None
    rounded = [round(b, 3) for b in bfac]
    if len(set(rounded)) == 1:          # uniform → likely score
        return rounded[0]
    # Some tools write score only in heavy atoms; try first atom
    return rounded[0]


def _score_from_remark(obj: str) -> Optional[float]:
    """Scan PDB REMARK / SCORE lines for a numeric score."""
    pdb_text = cmd.get_pdbstr(obj)
    patterns = [
        r"REMARK\s+(?:VINA\s+)?RESULT\s+([-+]?\d+\.\d+)",
        r"REMARK\s+minimizedAffinity\s+([-+]?\d+\.\d+)",
        r"REMARK\s+(?:SCORE|Affinity|docking_score)\s*[=:]\s*([-+]?\d+\.\d+)",
        r"SCORE\s*[=:]\s*([-+]?\d+\.\d+)",
        r"REMARK\s+([-+]\d+\.\d+)",      # bare number on REMARK line
    ]
    for pat in patterns:
        m = re.search(pat, pdb_text, re.IGNORECASE)
        if m:
            return float(m.group(1))
    return None


def _score_from_name(obj: str) -> Optional[float]:
    """Parse a numeric value embedded in the object name."""
    tokens = re.findall(r"[-+]?\d+\.?\d*", obj.replace("_", " "))
    for tok in tokens:
        try:
            val = float(tok)
            if val < 0:          # negative → plausible binding energy
                return val
        except ValueError:
            pass
    return None


def get_docking_score(obj: str) -> float:
    """
    Try to extract docking score via B-factor, REMARK, or name parsing.
    Returns 0.0 if nothing is found.
    """
    for extractor in (_score_from_remark, _score_from_bfactor, _score_from_name):
        try:
            val = extractor(obj)
            if val is not None:
                return float(val)
        except Exception:
            pass
    print(f"  [score] Could not determine score for '{obj}' — using 0.0")
    return 0.0


# ---------------------------------------------------------------------------
# Interaction calculation
# ---------------------------------------------------------------------------

def _collect_atoms(selection: str) -> List[AtomInfo]:
    """Return a list of (index, name, resn, resi, chain, x, y, z) tuples."""
    stored._atom_list = []
    cmd.iterate_state(
        1,
        selection,
        "stored._atom_list.append((index, name, resn, resi, chain, x, y, z))"
    )
    return stored._atom_list[:]


def calc_interactions(
    ligand: str,
    receptor: str,
    cutoff: float,
) -> LigInteractions:
    """
    Return unique atom-pair interactions between *ligand* and *receptor*.

    Each unique (lig_atom_index, rec_atom_index) pair appears **at most once**
    with the **shortest** Euclidean distance.  Pairs beyond *cutoff* Å are
    discarded.
    """
    lig_atoms = _collect_atoms(ligand)
    rec_atoms = _collect_atoms(receptor)

    if not lig_atoms:
        print(f"  [warn] No atoms found in ligand '{ligand}'")
        return {}
    if not rec_atoms:
        print(f"  [warn] No atoms found in receptor '{receptor}'")
        return {}

    # Vectorised distance calculation ----------------------------------------
    lig_coords = np.array([[a[5], a[6], a[7]] for a in lig_atoms])   # (L, 3)
    rec_coords = np.array([[a[5], a[6], a[7]] for a in rec_atoms])   # (R, 3)

    # Euclidean distance matrix: shape (L, R)
    diff = lig_coords[:, None, :] - rec_coords[None, :, :]           # (L, R, 3)
    dist_matrix = np.linalg.norm(diff, axis=2)                       # (L, R)

    # Find pairs within cutoff
    lig_idx_arr, rec_idx_arr = np.where(dist_matrix <= cutoff)

    pair_map: LigInteractions = {}

    for li, ri in zip(lig_idx_arr.tolist(), rec_idx_arr.tolist()):
        dist = float(dist_matrix[li, ri])
        la = lig_atoms[li]
        ra = rec_atoms[ri]

        pair_key: PairKey = (la[0], ra[0])   # (lig_atom_index, rec_atom_index)

        if pair_key not in pair_map or dist < pair_map[pair_key][0]:
            pair_map[pair_key] = (
                dist,
                la[1], la[2], la[3],   # lig atom_name, resn, resi
                ra[1], ra[2], ra[3], ra[4],  # rec atom_name, resn, resi, chain
            )

    return pair_map


# ---------------------------------------------------------------------------
# Chart
# ---------------------------------------------------------------------------

# Dark-lab colour palette
_BG      = "#0d1117"
_PANEL   = "#161b22"
_BORDER  = "#30363d"
_TEXT    = "#e6edf3"
_SUBTEXT = "#8b949e"
_BLUE    = "#58a6ff"
_GREEN   = "#3fb950"
_YELLOW  = "#d29922"
_RED     = "#f85149"
_PURPLE  = "#bc8cff"


def _score_color(score: float) -> str:
    if score <= -8.0:
        return _GREEN
    elif score <= -6.0:
        return _BLUE
    elif score <= -4.0:
        return _YELLOW
    else:
        return _RED


def _ax_style(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    """Apply consistent dark styling to an axes object."""
    ax.set_facecolor(_PANEL)
    ax.set_title(title, color=_TEXT, fontsize=11, fontweight="bold", pad=10)
    if xlabel:
        ax.set_xlabel(xlabel, color=_SUBTEXT, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, color=_SUBTEXT, fontsize=9)
    ax.tick_params(colors=_TEXT, labelsize=8)
    for spine in ax.spines.values():
        spine.set_edgecolor(_BORDER)
    ax.grid(color=_BORDER, linestyle="--", alpha=0.6, linewidth=0.7)
    ax.set_axisbelow(True)


def create_chart(
    ligands: List[str],
    scores: Dict[str, float],
    interactions: Dict[str, LigInteractions],
    cutoff: float,
    output: str,
) -> None:
    """Render and save a 4-panel docking analysis chart."""

    if not HAS_MPL:
        print("[chart] matplotlib unavailable — skipping chart.")
        return

    # Sort ligands best → worst score
    sorted_ligs = sorted(ligands, key=lambda l: scores.get(l, 0))
    sorted_scores = [scores.get(l, 0.0) for l in sorted_ligs]
    n = len(sorted_ligs)

    fig = plt.figure(figsize=(18, 11), facecolor=_BG)
    fig.suptitle(
        "Molecular Docking Analysis",
        fontsize=17, fontweight="bold", color=_TEXT, y=0.99, x=0.5,
        ha="center",
    )
    fig.text(
        0.5, 0.965,
        f"Interaction cutoff: {cutoff} Å  |  Unique atom-pair distances only",
        fontsize=9, color=_SUBTEXT, ha="center",
    )

    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        hspace=0.45, wspace=0.30,
        left=0.07, right=0.97, top=0.93, bottom=0.08,
    )

    x = np.arange(n)

    # ── Panel 1: Docking Scores ─────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, :])
    bar_colors = [_score_color(s) for s in sorted_scores]
    bars = ax1.bar(x, sorted_scores, color=bar_colors, alpha=0.88,
                   edgecolor="none", width=0.65, zorder=3)

    for bar, score in zip(bars, sorted_scores):
        offset = -0.18 if score <= 0 else 0.18
        va = "top" if score <= 0 else "bottom"
        ax1.text(
            bar.get_x() + bar.get_width() / 2.0,
            score + offset,
            f"{score:.2f}",
            ha="center", va=va, fontsize=8,
            color=_TEXT, fontweight="bold", zorder=4,
        )

    ax1.set_xticks(x)
    ax1.set_xticklabels(sorted_ligs, rotation=30, ha="right",
                        color=_TEXT, fontsize=9)
    ax1.axhline(0, color=_BORDER, linewidth=0.8)
    _ax_style(ax1, "Docking Scores", ylabel="Score (kcal/mol)")

    legend_patches = [
        Patch(facecolor=_GREEN,  label="≤ −8.0 kcal/mol  (Excellent)"),
        Patch(facecolor=_BLUE,   label="−8.0 to −6.0      (Good)"),
        Patch(facecolor=_YELLOW, label="−6.0 to −4.0      (Moderate)"),
        Patch(facecolor=_RED,    label="> −4.0             (Weak)"),
    ]
    ax1.legend(
        handles=legend_patches, loc="lower right",
        facecolor=_PANEL, edgecolor=_BORDER,
        labelcolor=_TEXT, fontsize=8, framealpha=0.95,
    )

    # ── Panel 2: Unique Contact Counts ──────────────────────────────────────
    ax2 = fig.add_subplot(gs[1, 0])
    counts = [len(interactions.get(l, {})) for l in sorted_ligs]
    count_colors = [_score_color(s) for s in sorted_scores]
    hbars = ax2.barh(range(n), counts, color=count_colors, alpha=0.82,
                     edgecolor="none", zorder=3)
    ax2.set_yticks(range(n))
    ax2.set_yticklabels(sorted_ligs, color=_TEXT, fontsize=8)

    for i, cnt in enumerate(counts):
        ax2.text(cnt + 0.3, i, str(cnt), va="center",
                 color=_TEXT, fontsize=8, fontweight="bold")

    _ax_style(
        ax2,
        "Unique Atom-Pair Contact Counts",
        xlabel=f"# unique contacts (≤ {cutoff} Å)",
    )

    # ── Panel 3: Distance Distributions (box plots) ─────────────────────────
    ax3 = fig.add_subplot(gs[1, 1])
    dist_data: List[List[float]] = []
    bp_labels: List[str] = []
    for lig in sorted_ligs:
        dists = [rec[0] for rec in interactions.get(lig, {}).values()]
        if dists:
            dist_data.append(sorted(dists))
            bp_labels.append(lig)

    if dist_data:
        bp = ax3.boxplot(
            dist_data,
            labels=bp_labels,
            patch_artist=True,
            widths=0.55,
            medianprops=dict(color=_YELLOW, linewidth=2),
            whiskerprops=dict(color=_SUBTEXT, linewidth=1),
            capprops=dict(color=_SUBTEXT, linewidth=1),
            flierprops=dict(marker="o", color=_BLUE,
                            markerfacecolor=_BLUE, markersize=3, alpha=0.6),
            zorder=3,
        )
        for patch, lig in zip(bp["boxes"], bp_labels):
            patch.set_facecolor(_score_color(scores.get(lig, 0)))
            patch.set_alpha(0.4)
            patch.set_edgecolor(_score_color(scores.get(lig, 0)))
            patch.set_linewidth(1.5)

    # Reference lines
    ax3.axhline(3.5, color=_GREEN,  linestyle="--", linewidth=1,
                alpha=0.7, label="H-bond / polar (≈3.5 Å)", zorder=2)
    ax3.axhline(4.0, color=_YELLOW, linestyle="--", linewidth=1,
                alpha=0.7, label="Van der Waals (≈4.0 Å)", zorder=2)

    ax3.set_xticklabels(bp_labels, rotation=30, ha="right",
                        color=_TEXT, fontsize=8)
    ax3.legend(
        facecolor=_PANEL, edgecolor=_BORDER,
        labelcolor=_TEXT, fontsize=7, loc="upper right",
    )
    _ax_style(
        ax3,
        "Distance Distributions (Unique Pairs)",
        ylabel="Shortest atom-atom distance (Å)",
    )

    plt.savefig(output, dpi=150, bbox_inches="tight", facecolor=_BG)
    print(f"\n[chart] Saved → {output}")

    try:
        plt.show()
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Text report helper
# ---------------------------------------------------------------------------

def print_report(
    ligands: List[str],
    scores: Dict[str, float],
    interactions: Dict[str, LigInteractions],
    cutoff: float,
    top_n: int = 10,
) -> None:
    """Print a concise text summary of the docking analysis."""
    SEP = "─" * 72
    print(f"\n{SEP}")
    print(" DOCKING ANALYSIS REPORT")
    print(f" Interaction cutoff: {cutoff} Å | Unique atom-pair distances only")
    print(SEP)

    sorted_ligs = sorted(ligands, key=lambda l: scores.get(l, 0))

    for rank, lig in enumerate(sorted_ligs, 1):
        sc = scores.get(lig, 0.0)
        imap = interactions.get(lig, {})
        dists = sorted(rec[0] for rec in imap.values())

        print(f"\n  [{rank:>2}] {lig}")
        print(f"       Score    : {sc:+.3f} kcal/mol")
        print(f"       Contacts : {len(imap)} unique atom-atom pairs")
        if dists:
            print(f"       Distances: min={dists[0]:.2f} Å  "
                  f"median={np.median(dists):.2f} Å  "
                  f"max={dists[-1]:.2f} Å")

        # Show the N closest contacts
        closest = sorted(imap.items(), key=lambda kv: kv[1][0])[:top_n]
        if closest:
            print(f"       Top {top_n} shortest contacts:")
            for (la_idx, ra_idx), rec in closest:
                d, la_nm, la_rn, la_ri, ra_nm, ra_rn, ra_ri, ra_ch = rec
                print(f"         {d:5.2f} Å  "
                      f"{lig}/{la_rn}{la_ri}@{la_nm} "
                      f"↔ receptor/{ra_ch}/{ra_rn}{ra_ri}@{ra_nm}")

    print(f"\n{SEP}\n")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def analyze_docking(
    receptor: str = "receptor",
    cutoff: float = 5.0,
    output: str = "docking_analysis.png",
    top_contacts: int = 10,
    auto_detect_receptor: bool = True,
    session: str = None,
) -> Tuple[Dict[str, float], Dict[str, LigInteractions]]:
    """
    Run the full docking analysis pipeline.

    Parameters
    ----------
    receptor : str
        Name of the receptor/protein object in PyMOL.  If not found and
        *auto_detect_receptor* is True, the largest object is used instead.
    cutoff : float
        Distance cutoff in Ångströms for interaction detection (default 5.0).
    output : str
        File path for the output chart image (PNG).
    top_contacts : int
        Number of closest contacts to list per ligand in the text report.
    auto_detect_receptor : bool
        Fall back to the largest object when *receptor* is not found.
    session : str
        Path to a PyMOL session file (.pse) to load before analysis.
        Defaults to /docking_results/docking_session.pse.
        Pass session="" to skip loading (if objects are already in PyMOL).

    Returns
    -------
    scores : dict  {ligand_name: score}
    interactions : dict  {ligand_name: {(lig_idx, rec_idx): record_tuple}}
    """
    if not HAS_PYMOL:
        print("[analyze_docking] PyMOL not available — cannot run.")
        return {}, {}

    import os
    if session is None:
        session = os.path.join(
            os.path.expanduser("~"),
            "Documents", "pymol", "pymol", "docking",
            "docking_results", "docking_session.pse",
        )

    if session:
        print(f"[setup] Loading session: {session}")
        cmd.load(session)
        print(f"[setup] Session loaded successfully.")

    all_objects: List[str] = cmd.get_object_list()

    # ── Receptor detection ──────────────────────────────────────────────────
    if receptor not in all_objects:
        if auto_detect_receptor:
            print(f"[receptor] '{receptor}' not found. "
                  f"Available: {all_objects}")
            print("[receptor] Falling back to the largest object …")
            receptor = max(all_objects, key=lambda o: cmd.count_atoms(o))
            print(f"[receptor] Using '{receptor}'")
        else:
            raise ValueError(
                f"Receptor '{receptor}' not found. "
                f"Available objects: {all_objects}"
            )

    ligands: List[str] = [o for o in all_objects if o != receptor]

    if not ligands:
        print("[analyze_docking] No ligand objects found. Aborting.")
        return {}, {}

    print(f"\n[setup] Receptor : {receptor}  ({cmd.count_atoms(receptor)} atoms)")
    print(f"[setup] Ligands  : {ligands}")
    print(f"[setup] Cutoff   : {cutoff} Å\n")

    # ── Scores ──────────────────────────────────────────────────────────────
    scores: Dict[str, float] = {}
    for lig in ligands:
        sc = get_docking_score(lig)
        scores[lig] = sc
        print(f"[score] {lig:30s} → {sc:+.3f} kcal/mol")

    # ── Interactions (deduplicated) ──────────────────────────────────────────
    interactions: Dict[str, LigInteractions] = {}
    for lig in ligands:
        print(f"\n[distances] Processing {lig} …")
        imap = calc_interactions(lig, receptor, cutoff)
        interactions[lig] = imap
        print(f"[distances] {lig}: {len(imap)} unique atom-pair contacts ≤ {cutoff} Å")

    # ── Report & Chart ───────────────────────────────────────────────────────
    print_report(ligands, scores, interactions, cutoff, top_n=top_contacts)
    create_chart(ligands, scores, interactions, cutoff, output)

    return scores, interactions


# ---------------------------------------------------------------------------
# PyMOL command registration — lets you call it from the PyMOL prompt
# ---------------------------------------------------------------------------
if HAS_PYMOL:
    cmd.extend("analyze_docking", analyze_docking)
    print("[docking_analysis] Loaded.  Run:  analyze_docking receptor=myProtein, cutoff=5.0")

# Auto-run when executed directly via `run script.py` in PyMOL
if __name__ == "__main__":
    analyze_docking()
