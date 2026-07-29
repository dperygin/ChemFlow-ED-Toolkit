# ChemFlow-ED: Suggested Student Assignment

**Course context:** Senior-level/advanced undergraduate or graduate medicinal chemistry
**Duration:** 4 weeks
**Format:** Independent project, individual therapeutic target, common computational methodology

## Assignment Overview

You will carry out a structure-based drug discovery investigation beginning
with only an experimentally determined protein structure — no known ligand or
inhibitor will be provided. Using the ChemFlow-ED workflow, you will identify
a candidate binding pocket, build a receptor-derived pharmacophore, run a
virtual screen, dock and analyze candidate compounds, and prioritize a
defensible shortlist for future experimental evaluation. Your final deliverable
is a publication-style scientific manuscript.

## Step 1 — Target Selection and Literature Review (Week 1)

Select a therapeutically relevant protein with a structure available in the
Protein Data Bank. Preference is given to proteins associated with inherited
metabolic disorders, neurodegenerative disease, or other areas in which a
small-molecule stabilizer or modulator would be of therapeutic interest, and
which currently lack a well-characterized small-molecule ligand.

Write a 1-2 page literature summary covering:
- The protein's biological function and structure
- The disease or condition it is associated with
- Why a computational, ligand-independent approach is appropriate for this
  target
- Your proposed strategy (e.g., stabilization vs. inhibition)

## Step 2 — Pocket Identification and Pharmacophore Generation (Week 2)

Using DoGSiteScorer and Pharmit Forge (see `docs/WORKFLOW_GUIDE.md`), identify
and justify a candidate binding pocket, then generate and refine a
receptor-derived pharmacophore. Document:
- The candidate pockets you considered and why you selected the one you did
- The initial and refined pharmacophore features, and your rationale for any
  changes

## Step 3 — Virtual Screening and Docking (Week 3)

Run a pharmacophore-based virtual screen in Pharmit, then dock the resulting
candidates using the ChemFlow-ED automation scripts. Document:
- The size and composition of your initial and refined hit lists
- Your top-ranked candidates by predicted binding affinity

## Step 4 — Interaction Analysis and Prioritization (Week 3–4)

Run the interaction-analysis scripts on your top candidates. Compare
interaction patterns across compounds with similar docking scores. Prioritize
a shortlist (typically 3–5 compounds) using multiple lines of evidence, and
explicitly justify your reasoning.

## Step 5 — Manuscript and Presentation (Week 4)

Prepare a publication-style manuscript (ACS format) including:
- Introduction and biological background
- Computational methodology
- Results (pocket selection, pharmacophore, screening, docking, interaction
  analysis)
- Discussion, including limitations of the computational predictions
- A brief oral presentation (10–15 minutes) summarizing your findings

## Important Framing

Your goal is **not** to identify a single "correct" compound. Computational
docking predicts relative binding affinity and plausible orientation; it does
not establish biological activity. Use language such as "candidate,"
"predicted," "hypothesized," and "prioritized" throughout your report, and
discuss what experimental follow-up would be needed to test your top
candidates.
