# Release Instructions

These steps produce a versioned, citable release of ChemFlow-ED that matches
the version and DOI referenced in the manuscript's Data and Software
Availability section (currently: **v1.0.2**, DOI **10.5281/zenodo.21302160**).
If you change the version number or generate a new DOI, update the manuscript,
`CITATION.cff`, and `README.md` to match — these must all agree before
submission.

## 1. Pre-release checklist

- [ ] All scripts run cleanly against `example_data/` (spot-check at least one
      full pass of `docs/WORKFLOW_GUIDE.md`, Steps 1–9).
- [ ] `docking_affinities_summary.csv` files under `example_data/` contain
      only runs for the intended target/pocket (no mixed-in results from other
      projects or targets — check `receptor` and `subfolder` columns).
- [ ] `requirements.txt` versions match what you actually tested with.
- [ ] No personal file paths, usernames, or informal/scratch notes remain in
      files intended for public release (check `scripts/*.txt` and any
      `workflow.txt`/notes files under `example_data/` — see the consistency
      review for specific files to check).
- [ ] Terminology is consistent across `docs/WORKFLOW_GUIDE.md`, the
      manuscript, `supporting_information/Project_Timeline.md`, and any other
      document that names a specific tool (e.g., confirm "PharmacoForge" vs.
      "Pharmit Forge" is used correctly and consistently everywhere it
      appears — these may be two different, real tools; don't conflate them).
- [ ] `CITATION.cff` version, date, and DOI match what you intend to publish.
- [ ] `LICENSE`, `README.md`, and the manuscript's author name/affiliation/
      ORCID all agree.

## 2. Tag and push the release commit

```bash
git add -A
git commit -m "Prepare v1.0.2 release"
git tag -a v1.0.2 -m "ChemFlow-ED v1.0.2"
git push origin main --tags
```

## 3. Create the GitHub release

1. On GitHub, go to your repository → **Releases** → **Draft a new release**.
2. Choose the tag you just pushed (`v1.0.2`).
3. Title: `ChemFlow-ED v1.0.2`.
4. Description: briefly summarize what changed since the last version (or,
   for a first release, summarize the contents — workflow guide, automation
   scripts, example dataset, supporting instructional materials).
5. Click **Publish release**.

## 4. Archive the release on Zenodo

**If this is your first release of this repository:**

1. Log in to [Zenodo](https://zenodo.org) with your GitHub account.
2. Go to **GitHub** under your Zenodo account settings and flip the toggle
   "on" for the `ChemFlow-ED` repository *before* publishing the GitHub
   release above — Zenodo only archives releases published after the
   repository is enabled.
3. Publishing the GitHub release (Step 3) will trigger Zenodo to
   automatically create an archive and mint a DOI.
4. On Zenodo, open the newly created record and fill in/confirm: title,
   author name and ORCID, license (MIT), and a description matching the
   README.
5. Note the DOI Zenodo assigns. If it does not match
   `10.5281/zenodo.21302160`, update the manuscript, `README.md`, and
   `CITATION.cff` to the DOI Zenodo actually issued — **do not leave a
   placeholder or guessed DOI in the manuscript**.

**If you already have a Zenodo record for a previous version:**

1. Ensure GitHub-Zenodo sync is still enabled for the repository.
2. Publish the GitHub release (Step 3) — Zenodo will automatically create a
   new version linked to the same "concept DOI," with its own
   version-specific DOI.
3. Confirm the new version-specific DOI and update the manuscript and
   `CITATION.cff` accordingly.

## 5. Verify before submitting the manuscript

- [ ] Visit `https://doi.org/<your-DOI>` in a private/incognito window and
      confirm it resolves to the correct Zenodo record.
- [ ] Confirm the Zenodo record's license, author, and version match the
      manuscript exactly.
- [ ] Confirm the GitHub repository is public (not private) and the license
      badge/file is visible.
- [ ] Add a full reference-list entry for the dataset/software deposit in the
      manuscript's References section (ACS's Research Data Policy requires
      citing publicly available underlying data in-text with a full citation
      in the reference list, not just a mention in prose). Suggested format,
      per the ACS Style Quick Guide:

  > Perygin, D. H. ChemFlow-ED (version v1.0.2). *Zenodo*. American Chemical
  > Society; Jacksonville State University, 2026. DOI:
  > 10.5281/zenodo.21302160.

  and cite it with a superscript reference number at its mention in the Data
  and Software Availability section, the same way every other source in the
  manuscript is cited.

## 6. After publication

- If the manuscript is later revised and the DOI or version changes, repeat
  Steps 2–5 and update the manuscript before resubmission.
- Do not delete or unpublish a Zenodo record once it has been cited in a
  published article — mint a new version instead if changes are needed.
