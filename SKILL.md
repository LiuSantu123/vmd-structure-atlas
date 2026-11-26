---
name: vmd-structure-atlas
description: Render existing protein–small molecule PDB or mmCIF structures with VMD and Tachyon as three-view cards and multi-page review PDFs. Use for binding-site or pocket figures, or for visual comparisons across ligands, seeds, or structure predictors.
---

# VMD Structure Atlas

Turn existing coordinates into a review PDF. Each card shows the complete complex from the binding face and 82° to the side, then a pocket view with the local protein surface, the cartoon of the same residues, selected residue sticks, and numbered heavy atom distances. The two overviews use protein cartoon plus a 75% transparent surface and show the full ligand as sticks. The pocket defaults to an 8.5 Å shell around the ligand atoms selected for display. The A3 layout has two cards (six views) per row and up to eight rows per page.

1. Make a CSV with one row per structure and unique `id` and `structure` columns. Read [the manifest reference](references/manifest.md) for selectors, head/tail handling, confidence sources, and a worked CSV example. For full lipids, specify the head atoms explicitly so that the two overviews keep the tails and the pocket shows the head. Confirm the ligand has an embedded CCD graph or PDB `CONECT` records, or provide a matching CCD.
2. Pair comparable variants in CSV row order. Set `align_group` on rows that should use the first structure's Cα alignment, viewing axes and common scale. Leave it blank for independent framing or when residue correspondence has not been established. Use confidence data for the actual coordinates shown. Missing values render as `—`.
3. Run the bundled script with an environment that imports `gemmi`, `numpy`, `scipy`, `PIL`, and `pymupdf`. Use VMD in text mode and the CPU Tachyon executable distributed with it. Start with a small preview when choosing a style for a new dataset, then render the full manifest to a separate directory.

```bash
python "$SKILL_DIR/scripts/render_atlas.py" manifest.csv --out atlas_preview --preview 2 --rows-per-page 2
python "$SKILL_DIR/scripts/render_atlas.py" manifest.csv --out atlas_final --title 'Complex structure review' --workers 4
```

Set `SKILL_DIR` to the directory containing this `SKILL.md`. Set `--vmd` and `--tachyon` to executable paths if they are outside `PATH`. See [the manifest reference](references/manifest.md) for existing environment paths and troubleshooting. Use a new output directory for each run.

Deliver `atlas.pdf` with `structure_index.csv` and `heavy_atom_contacts.csv`; `views/` contains the individual PNGs. For a layout review, check the intended page size, centering, full ligand visibility, readable contact numbers, and agreement between those numbers and the contact table. Polar or basic–phosphate labels describe geometric proximity. Treat possible overlaps and inferred ligand bonds as items to review before interpreting the contacts.
