# Manifest and execution

## CSV columns

Use a UTF-8 CSV. Paths may be absolute or relative to the CSV directory. A row is one card; row order determines pairings and page order. The tool reads model 1 (the first record) of each PDB/mmCIF and leaves the source files intact.

| Column | Value |
| --- | --- |
| `id`, `structure` | Required. Unique card ID and path to the existing `.pdb`, `.cif`, or `.mmcif`. |
| `label`, `model`, `seed` | Optional display values, e.g. `CDL full`, `AF3`, `42`. Keep unknown seeds blank. Use short Latin labels in the PDF; UTF-8 text is preserved in the CSV. |
| `protein_chains` | Source chain IDs separated by spaces. Blank selects amino acid residues with Cα from all chains. Protein chains stay separate in the staging PDB; contact labels retain source numbering. |
| `ligand_code`, `ligand_chain`, `ligand_resid` | Source residue name, chain ID and author residue ID (include an insertion code if present). All provided selectors must match exactly. Auto selection works only when one non-water, non-protein candidate matches. Each card supports one ligand residue. |
| `head_atoms` | Space separated names from the selected ligand residue, e.g. `P O1 O2 C1`. Blank selects the entire ligand for the pocket. For a full phospholipid, use its actual head atom names (inspect the CCD and atom mapping). The remainder is counted separately and appears only in both overviews. |
| `ccd` | Optional component mmCIF containing `_chem_comp_bond` records with atom names matching the coordinates. Embedded mmCIF component bonds or PDB `CONECT` records are used when `ccd` is blank. Rendered sticks use the graph without encoding bond order. |
| `infer_bonds` | Set `true` only to allow covalent-distance inference when embedded records are unavailable. An inferred connected graph can still have wrong edges. Visually inspect it or supply a CCD for delivery. |
| `align_group` | Optional comparison group. All members must share the same ordered protein chain IDs, author residue IDs and residue names. Uses Cα alignment to the first row in the group; its binding face defines both views and the group shares their scale and crop. Use the desired reference row first. The preview computes framing within its selected subset. |
| `iptm`, `ipsae`, `mean_pae`, `protein_plddt`, `ligand_plddt` | Optional precomputed metrics in 0–1, 0–1, Å, 0–100 and 0–100, respectively. Values override extras derived from the optional AF3 files. Use consistent confidence definitions for comparison. For Boltz, Protenix, OpenDDE or other sources provide the precomputed values and their method in the accompanying report. |
| `af3_confidence`, `af3_summary` | Optional AF3 full confidence and summary JSONs for the exact coordinate sample. See below. |

Example (replace with real source paths and CCD atom names):

```csv
id,structure,label,model,seed,protein_chains,ligand_code,ligand_chain,head_atoms,ccd,align_group,af3_confidence,af3_summary
X1-head-s42,inputs/X1_head.cif,Head,AF3,42,A,HDG,B,,inputs/HDG.cif,X1,inputs/X1_head_confidences.json,inputs/X1_head_summary_confidences.json
X1-full-s42,inputs/X1_full.cif,Full ligand,AF3,42,A,LIP,B,P O1 O2 C1,inputs/LIP.cif,X1,inputs/X1_full_confidences.json,inputs/X1_full_summary_confidences.json
```

To migrate an existing table, join metadata by design, model, ligand form and seed (or a unique structure hash). For example map `wet_structure` to `structure`, attach the accession or form as `label`, and map `ipSAE_max_pae10`, `cross_pae_mean_A`, `protein_plddt_mean`, and `ligand_plddt_mean` to `ipsae`, `mean_pae`, `protein_plddt`, and `ligand_plddt`. Set `align_group` to the design ID only if the ordered residue IDs and names correspond. The script does not depend on those legacy column names or paths.

### Reading AF3 confidence

The full confidence file must map the selected protein chains to residue tokens and the ligand to heavy atom tokens on a separate chain. Chain IDs and counts are checked against the coordinate file; the caller must ensure both files are from the same sample. `mean_pae` is the mean of both cross-chain blocks, outside of any PAE cutoff. `ipsae` follows the [Dunbrack d0res expression](https://github.com/DunbrackLab/IPSAE) with PAE < 10 Å, taking the maximum of both directions and using ligand heavy atoms as tokens. This ligand adaptation is not calibrated as a measure of ligand specificity. pLDDT is the atom mean from AF3, on 0–100. Summary iPTM is the overall value; supply a pair-specific value explicitly if needed for a complex with multiple partners. Crystallographic B factors are never used as confidence.

## Execution and output

Use Python 3.10+ with `gemmi`, `numpy`, `scipy`, `Pillow` and `PyMuPDF` (`import pymupdf`). VMD must support text mode, NewCartoon, QuickSurf, and Tachyon scene export; use the matching external CPU Tachyon executable. Common installs expose `vmd` on `PATH` and put `tachyon_*` in `$VMDDIR` or `<prefix>/lib/vmd/`. Other layouts can set `VMD_ATLAS_VMD`, `VMD_ATLAS_TACHYON`, `--vmd`, or `--tachyon`. On the original workstation, the `rfdf4` Conda environment provides Python packages, and `cdl_vmd_review` provides `bin/vmd` and `lib/vmd/tachyon_LINUXAMD64`. Those environments are optional; the code does not require them.

`--preview N` renders the first N rows, `--rows-per-page 1–8` controls density, and `--contact` (default `4.0 Å`) and `--pocket-radius` (default `8.5 Å`) define inclusion. `--workers` controls the number of active cards per run. One worker creates three temporary Tachyon scenes, traces them on the CPU, and keeps only the finished images; start low on a shared machine or with very large complexes. Each output directory must be new. Use the named card and `logs/` to diagnose VMD errors or blank/clipped renders.

The PDF includes page bookmarks. `structure_index.csv` maps each card to its source, source SHA256, model/seed, PDF page, pocket residues, score provenance, and topology provenance. `heavy_atom_contacts.csv` lists all pairs within the contact cutoff with original atom names and distances, including atoms outside the selected head. The picture numbers are in `badge`; at most three head pairs are shown, from different protein residues. Priorities are basic–phosphate, polar, carbon, then other proximity; pairs shorter than 2.2 Å are flagged as overlaps. Labels do not establish hydrogen bonds or binding preferences. Ribbon ends in the pocket view are selection boundaries, with secondary structure assigned on the intact staged protein. `views/` holds unannotated PNGs (the contact dashes are rendered; numbers and labels are added in the PDF).
