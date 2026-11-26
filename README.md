# VMD Structure Atlas

A [Codex skill](SKILL.md) for reviewing existing protein–small molecule structures. Each card has two complete views and a close-up of the binding pocket, with a local protein surface, cartoon, selected side chains, and numbered heavy atom distances. The tool builds an A3 PDF and exports the individual PNGs, a structure index, and a contact table.

## Install for Codex

```bash
skill_dir="${CODEX_HOME:-$HOME/.codex}/skills"
mkdir -p "$skill_dir"
git clone https://github.com/LiuSantu123/vmd-structure-atlas.git "$skill_dir/vmd-structure-atlas"
```

For a private repository, authenticate Git before cloning. In Codex, ask:

> Use `$vmd-structure-atlas` to turn these AF3 structures into a review PDF. Put headgroup and full ligand structures side by side, label the model and seed, and make a preview first.

## Run directly

You need Python 3.10+ with `gemmi`, `numpy`, `scipy`, `Pillow`, and `PyMuPDF`, plus a local installation of VMD with its CPU Tachyon executable. Prepare a [CSV manifest](references/manifest.md), then run:

```bash
python scripts/render_atlas.py manifest.csv --out atlas_preview --preview 2 --rows-per-page 2
```

Use `--vmd` and `--tachyon` if the executables are outside `PATH`. For a full run, omit `--preview` and use a new output directory. The PDF is saved as `atlas.pdf` inside that directory.

The repository contains the skill and renderer. Supply your own PDB/mmCIF files and any matching confidence files. This packaging has passed Codex skill format validation; an end-to-end render with the packaged CLI has not yet been verified.
