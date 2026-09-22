# receptor-configuration — analysis and figure code

Code behind the figures of *"Fixed installation, inducible gate: a human B-cell receptor
configuration separates lipid A sensing from inflammatory licensing"* (H. Nagashima).

Every panel is drawn by a script that reads data and writes both the figure
(`png`/`pdf`/`svg`) and the source table behind it (`fig*_source.csv`). No numbers are typed
into the plotting code. `numpy` seed is fixed (`SEED = 0`, `common.py`). All input data are
public and were treated as read-only; no new data were generated for this study.

## Quick start

```bash
export DATA_ROOT=/path/where/you/put/the/inputs   # default: /Volumes/M4_SSD/projects
python3 fig1_module_placement.py
python3 fig2_fixed_installation.py
python3 fig3_gate_not_receiver.py
python3 fig4_organ_placement.py
```

Figures and source tables are written next to the scripts. A missing or empty input stops the
script with a message naming the file it wants (`common.need`).

Requires Python 3 with `numpy`, `pandas`, `scipy`, `matplotlib`, and — for Fig. 4 only —
`scanpy` (to read `b-cells.h5ad` in backed mode).

## Data manifest

Paths are relative to `DATA_ROOT`. Inputs fall into three classes; please read class C before
assuming the repository is self-contained.

### A. Public raw files — download these yourself

| Path under `DATA_ROOT` | Source | Used by |
|---|---|---|
| `tlr10_crosscelltype/data/GSE107011_Processed_data_TPM.txt.gz` | GEO GSE107011 (Monaco, 29 sorted immune subsets) | Fig. 1a |
| `bcell_lipidA_whitespace/GSE327707/data/GSE327707_raw_counts.csv.gz`, `.../coldata.csv` | GEO GSE327707 (sorted PBMC populations across stimuli) | Fig. 1b, Fig. 2a |
| `spencer_galt/geo/GSE327417/extracted/GSM9656505_318B-2_features.tsv.gz` | GEO GSE327417 (used only as an Ensembl→symbol table) | Fig. 1a, Fig. 2a |
| `tica_pan_immune_bcell/data/b-cells.h5ad` | tissueimmunecellatlas.org (cross-tissue immune atlas, B cells) | Fig. 4a |
| `md1_ladder/rna_single_cell_type.tsv` | Human Protein Atlas | Fig. 4b |
| `tlr10_uc_scrna/data/smillie_imm_export/` (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`, cell metadata) | Broad Single Cell Portal SCP259 (Smillie 2019, UC lamina propria immune cells) | `prep_smillie_lineage.py` → Fig. 1c |
| SCP3755 (UC single-cell atlas) Seurat object and its cell metadata CSV | Broad Single Cell Portal SCP3755 | `pb_counts.csv`, `pb_meta.csv`, `scp3755_metadata.csv` below |

### B. Small derived tables — included in this repository

These are intermediates, not raw data. They are small enough to ship, so Fig. 1c and the
patient-level panels run without rebuilding them.

| File | What it is |
|---|---|
| `fig1c_smillie.csv` | lineage means of the 7 configuration genes in SCP259, written by `prep_smillie_lineage.py` |
| `fig*_source.csv`, `fig3a_delta.csv`, `fig3a_p.csv`, `fig3a_nullp.csv`, `fig4a_subset.csv` | the source table behind each published panel |
| `fig4a_cells.csv.gz` | per-cell table from `b-cells.h5ad`, three columns only (organ, curated cell type, TLR10 value) — no barcodes, no donor identifiers, no clinical metadata |

### C. Large derived inputs — NOT in this repository

Fig. 2b and Fig. 3 read pseudobulk objects built from the SCP3755 atlas. They exceed sensible
repository size, so they are not shipped, and **a clone of this repository alone will not
reproduce those two figures.** Rebuild them from the SCP3755 download with
`build_pseudobulk.R` in this repository, or ask the author for the intermediates.

| Path under `DATA_ROOT` | Size | Used by |
|---|---|---|
| `scp3755_uc_atlas/pb_counts.csv` | 84 MB | Fig. 3 |
| `scp3755_uc_atlas/pb_meta.csv` | 68 KB | Fig. 2b, Fig. 3 |
| `scp3755_uc_atlas/E_logcpm.pkl` | 144 KB | Fig. 2b, Fig. 3 |
| `scp3755_uc_atlas/scp3755_metadata.csv` | 103 MB | Fig. 2b |
| `gbp2_cytosolic_arm/patient_lineage.csv` | 28 KB | Fig. 1c, Fig. 2b, Fig. 3 |

Note for anyone reproducing this on the development machine: `scp3755_metadata.csv` and
`scp3755.RDS` are symbolic links to an external disk. If that disk is not mounted they read as
empty and `common.need` will say so rather than failing silently.

## What each script does

- `fig1_module_placement.py` — (a) MD-1(LY86):MD-2(LY96) ratio across 29 sorted subsets;
  (b) baseline CPM in cleanly sorted populations: B cells carry TLR4 but not MD-2 or CD14;
  (c) the three-tier layout in colonic single cells.
- `fig2_fixed_installation.py` — (a) TLR10 across stimuli in sorted cells; (b) within-patient
  inflamed/uninflamed pairs in SCP3755. TLR10 does not move in either.
- `fig3_gate_not_receiver.py` — the interferon-inducible arm rises lineage-selectively while
  CASP4 does not; surface reader and cytosolic gate are uncoupled at single-cell resolution.
  Expression-matched null calibration included (`fig3a_nullp.csv`).
- `fig4_organ_placement.py` — configuration by organ and by cell type (cross-tissue atlas, HPA).
- `build_pseudobulk.R` — rebuilds the SCP3755 pseudobulk inputs of class C (`pb_counts.csv`, `pb_meta.csv`) from the downloaded Seurat object; set `SCP3755_DIR` or `DATA_ROOT`.
- `prep_smillie_lineage.py` — one pass over the SCP259 matrix to cache lineage means.
- `common.py` — paths (`DATA_ROOT`, `FIGDIR`, `need`), colour-blind-safe palette, fixed seed.

## Licence

Code: MIT. The input datasets keep the licences of their original repositories.
