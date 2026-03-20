# Homopolymer analysis (Nanopore)

Scripts for **read-level QC**, **homopolymer accuracy**, **mitochondrial annotation context**, **circos layouts**, and **comparative homopolymer statistics**, using outputs from alignment / quality tools (e.g. Giraffe-style `Observed_information.txt`, `Homoploymer_summary.txt`, `homopolymer_in_reference` tables), **GFF3 / BED**, **Homopolish** site lists, and optional **NCBI taxonomy** dumps.

---

## Table of contents

- [Overview](#overview)
- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Installation](#installation)
- [Scripts](#scripts)
  - [Python](#python)
  - [R](#r)
- [Usage notes](#usage-notes)
- [Reproducibility](#reproducibility)
- [Data availability](#data-availability)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

This folder contains **2 Python** and **5 R** scripts that together support:

- **Read-level QC**: observed vs estimated per-read accuracy (Giraffe-style tables).
- **Homopolymer accuracy**: species- and library-type summaries from `Homoploymer_summary.txt`; point accuracy along the genome (e.g. T homopolymers).
- **Mitochondrial context**: homopolymer intervals overlapped with CDS / tRNA / rRNA (GFF3); circos gene + homopolymer tracks for multiple references; Homopolish sites on a gene model.
- **Error / indel panels**: mean insertion and deletion per read; homopolymer-class and position-related panels for selected datasets.
- **Taxonomy-linked summaries**: family–genus lists from NCBI `names.dmp` / `nodes.dmp`; phylum-level stacked bars from species homopolymer count spreadsheets plus NCBI lineage tables.

> **Note:** Paths use `setwd`, `~/Desktop/`, and other fixed locations. Before publishing or sharing, point these to your data layout or use project-relative paths (e.g. [`here`](https://here.r-lib.org/)).

Overall pipeline
[Pipeline.tif](https://github.com/user-attachments/files/26134399/Pipeline.tif)

---

## Repository layout

```text
code/
├── README.md
├── requirements.txt
├── extract_family_genus.py       # NCBI taxonomy → family/genus Excel
├── extract_number of homo.py     # FASTA → homopolymer intervals (≥4 bp); see filename note below
├── Figure1.R                     # QC boxplots + homopolymer accuracy scatter (multi-species)
├── Figure2.R                     # Circos: gene track + homopolymer length bands (6 genomes)
├── Figure3.R                     # Ins/del bars + homopolymer / error / accuracy panels (HK2025, fly)
├── Figure4.R                     # Feature–homopolymer bars + duplex T accuracy + Homopolish / BED
└── Figure5.R                     # Phylum (and Nematoda order) vs long-homopolymer stacked bars
```

**Filename note:** `extract_number of homo.py` contains a space. For Git and command-line use, consider renaming to e.g. `extract_homopolymer_positions.py` and updating any references.

---

## Requirements

### Python

- **Version:** 3.8+
- **Packages:** see `requirements.txt`
  - `extract_family_genus.py`: **pandas**, **openpyxl** (Excel export)
  - `extract_number of homo.py`: standard library only

### R

| Package | Typical use |
|--------|-------------|
| `ggplot2`, `dplyr`, `tidyr`, `readr`, `stringr`, `forcats` | Plotting and wrangling |
| `tidyverse` | Meta-package (Figure1) |
| `ggdist` | Distribution / half-eye plots (Figure1) |
| `patchwork` | Multi-panel layout (Figure3) |
| `scales` | Axis and percent labels |
| `circlize` | Circos (Figure2) |
| `GenomicRanges`, `rtracklayer` | GFF / interval overlap (Figure4) |
| `readxl` | Excel input (Figure5) |

Install in R:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr", "forcats",
  "tidyverse", "ggdist", "patchwork", "scales", "circlize", "readxl"
))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

---

## Installation

```bash
git clone https://github.com/<your-username>/homopolymer-analysis.git
cd homopolymer-analysis/code   # or your repo path containing these scripts
python3 -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

---

## Scripts

### Python

#### `extract_family_genus.py`

Reads NCBI taxonomy **`names.dmp`** and **`nodes.dmp`** (pipe-separated), keeps scientific names, walks the tree to attach each genus to a family, and writes **`family_genus_list.xlsx`** (family, genus list, genus count). Edit the input/output paths at the top of the script.

#### `extract_number of homo.py`

Reads a **FASTA**, finds homopolymer runs of **A/T/C/G with length ≥ 4** (regex), writes a tab-separated file with columns **`base`** (e.g. `12T`), **`start`**, **`end`** (**1-based** start). Prints run-length counts and highlights runs **≥ 10**. Default `fasta_file` / `output_file` paths are set in `__main__`—change before batch use.

---

### R

#### `Figure1.R`

- **Panel A:** Boxplots of per-read **observed vs estimated** accuracy (`Observed_information.txt`, `Estimated_information.txt`) for human, zebrafish, *K. azureus*, *C. elegans*, *D. melanogaster* (simplex/duplex), *D. hongkongensis* (simplex/duplex).
- **Panel B:** Scatter of **homopolymer overall accuracy** by base from **`Homoploymer_summary.txt`** for the same groups.

Expects `setwd("~/Desktop/benchmark/")` and a Giraffe-style directory tree under that root.

#### `Figure2.R`

Single script driving **six circos PDFs** (*Danio*, human, *K. azureus*, *C. elegans* CGC1, fly, HK2025): gene-type coloring on one track and three homopolymer length bands (4–7, 8–10, &gt;10) on separate tracks. Reads species-specific **BED-like** tables and homopolymer interval files under `circos/`. `setwd("~/Desktop/benchmark/")`.

#### `Figure3.R`

Combined figure workflow: mean **insertion** and **deletion** per read from **`Observed_information.txt`** (HK2025 simplex/duplex, *Drosophila* simplex/duplex); additional sections for homopolymer-related bars, error curves, and accuracy vs homopolymer count (see in-script comments and paths). `setwd("~/Desktop/benchmark/")`.

#### `Figure4.R`

- **Part 1:** Reads homopolymer intervals (e.g. from the Python homopolymer extractor as **`nc031365.txt`**) and a **GFF3**, overlaps with CDS / tRNA / rRNA, optional split of a long Arg-run; stacked bar of homopolymer length class by feature.
- **Part 2:** Duplex **T-homopolymer accuracy vs position** from **`duplex.homopolymer_in_reference.txt`**.
- **Part 3:** HK2025 **BED** + **Homopolish** CSV → genomic ranges, site categories, and summary plots.

Coordinate handling follows comments in-file (e.g. BED start **+1** where noted).

#### `Figure5.R`

Joins **`homopolymer_stats.xlsx`** (species-level A/T/C/G homopolymer count strings) with **`ncbi_lineages_*.csv`** (eukaryote lineages), classifies whether any run is **&gt; 10 bp**, and builds **stacked bar plots** by **phylum** (rare phyla collapsed to “Others”) and follow-on **Nematoda by order**. `setwd("~/Desktop/")`. Uses **Times New Roman** in theme for some panels.

---

## Usage notes

1. **Coordinates:** Homopolymer positions from the Python script are **1-based** on the concatenated FASTA sequence. GFF coordinates follow the file; BED may need **+1** on start where scripts comment.
2. **Giraffe-style inputs:** Figure1 and Figure3 expect paths like `.../2_Observed_quality/Observed_information.txt`, `Homoploymer_summary.txt`, etc.
3. **Homopolish:** Figure4 expects a CSV export with columns matching the script’s `read_csv` / rename logic (`After Homopolish`, `Gene`, etc.); adjust if your export differs.
4. **NCBI dumps:** For `extract_family_genus.py`, obtain current `names.dmp` and `nodes.dmp` from the NCBI taxonomy dump.

---

## Reproducibility

- **Figure4** deduplicates and special-cases specific long homopolymers as documented in code (e.g. split 18T across genes).
- **Figure5** matches taxa by **species** first, then **genus** for unmatched rows.
- Record **Git commit**, input file versions, and any manual edits to path constants when reproducing figures.

---

## Contact

**ZHOU You** — [yzhou799-c@my.cityu.edu.hk](mailto:yzhou799-c@my.cityu.edu.hk)

Issues: `https://github.com/<your-username>/homopolymer-analysis`

---

**Reminder:** Check column names, separators, and **0-based vs 1-based** conventions before batch runs; mismatches are the usual cause of interval errors.
