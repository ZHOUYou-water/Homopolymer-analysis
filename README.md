# Homopolymer analysis (Nanopore)

Scripts for analyzing homopolymer-associated errors and related QC in **mitochondrial and comparative genomics** workflows, using outputs from alignments / quality tools (e.g. Giraffe-style `Observed_information.txt`, `Homoploymer_summary.txt`) and **Homopolish** site lists.

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

This collection supports:

- **Read-level QC**: observed vs estimated accuracy; insertion/deletion summaries from per-read tables.
- **Homopolymer accuracy**: accuracy vs homopolymer length across species or library type (simplex / duplex); summaries from `Homoploymer_summary.txt`-style files.
- **Mitochondrial context**: overlap of homopolymers with CDS / tRNA / rRNA (GFF3 + BED); lollipop-style site maps.
- **Genome layout**: circos tracks for genes + homopolymer segments; multi-genome length / feature comparison from FASTA + GFF3.
- **Assembly / error panels**: homopolymer-class error rates (e.g. HK2025-style count tables) for figures.

> **Note:** Scripts were developed with local absolute paths (`setwd`, input paths). Before publishing or sharing, factor paths into a small config section or RStudio project-relative paths.

---

## Repository layout

Copy or symlink your analysis scripts into a structure similar to:

```text
homopolymer-analysis/
├── README.md
├── requirements.txt
├── python/
│   ├── count_homopolymers.py      # FASTA-level homopolymer enumeration
│   ├── count_cds_homo.py          # Homopolish sites → CDS mapping + summaries
│   └── stat_homo_in_Cds.py        # Per-CDS short vs long homopolymer counts
└── r/
    ├── figure1_homo.R
    ├── figure2_homo.R
    ├── new_revise_figure1_homo.R
    ├── new_revise_figure1_qc.R
    ├── figure1_qc.R
    ├── figure2_qc.R
    ├── figure4_assembly_acc.R
    ├── point_accurcy.R
    ├── read_insertion_deletion.R
    ├── count_mito_homo.R
    ├── lolipop.R
    ├── genome_comparison.R
    ├── error.R
    ├── count_length.R
    ├── human_circos_revised.R
    ├── fly_circos_revised.R
    ├── danio_circos_revised.R
    ├── k_circos_revised.R
    ├── cgc1_circos_revised.R
    ├── hk2025_circos_revised.R
    └── code/
        ├── Figure1.R
        ├── Figure2.R
        ├── Figure3.R
        ├── Figure4.R
        └── Figure5.R
```

Adjust names to match what you actually commit; the descriptions below follow the current script logic.

---

## Requirements

### Python

- **Version:** 3.8+
- **Packages:** see `requirements.txt` (at minimum `pandas` for CDS / Homopolish integration).

### R

Typical packages used across scripts (versions may vary):

| Package        | Typical use                          |
|----------------|--------------------------------------|
| `ggplot2`      | Core plotting                        |
| `dplyr`, `tidyr`, `readr`, `stringr` | Data wrangling (`tidyverse` meta) |
| `patchwork`    | Multi-panel figures                  |
| `scales`       | Axis / label formatting              |
| `ggdist`       | Distribution / raincloud-style plots |
| `GenomicRanges`, `rtracklayer` | Interval overlap (GFF/BED/sites) |
| `circlize`     | Circos genome + homopolymer tracks   |

Install missing packages in R:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr",
  "patchwork", "scales", "tidyverse", "ggdist", "circlize"
))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

---

## Installation

```bash
git clone https://github.com/<your-username>/homopolymer-analysis.git
cd homopolymer-analysis
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

---

## Scripts

### Python

#### `count_homopolymers.py`

Scans one or more sequences in a **FASTA** file, reports homopolymer counts by base and by run length (configurable minimum length).

```bash
python3 count_homopolymers.py reference.fasta -m 2
python3 count_homopolymers.py reference.fasta -m 2 -v   # per-site lines
```

#### `count_cds_homo.py` (Homopolish → CDS)

Maps **Homopolish** correction sites to **CDS** coordinates using a CDS FASTA whose headers encode genomic intervals, e.g. `>gene|start-end`. Duplicate genome positions are collapsed (**one row per `Before_Pos`**). Writes:

- Per-site table with CDS-relative position and local homopolymer length.
- Per-gene summary (counts in homopolymer, length &gt; 10, etc.).

Edit the file-header paths for: Homopolish site list, CDS FASTA, and output CSV paths.

#### `stat_homo_in_Cds.py`

Reads a CDS multi-FASTA and, per sequence, counts homopolymer runs of length ≥ 4, split into **&lt; 10** vs **≥ 10** bp. Prints a TSV-style table to stdout. Paths are set inside the script—parameterize or use `argparse` if you generalize for the repo.

---

### R

#### Homopolymer accuracy and read QC

| Script | Role |
|--------|------|
| `figure1_homo.R` | Multi-species lines: homopolymer length vs accuracy (%), and vs log2 count; facet by species. |
| `figure2_homo.R` | Simplex/duplex (e.g. *Drosophila*, *Dirofilaria*): A/T-focused homopolymer panels. |
| `new_revise_figure1_homo.R` | Revised multi-dataset figure from `Homoploymer_summary.txt` inputs + `ggdist` / `patchwork`. |
| `figure1_qc.R`, `figure2_qc.R`, `new_revise_figure1_qc.R` | Read accuracy: observed vs estimated from Giraffe-style `Observed_information.txt` / `Estimated_information.txt`. |
| `code/Figure1.R`–`Figure5.R` | Publication-style composite panels (QC, homopolymer, etc.—see in-file comments and `setwd`). |
| `point_accurcy.R` | Point-level accuracy from `homopolymer_in_reference`-style tables (e.g. T homopolymers, genomic segments). |
| `read_insertion_deletion.R` | Insertion/deletion summaries from per-read observed quality tables. |
| `error.R`, `count_length.R` | Supporting error / length summaries (paths inside script). |

#### Mitochondrial genome, genes, and Homopolish sites

| Script | Role |
|--------|------|
| `count_mito_homo.R` | Homopolymer intervals + GFF3 features → overlap with CDS/tRNA/rRNA; optional special-case splitting of long runs. |
| `lolipop.R` | BED gene model + Homopolish sites CSV → genomic ranges and lollipop-style visualization. |

#### Circos (genome + homopolymers)

Species- or dataset-specific drivers (same pattern): `human_circos_revised.R`, `fly_circos_revised.R`, `danio_circos_revised.R`, `k_circos_revised.R`, `cgc1_circos_revised.R`, `hk2025_circos_revised.R`. Each reads a `.bed` (gene track) and a homopolymer interval file, then draws circos PDFs via `circlize`.

#### Assembly / comparison

| Script | Role |
|--------|------|
| `figure4_assembly_acc.R` | Homopolymer-class error rates (e.g. simplex vs duplex) from summarized count tables. |
| `genome_comparison.R` | Compare genome lengths and GFF features across references (horizontal layout figure). |

---

## Usage notes

1. **Coordinates:** Scripts assume **1-based** genomic positions where noted (e.g. Homopolish `After_Pos`, GFF-derived intervals). BED inputs may be converted with `+1` on start where the script comments specify 0-based BED.
2. **Giraffe-style inputs:** Many R scripts expect directory trees like `.../2_Observed_quality/Observed_information.txt`, `Homoploymer_summary.txt`, etc. Mirror that layout or change `read.table` paths.
3. **Homopolish:** Site files may be whitespace-separated columns (`Before_Pos`, bases, `After_Pos`) or CSV exports—use the Python script that matches your format (`count_cds_homo.py` vs workflows in `lolipop.R`).
4. **Paths:** Replace `setwd("/Users/...")` and hard-coded paths with your machine or use [`here`](https://here.r-lib.org/) / project-relative paths.

---

## Reproducibility

- **Duplicate genomic positions:** Homopolish → CDS mapping deduplicates on `Before_Pos` (one correction per position).
- **Homopolymer length** in CDS mapping is computed **locally** along the CDS sequence at the mapped coordinate.
- **Analysis parameters** (length thresholds, species subsets, facet layouts) are set explicitly in each script; record Git commit hashes and random seeds if you add stochastic steps later.

---

## Data availability

Raw sequencing data are not included (size / ethics / access). Processed tables compatible with these scripts (e.g. Giraffe outputs, homopolymer summary TSVs) can be shared per publication policy or on request.

---

## Citation

If you use this repository in a publication, please cite the associated paper (add DOI / bioRxiv link when available) and, where applicable, tools such as **Homopolish**, **Giraffe**, and **Oxford Nanopore** basecalling / chemistry documentation.

---

## Contact

**ZHOU You** — [yzhou799-c@my.cityu.edu.hk](mailto:yzhou799-c@my.cityu.edu.hk)

For bugs or feature requests, please open an issue on GitHub (replace with your real repo URL after upload):

`https://github.com/<your-username>/homopolymer-analysis`

---

**Reminder:** Verify input column names, separators, and coordinate systems (0-based vs 1-based) before batch runs; mismatches are the most common source of off-by-one overlaps in genomics scripts.
