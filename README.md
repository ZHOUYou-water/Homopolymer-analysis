# Homopolymer analysis (Nanopore)

Scripts for **read-level QC**, **ONT mitochondrial assembly and polishing**, **homopolymer accuracy**, **mitochondrial annotation context**, **circos layouts**, and **comparative homopolymer statistics**, using outputs from alignment / quality tools (e.g. Giraffe-style `Observed_information.txt`, `Homoploymer_summary.txt`, `homopolymer_in_reference` tables), **GFF3 / BED**, **Homopolish** site lists, and optional **NCBI taxonomy** dumps.

---

## Table of contents

- [Overview](#overview)
- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Installation](#installation)
- [Scripts](#scripts)
  - [Assembly pipeline (Shell + Python)](#assembly-pipeline-shell--python)
  - [Python](#python)
  - [R](#r)
- [Usage notes](#usage-notes)
- [Reproducibility](#reproducibility)
- [Data availability](#data-availability)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

This repository contains **1 Shell**, **4 Python**, and **5 R** scripts that together support:

- **Read-level QC**: observed vs estimated per-read accuracy (Giraffe-style tables).
- **Mitochondrial assembly**: FASTQ filtering, Flye assembly, Racon and Medaka polishing, optional replicate selection by indel count vs a reference (e.g. PX108868).
- **Homopolymer accuracy**: species- and library-type summaries from `Homoploymer_summary.txt`; point accuracy along the genome (e.g. T homopolymers).
- **Mitochondrial context**: homopolymer intervals overlapped with CDS / tRNA / rRNA (GFF3); circos gene + homopolymer tracks for multiple references; Homopolish sites on a gene model.
- **Error / indel panels**: mean insertion and deletion per read; homopolymer-class and position-related panels for selected datasets.
- **Taxonomy-linked summaries**: family–genus lists from NCBI `names.dmp` / `nodes.dmp`; phylum-level stacked bars from species homopolymer count spreadsheets plus NCBI lineage tables.

> **Note:** Paths use `setwd`, `~/Desktop/`, and other fixed locations. Before publishing or sharing, point these to your data layout or use project-relative paths (e.g. [`here`](https://here.r-lib.org/)).

### Overall pipeline

```
raw FASTQ / BAM
      │
      ├─► [assembly/] QC → Flye → Racon → Medaka → draft.fasta
      │         │
      │         └─► homopolymer extraction → GFF overlap → figures
      │
      └─► [Giraffe] Observed_information / homopolymer_in_reference → Figure1–5
```

[Pipeline.tif](https://github.com/user-attachments/files/26134435/Pipeline.tif)

---

## Repository layout

```text
code/
├── README.md
├── requirements.txt
├── assembly/
│   ├── mt_pipeline.sh            # Main: QC → Flye → Racon → Medaka
│   ├── filter_fastq.py           # FASTQ filter (length + mean Q)
│   └── score_vs_ref.py           # MUMmer indel/SNP scoring vs reference
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

### Assembly pipeline (external binaries)

| Tool | Purpose | Install |
|------|---------|---------|
| [Flye](https://github.com/fenderglass/Flye) ≥ 2.9 | Assembly | `conda install -c bioconda flye` |
| [minimap2](https://github.com/lh3/minimap2) | Read mapping | `conda install -c bioconda minimap2` |
| [Racon](https://github.com/lbcb-sci/racon) | Polishing | `conda install -c bioconda racon` |
| [Medaka](https://github.com/nanoporetech/medaka) | Polishing | `conda install -c bioconda medaka` |
| [MUMmer](https://github.com/mummer4/mummer) | Optional: score vs reference | `conda install -c bioconda mummer` |

Recommended conda environment:

```bash
conda create -n ont-mt -c bioconda -c conda-forge \
  flye minimap2 racon medaka mummer python=3.10 pandas openpyxl
conda activate ont-mt
```

### Python

- **Version:** 3.8+
- **Packages:** see `requirements.txt`
  - `extract_family_genus.py`: **pandas**, **openpyxl** (Excel export)
  - `extract_number of homo.py`: standard library only
  - `assembly/filter_fastq.py`, `assembly/score_vs_ref.py`: standard library only

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
cd homopolymer-analysis/code
python3 -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
chmod +x assembly/mt_pipeline.sh assembly/*.py
```

---

## Scripts

### Assembly pipeline (Shell + Python)

End-to-end ONT mitochondrial genome workflow used in *D. asiatica* projects (e.g. PX108868 as reference).

#### `assembly/mt_pipeline.sh`

Main pipeline:

```
raw FASTQ → [QC] → [Flye ×N] → [Racon ×2] → [Medaka] → final assembly
```

| Step | Tool | Default |
|------|------|---------|
| QC | `filter_fastq.py` | length ≥ 1000 bp, mean Q ≥ 7 |
| Assembly | Flye | `--genome-size 15k --min-overlap 1000 --meta` |
| Polish | Racon | 2 rounds (`minimap2` + `racon`) |
| Polish | Medaka | model `r1041_e82_400bps_sup_v5.0.0` |

**Example** (6 Flye replicates, pick draft with fewest indels vs PX108868):

```bash
cd assembly
./mt_pipeline.sh \
  -i /path/to/reads.fastq \
  -o /path/to/output \
  -r /path/to/PX108868.fasta \
  --flye-runs 6 \
  -t 16
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | required | Input FASTQ |
| `-o, --outdir` | required | Output directory |
| `-r, --reference` | none | Reference FASTA for picking best Flye replicate |
| `--min-length` | 1000 | Minimum read length |
| `--min-mean-q` | 7 | Minimum mean Phred quality |
| `--flye-runs` | 1 | Number of Flye replicates |
| `--racon-rounds` | 2 | Racon polishing rounds |
| `--skip-qc / --skip-flye / --skip-racon / --skip-medaka` | | Skip individual steps |

**Output:**

```text
output/
├── filter_q7_1000.fastq
├── flye_run_1/ ... flye_run_N/
├── draft.fasta
├── racon.fasta
├── medaka.fasta          # final assembly
├── flye_scores.tsv       # per-replicate indel counts (if --reference)
└── pipeline.log
```

> **Note:** For Dorado `calls.bam` (HAC/SUP basecalls), extract duplex reads (`dx:i:1`) and convert to FASTQ before running this pipeline. For manually corrected drafts, consider `--racon-rounds 1` to avoid over-truncation.

#### `assembly/filter_fastq.py`

Standalone FASTQ filter (same logic as pipeline step 1):

```bash
python3 filter_fastq.py -i reads.fastq -o filter_q7_1000.fastq -l 1000 -q 7
```

#### `assembly/score_vs_ref.py`

Score an assembly against a reference with MUMmer (`nucmer` + `show-snps`):

```bash
python3 score_vs_ref.py \
  --ref PX108868.fasta \
  --assembly draft.fasta \
  --mummer /path/to/mummer-4.0.0rc1
```

Outputs indel / SNP counts to `*.scores.tsv`. Used internally to select the best Flye replicate and to compare polishing steps (see Supplementary Table 6 in the manuscript).

---

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

Single script driving **six circos PDFs** (*Danio*, human, *K. azureus*, *C. elegans* CGC1, fly, HK2025): gene-type coloring on one track and three homopolymer length bands (4–7, 8–10, >10) on separate tracks. Reads species-specific **BED-like** tables and homopolymer interval files under `circos/`. `setwd("~/Desktop/benchmark/")`.

#### `Figure3.R`

Combined figure workflow: mean **insertion** and **deletion** per read from **`Observed_information.txt`** (HK2025 simplex/duplex, *Drosophila* simplex/duplex); additional sections for homopolymer-related bars, error curves, and accuracy vs homopolymer count (see in-script comments and paths). `setwd("~/Desktop/benchmark/")`.

#### `Figure4.R`

- **Part 1:** Reads homopolymer intervals (e.g. from the Python homopolymer extractor as **`nc031365.txt`**) and a **GFF3**, overlaps with CDS / tRNA / rRNA, optional split of a long Arg-run; stacked bar of homopolymer length class by feature.
- **Part 2:** Duplex **T-homopolymer accuracy vs position** from **`duplex.homopolymer_in_reference.txt`**.
- **Part 3:** HK2025 **BED** + **Homopolish** CSV → genomic ranges, site categories, and summary plots.

Coordinate handling follows comments in-file (e.g. BED start **+1** where noted).

#### `Figure5.R`

Joins **`homopolymer_stats.xlsx`** (species-level A/T/C/G homopolymer count strings) with **`ncbi_lineages_*.csv`** (eukaryote lineages), classifies whether any run is **> 10 bp**, and builds **stacked bar plots** by **phylum** (rare phyla collapsed to “Others”) and follow-on **Nematoda by order**. `setwd("~/Desktop/")`. Uses **Times New Roman** in theme for some panels.

---

## Usage notes

1. **Assembly → homopolymer analysis:** Run `assembly/mt_pipeline.sh` first to obtain `medaka.fasta`, then feed the assembly into `extract_number of homo.py` and downstream R scripts (Figure2–4).
2. **Coordinates:** Homopolymer positions from the Python script are **1-based** on the concatenated FASTA sequence. GFF coordinates follow the file; BED may need **+1** on start where scripts comment.
3. **Giraffe-style inputs:** Figure1 and Figure3 expect paths like `.../2_Observed_quality/Observed_information.txt`, `Homoploymer_summary.txt`, etc.
4. **Homopolish:** Figure4 expects a CSV export with columns matching the script’s `read_csv` / rename logic (`After Homopolish`, `Gene`, etc.); adjust if your export differs.
5. **NCBI dumps:** For `extract_family_genus.py`, obtain current `names.dmp` and `nodes.dmp` from the NCBI taxonomy dump.
6. **Polishing benchmarks:** Use `score_vs_ref.py` to quantify indels vs PX108868 at each polishing step (draft, racon1, racon2, medaka) for supplementary tables.

---

## Reproducibility

- **Assembly:** Record Flye version, Medaka model, read filter thresholds, and number of replicates. Log written to `output/pipeline.log`.
- **Figure4** deduplicates and special-cases specific long homopolymers as documented in code (e.g. split 18T across genes).
- **Figure5** matches taxa by **species** first, then **genus** for unmatched rows.
- Record **Git commit**, input file versions, and any manual edits to path constants when reproducing figures.

---

## Data availability

*(Add accession numbers, Zenodo DOI, or supplementary file links here.)*

---

## Citation

If you use this pipeline, please cite the underlying tools: **Flye**, **Racon**, **Medaka**, **minimap2**, and any Giraffe / Homopolish references used in your analysis.

---

## Contact

**ZHOU You** — [yzhou799-c@my.cityu.edu.hk](mailto:yzhou799-c@my.cityu.edu.hk)

Issues: `https://github.com/<your-username>/homopolymer-analysis`

---

**Reminder:** Check column names, separators, and **0-based vs 1-based** conventions before batch runs; mismatches are the usual cause of interval errors.
