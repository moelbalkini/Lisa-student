# Practical · Part 1 — Setup & Installation

*Part of the [LISA scRNA-seq practical series](integration_comparison.md) · **Part 1 of 3** · next → [Part 2: Loading & downstream analysis](practical_part2_analysis.md)*

This part gets your environment ready and downloads the public dataset. By the end you will have every package installed and a Seurat object loaded, ready for analysis in Part 2.

---

## What you need first

- **R ≥ 4.3** and **RStudio** (recommended).
- Internet access for the one-time package and data download.
- ~4 GB free RAM for this teaching dataset (real projects need much more).

We use **Seurat v5**, whose layer-based workflow is required for the integration steps in Part 2.

---

## 1. Install the necessary packages

Run this block **once**. It pulls from CRAN, plus one GitHub-only package (`SeuratData`).

```r
# --- CRAN packages ---
install.packages(c(
  "Seurat",      # core single-cell toolkit (v5+)
  "remotes",     # to install packages from GitHub
  "ggplot2",     # plotting
  "patchwork",   # arrange plots side by side
  "dplyr",       # data wrangling
  "harmony"      # fast integration method (used in Part 2)
))

# --- SeuratData: public example datasets (GitHub only) ---
remotes::install_github("satijalab/seurat-data")

# --- Optional: much faster differential expression (used in Part 3) ---
install.packages("presto")
# fallback if the CRAN build fails:
# remotes::install_github("immunogenomics/presto")
```

> **If an install fails**, read the first error line. The usual causes are a missing system library (install it via your OS package manager) or an out-of-date R. On institutional machines you may need to set a personal library with `.libPaths()`.

---

## 2. Load the packages and confirm versions

```r
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

packageVersion("Seurat")   # must be 5.x for this series
set.seed(1234)             # makes clustering / UMAP reproducible
```

If `packageVersion("Seurat")` shows a 4.x version, update with `install.packages("Seurat")` and restart R before continuing.

---

## 3. Download the public data

We use **`ifnb`** (Kang et al., 2018): ~14,000 human PBMCs in two conditions — resting **control (CTRL)** and **interferon-β–stimulated (STIM)**. It ships through `SeuratData`, so there is no manual download or Cell Ranger run.

```r
InstallData("ifnb")                       # downloads once, caches locally
ifnb <- LoadData("ifnb")                  # load into memory
ifnb <- UpdateSeuratObject(ifnb)          # ensure the v5 object structure

ifnb
table(ifnb$stim)                          # cells per condition: CTRL vs STIM
```

You should see two conditions with several thousand cells each. The `stim` metadata column is the **batch/condition label** we will compare across in Part 2.

> **Using your own data instead?** Replace the three lines above with your own object — e.g. read each sample with `Read10X()`, `CreateSeuratObject()` for each, then `merge()` them, making sure a metadata column records which sample each cell came from. Everything in Parts 2–3 then works the same, using your sample column in place of `stim`.

---

## 4. Set up a project folder (optional but recommended)

Keeping inputs, scripts and results separate makes the analysis reproducible.

```r
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Save the freshly loaded object so Part 2 can start from a clean checkpoint
saveRDS(ifnb, "results/01_ifnb_loaded.rds")
```

---

## Checkpoint

You now have:

- All packages installed and loading without error.
- Seurat reporting version 5.x.
- The `ifnb` object in memory (and optionally saved to `results/01_ifnb_loaded.rds`).

**Next:** [Part 2 — Loading & downstream analysis](practical_part2_analysis.md), where we run QC, the standard Seurat workflow, and the with/without-integration comparison.

---

### Sources
- [SeuratData package](https://github.com/satijalab/seurat-data)
- [Seurat installation](https://satijalab.org/seurat/articles/install)
- Kang et al. (2018), *Multiplexed droplet single-cell RNA-sequencing using natural genetic variation*, Nat. Biotechnol. — the `ifnb` dataset.
