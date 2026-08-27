# Group Exercise — Analyse Your Own PBMC Dataset

*Companion hands-on for the [LISA scRNA-seq workshop](../README.md). Work in small groups; each group gets its own dataset, analyses it, and presents.*

---

## The idea

You have spent the workshop learning how raw sequencing becomes a gene-by-cell matrix and how Seurat turns that matrix into clusters and cell types. Now you do it yourselves.

Each group receives **one small PBMC dataset** in the standard 10x Genomics format. Your job is to take it through the whole pipeline — QC, clustering, marker genes — decide **what cell types are in your sample**, and **present your answer with the evidence**. Every dataset is a little different, so no two groups will present the same thing.

> **About the data.** These are **simulated** PBMC datasets built for teaching — not real patient samples. They use the real canonical human PBMC marker genes (CD3D, MS4A1, LYZ, NKG7, …), and each one has some deliberately bad cells mixed in for you to remove at the QC step. Everything you learn here transfers directly to real 10x data, which you load exactly the same way.

---

## What each group gets

A folder under `datasets/`:

```
datasets/group_A/   barcodes.tsv.gz   features.tsv.gz   matrix.mtx.gz
datasets/group_B/   ...
datasets/group_C/   ...
datasets/group_D/   ...
```

That trio of files is exactly what `cellranger count` produces in `filtered_feature_bc_matrix/`, so you load it the same way you would load real data: `Read10X("datasets/group_X")`.

Your instructor will tell you which folder is yours.

---

## Timing (about 45 minutes)

| Time | Step |
|---|---|
| 5 min | Load the data, get oriented (how many cells? how many genes?) |
| 10 min | QC — look at the plots, choose thresholds, filter |
| 10 min | Normalize → variable features → PCA → cluster → UMAP |
| 10 min | Marker genes → match to canonical PBMC markers → **name the clusters** |
| 10 min | Build your short presentation |

---

## How to run it

Open **`analyse_your_dataset.R`** in RStudio, change the one line at the top to your group's folder, and run it block by block. **Read the comments and stop at every plot** — the point is to interpret, not just execute. The same steps are summarised below.

### 1. Load

```r
library(Seurat); library(dplyr); library(ggplot2)
counts <- Read10X("datasets/group_A")        # <- your folder
pbmc   <- CreateSeuratObject(counts, project = "group_A",
                             min.cells = 3, min.features = 100)
pbmc
```

*Questions:* How many genes and cells? Is that a plausible number of cells?

### 2. Quality control

```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(pbmc, "nCount_RNA", "percent.mt")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
```

*Your dataset has bad cells seeded in on purpose:* some **damaged** cells (high `percent.mt`) and some **empty-ish** barcodes (very low counts / few genes). A good filter should remove them. *Question:* how many cells did you lose, and does the filtered violin plot look cleaner?

### 3. Standard workflow

```r
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 30)
ElbowPlot(pbmc, ndims = 30)
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)   # also try 0.3 and 0.8
pbmc <- RunUMAP(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
```

*Question:* how many clusters at resolution 0.5? Does changing the resolution merge or split any?

### 4. Markers → identity

```r
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)

DotPlot(pbmc, features = c("CD3D","IL7R","CD8A","NKG7","GNLY","MS4A1","CD79A",
                           "LYZ","S100A8","CD14","FCGR3A","MS4A7",
                           "FCER1A","PPBP","ISG15","IFIT1")) + RotatedAxis()
```

Match each cluster to a lineage using **several markers**, not one:

| Lineage | Look for |
|---|---|
| **CD4 T** | CD3D, IL7R, CCR7 |
| **CD8 T** | CD3D, CD8A, GZMK |
| **NK** | NKG7, GNLY, KLRD1 (and *no* CD3D) |
| **B** | MS4A1, CD79A |
| **CD14 monocyte** | LYZ, S100A8, CD14 |
| **CD16 monocyte** | FCGR3A, MS4A7 (and low CD14) |
| **Dendritic cell** | FCER1A, CLEC10A |
| **Platelet** | PPBP, PF4 |
| **Interferon-activated** | ISG15, IFIT1, IFI6 (on top of a lineage) |

Then name them:

```r
# edit to match YOUR cluster order and evidence
# new_ids <- c("CD4 T","CD14 Mono","CD8 T","B","NK")
# names(new_ids) <- levels(pbmc); pbmc <- RenameIdents(pbmc, new_ids)
```

> **Golden rule (from the main workshop):** annotation is an evidence argument, not a renaming exercise. Use several positive markers, check that *incompatible* markers are absent, and keep uncertain labels broad.

---

## Watch out — your dataset may have a twist

Not every dataset is a plain textbook mix. As you work, stay alert for:

- a **very small cluster** that is easy to overlook but has a crystal-clear marker,
- **two clusters that look almost the same** and differ by only one or two genes,
- **one lineage that appears twice** — same identity, but one version is switched into an activated state.

Part of your presentation is telling us whether your sample had a twist like this and how you spotted it.

---

## Your presentation (5 minutes)

Prepare **five short slides** (or five points on one page). Use the UMAP and DotPlot you saved to `figures/` and `results/`.

1. **Your sample in one line** — how many cells passed QC, how many clusters you found.
2. **QC** — what you filtered and why (show the before/after violin or the mito scatter).
3. **The map** — your annotated UMAP: which cell types, and roughly what proportions.
4. **Your evidence** — for **two** cell types, the markers that convinced you (DotPlot / FeaturePlot).
5. **The twist** — anything surprising: a rare population, two similar subsets, an activation state — and one thing you're *not* sure about.

> A great presentation is honest about uncertainty. "We think cluster 4 is NK because GNLY and NKG7 are high and CD3D is absent, but we're not fully sure it isn't cytotoxic CD8 T" is exactly the right scientific tone.

---

## If you finish early

- Re-run `FindClusters` at resolution 0.3 and 0.8 — what merges or splits? Which resolution best matches the marker evidence?
- Compare your cell-type proportions with another group that had a different dataset.
- For an activation twist: within the affected lineage, run `FindMarkers` between the two states and read off the top genes.

---

## Before you leave — give us feedback

When your group has presented, please take **2 minutes** to fill in the anonymous feedback form. It genuinely shapes the next run of the workshop.

<p align="center">
  <img src="feedback/feedback_qr.png" alt="QR code linking to the feedback form" width="230"><br>
  <strong>Scan the QR code</strong> or open <a href="https://form.jotform.com/262382549939069">form.jotform.com/262382549939069</a>
</p>

*(Instructors: a printable full-page version is in `feedback/feedback_card.html` — open it and put it on the projector as the closing slide.)*

---

## Files in this folder

| File | What it is |
|---|---|
| `README.md` | This worksheet (student-facing). |
| `analyse_your_dataset.R` | Ready-to-run analysis template — change one line to your group. |
| `datasets/group_A … _D/` | The four group datasets in 10x format. |
| `feedback/` | Feedback QR code (`.png`/`.svg`) and a printable closing slide (`feedback_card.html`). |

*(The answer key and solution code live separately, in an instructor-only folder — not here.)*
