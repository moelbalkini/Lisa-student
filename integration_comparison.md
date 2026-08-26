# scRNA-seq Integration Practical — With vs. Without Integration

### A hands-on comparison on public data

This practical extends the [LISA scRNA-seq workshop](README.md). You will download a **public** two-condition PBMC dataset, run the standard Seurat workflow **twice** — once on the merged data with **no integration**, once **with integration** — and compare the results side by side.

The goal is to *see* what a batch/condition effect looks like, and to understand what integration does and does not fix.

> **Uses Seurat v5.** The layer-based workflow (`split`, `IntegrateLayers`, `JoinLayers`) below requires Seurat ≥ 5.0. Check with `packageVersion("Seurat")`.

---

## The dataset

We use **`ifnb`** (Kang et al., 2018): ~14,000 human PBMCs in two conditions — resting **control (CTRL)** and **interferon-β–stimulated (STIM)**. It is a public teaching dataset distributed through the `SeuratData` package, so no manual download or Cell Ranger run is needed.

Because the two conditions were processed separately and stimulation changes many genes, cells cluster **by condition** if you do nothing — a clear, reproducible example of the effect integration is designed to address.

> **Swap in your own data.** Anywhere below, replace the `LoadData("ifnb")` block with your own object (e.g. two samples read with `Read10X` and merged). The rest of the workflow is identical — just make sure the variable that labels batch/sample is stored in the metadata (here it is `stim`).

---

## 1. Install the necessary packages

Run this **once**. It installs from CRAN, Bioconductor and GitHub. `SeuratData` is GitHub-only.

```r
# --- CRAN packages ---
install.packages(c(
  "Seurat",      # core single-cell toolkit (v5+)
  "remotes",     # to install from GitHub
  "ggplot2",     # plotting
  "patchwork",   # arrange plots side by side
  "dplyr",       # data wrangling
  "harmony"      # fast integration method (alternative to CCA)
))

# --- SeuratData: public example datasets (GitHub only) ---
remotes::install_github("satijalab/seurat-data")

# --- Optional but recommended: faster differential expression ---
install.packages("presto")
# if the CRAN install fails, use:
# remotes::install_github("immunogenomics/presto")
```

Load everything and confirm the Seurat version:

```r
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

packageVersion("Seurat")   # should be 5.x
set.seed(1234)             # reproducible clustering/UMAP
```

---

## 2. Get the public data

```r
# Download the dataset once, then load it into memory
InstallData("ifnb")
ifnb <- LoadData("ifnb")

# Older cached objects may need updating to the v5 structure
ifnb <- UpdateSeuratObject(ifnb)

ifnb
table(ifnb$stim)     # how many cells per condition: CTRL vs STIM
```

Split the RNA measurements into **one layer per condition**. This is the v5 mechanism that lets Seurat normalize each batch on its own and later integrate them.

```r
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
```

A light QC filter (thresholds are illustrative — see the QC caution in the main workshop):

```r
ifnb[["percent.mt"]] <- PercentageFeatureSet(ifnb, pattern = "^MT-")
ifnb <- subset(ifnb, subset = nFeature_RNA > 200 &
                              nFeature_RNA < 2500 &
                              percent.mt < 5)
```

---

## 3. Standard preprocessing (shared by both branches)

Normalization, variable-feature selection, scaling and PCA are run the same way regardless of whether you integrate. Integration happens *after* PCA.

```r
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst", nfeatures = 2000)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb, npcs = 30)

ElbowPlot(ifnb, ndims = 30)
```

---

## 4. Branch A — clustering WITHOUT integration

Here we cluster straight off the PCA of the merged data. Nothing corrects for the difference between CTRL and STIM.

```r
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 0.5, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca",
                reduction.name = "umap.unintegrated")

# Colour the SAME embedding two ways: by condition, and by cluster
DimPlot(ifnb, reduction = "umap.unintegrated",
        group.by = c("stim", "unintegrated_clusters"))
```

**What to look for.** In the plot coloured by `stim`, the CTRL and STIM cells sit in **separate regions**. Most clusters are made almost entirely of one condition. The biology (a T cell is a T cell in both conditions) is hidden because the condition difference dominates the top principal components.

```r
# Quantify the split: how mixed is each cluster by condition?
table(ifnb$unintegrated_clusters, ifnb$stim)
```

A cluster that is ~100% CTRL or ~100% STIM is a warning sign that structure is being driven by condition rather than cell identity.

---

## 5. Branch B — clustering WITH integration

`IntegrateLayers` learns a corrected, batch-aware embedding from the per-condition layers. We keep the un-integrated results from Branch A intact and write the integrated embedding to a **new** reduction, so both live in the same object for comparison.

```r
# CCA integration (Seurat's canonical method)
ifnb <- IntegrateLayers(
  object         = ifnb,
  method         = CCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.cca",
  verbose        = FALSE
)
```

Prefer a faster method on larger data? **Harmony** is a drop-in alternative:

```r
# Alternative: Harmony integration
ifnb <- IntegrateLayers(
  object         = ifnb,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.harmony",
  verbose        = FALSE
)
```

Now cluster and embed **on the integrated reduction** (use `"integrated.harmony"` instead if you ran Harmony):

```r
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "integrated.cca")
ifnb <- FindClusters(ifnb, resolution = 0.5, cluster.name = "integrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca",
                reduction.name = "umap.integrated")

DimPlot(ifnb, reduction = "umap.integrated",
        group.by = c("stim", "integrated_clusters"))
```

**What to look for.** Now CTRL and STIM cells **overlap** within each cell-type cluster. Clusters are defined by cell identity, and each contains both conditions — which is what lets you ask *"how does this cell type respond to stimulation?"* rather than just rediscovering that the two samples differ.

```r
table(ifnb$integrated_clusters, ifnb$stim)
```

---

## 6. Compare the two side by side

```r
p1 <- DimPlot(ifnb, reduction = "umap.unintegrated", group.by = "stim") +
        ggtitle("Without integration")
p2 <- DimPlot(ifnb, reduction = "umap.integrated",   group.by = "stim") +
        ggtitle("With integration (CCA)")

p1 + p2      # patchwork places them side by side
```

| | Without integration | With integration |
|---|---|---|
| **UMAP by condition** | CTRL and STIM separate | CTRL and STIM overlap |
| **What drives clusters** | Condition / batch | Cell identity |
| **Cluster composition** | Mostly one condition each | Both conditions in each |
| **Good for** | Seeing the raw batch effect; QC | Cross-condition cell-type comparison, DE |

> **Key message.** Integration aligns shared cell states across samples so that clustering reflects **biology, not batch**. It does not delete the treatment effect — that signal is recovered afterwards through differential expression *within* each aligned cell type.

---

## 7. After integration: recover the biology

Integration is a means, not the end. Re-join the layers and run differential expression **within** a cell type, across conditions — this is where the actual stimulation response is measured.

```r
# Re-join the per-condition layers before expression analysis
ifnb <- JoinLayers(ifnb)

# Example: response to stimulation inside one cell type
# (annotate clusters first, as in the main workshop; here we use a cluster id)
Idents(ifnb) <- "integrated_clusters"
subset_cluster <- subset(ifnb, idents = "0")

response <- FindMarkers(
  subset_cluster,
  group.by = "stim",
  ident.1  = "STIM",
  ident.2  = "CTRL"
)
head(response, 10)   # genes up/down on stimulation within this cell type
```

Interferon-response genes such as `ISG15`, `IFI6` and `IFIT1` typically rise in the STIM group — the biological signal, now measured *within* aligned cell types instead of confounding the clustering.

---

## 8. When should you integrate? (Read before applying to real data)

Integration is powerful and easy to over-apply. From the workshop's cautions:

- **Integrate when** batch/technical/sample differences are burying shared cell types you want to compare across conditions.
- **Be careful when** "batch" and "biology" are the same thing. If each condition genuinely contains different cell populations, aggressive integration can **erase the real difference** you care about.
- **Always inspect the design first.** Look at the un-integrated UMAP (Branch A) *before* deciding — if conditions already mix well, you may not need integration at all.
- **Compare methods.** CCA, RPCA and Harmony can give different results; RPCA/Harmony are more conservative and faster on large data.
- **Cells are not replicates.** For condition-level statistics, use sample-aware or pseudobulk methods — not per-cell tests treating thousands of cells from one patient as independent.

> **Key message.** The un-integrated view is diagnostic, not a mistake to skip. Run both, look at both, and integrate only when the design justifies it.

---

## Appendix — Full script (copy-ready)

```r
library(Seurat); library(SeuratData)
library(ggplot2); library(patchwork); library(dplyr)
set.seed(1234)

# Public data
InstallData("ifnb")
ifnb <- UpdateSeuratObject(LoadData("ifnb"))
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)

# Shared preprocessing
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb, nfeatures = 2000)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb, npcs = 30)

# WITHOUT integration
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 0.5, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# WITH integration (CCA)
ifnb <- IntegrateLayers(ifnb, method = CCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "integrated.cca")
ifnb <- FindClusters(ifnb, resolution = 0.5, cluster.name = "integrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap.integrated")

# Compare
(DimPlot(ifnb, reduction = "umap.unintegrated", group.by = "stim") + ggtitle("Without integration")) +
(DimPlot(ifnb, reduction = "umap.integrated",   group.by = "stim") + ggtitle("With integration"))

# Recover biology
ifnb <- JoinLayers(ifnb)
```

---

## Sources

- [Seurat v5 — Integrative analysis (introduction)](https://satijalab.org/seurat/articles/integration_introduction)
- [Seurat v5 — Data integration methods](https://satijalab.org/seurat/articles/seurat5_integration)
- [SeuratData package](https://github.com/satijalab/seurat-data)
- [Harmony (Korsunsky et al., 2019)](https://github.com/immunogenomics/harmony)
- Kang et al. (2018), *Multiplexed droplet single-cell RNA-sequencing using natural genetic variation*, Nat. Biotechnol. — the `ifnb` dataset.
