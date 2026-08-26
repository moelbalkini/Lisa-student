# Practical · Part 2 — Loading & Downstream Analysis

*Part of the [LISA scRNA-seq practical series](integration_comparison.md) · **Part 2 of 3** · prev → [Part 1: Setup](practical_part1_setup.md) · next → [Part 3: Differential expression](practical_part3_differential_expression.md)*

Here you take the loaded object from Part 1 through the standard Seurat workflow, then run it **twice** — **without** integration and **with** integration — and compare. This is where a batch/condition effect becomes visible and where you decide whether to correct it.

```r
library(Seurat); library(SeuratData)
library(ggplot2); library(patchwork); library(dplyr)
set.seed(1234)

# start from the Part 1 checkpoint (or re-load the data)
ifnb <- readRDS("results/01_ifnb_loaded.rds")
```

---

## 1. Split into per-condition layers and QC

Splitting the RNA measurements into **one layer per condition** is the Seurat v5 mechanism that lets each batch be normalized on its own and later integrated.

```r
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
```

A light QC filter (thresholds are illustrative — see the QC caution in the [main workshop](README.md)):

```r
ifnb[["percent.mt"]] <- PercentageFeatureSet(ifnb, pattern = "^MT-")

VlnPlot(ifnb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, group.by = "stim")

ifnb <- subset(ifnb, subset = nFeature_RNA > 200 &
                              nFeature_RNA < 2500 &
                              percent.mt < 5)
```

---

## 2. Standard preprocessing (shared by both branches)

Normalization, variable features, scaling and PCA are identical whether or not you integrate. **Integration happens *after* PCA.**

```r
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst", nfeatures = 2000)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb, npcs = 30)

ElbowPlot(ifnb, ndims = 30)   # how many PCs carry real signal?
```

---

## 3. Branch A — WITHOUT integration

Cluster straight off the merged PCA. Nothing corrects for the CTRL vs STIM difference.

```r
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 0.5, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca",
                reduction.name = "umap.unintegrated")

# same embedding, coloured two ways
DimPlot(ifnb, reduction = "umap.unintegrated",
        group.by = c("stim", "unintegrated_clusters"))
```

**What to look for.** Coloured by `stim`, CTRL and STIM cells sit in **separate regions**; most clusters are almost entirely one condition. The biology (a T cell is a T cell in both) is hidden because the condition difference dominates the top PCs.

```r
# how mixed is each cluster by condition?
table(ifnb$unintegrated_clusters, ifnb$stim)
```

A cluster that is ~100% one condition is a warning sign that structure is driven by condition, not cell identity.

---

## 4. Branch B — WITH integration

`IntegrateLayers` learns a batch-aware embedding from the per-condition layers and writes it to a **new** reduction, so Branch A stays intact for comparison.

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

Faster alternative on larger data — **Harmony** (drop-in):

```r
ifnb <- IntegrateLayers(
  object = ifnb, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)
```

Cluster and embed **on the integrated reduction** (swap in `"integrated.harmony"` if you used Harmony):

```r
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "integrated.cca")
ifnb <- FindClusters(ifnb, resolution = 0.5, cluster.name = "integrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca",
                reduction.name = "umap.integrated")

DimPlot(ifnb, reduction = "umap.integrated",
        group.by = c("stim", "integrated_clusters"))
```

**What to look for.** CTRL and STIM cells now **overlap** within each cell-type cluster. Clusters are defined by identity, and each holds both conditions — which is what lets you later ask *"how does this cell type respond to stimulation?"*

```r
table(ifnb$integrated_clusters, ifnb$stim)
```

---

## 5. Compare the two side by side

```r
p1 <- DimPlot(ifnb, reduction = "umap.unintegrated", group.by = "stim") +
        ggtitle("Without integration")
p2 <- DimPlot(ifnb, reduction = "umap.integrated",   group.by = "stim") +
        ggtitle("With integration (CCA)")
p1 + p2
```

| | Without integration | With integration |
|---|---|---|
| **UMAP by condition** | CTRL and STIM separate | CTRL and STIM overlap |
| **What drives clusters** | Condition / batch | Cell identity |
| **Cluster composition** | Mostly one condition each | Both conditions in each |
| **Good for** | Seeing the raw batch effect; QC | Cross-condition cell-type comparison, DE |

> **Key message.** Integration aligns shared cell states across samples so clustering reflects **biology, not batch**. It does not delete the treatment effect — that signal is recovered afterwards (Part 3) through differential expression *within* each aligned cell type.

---

## 6. When should you integrate?

- **Integrate when** batch/technical/sample differences bury shared cell types you want to compare.
- **Be careful when** "batch" and "biology" coincide — aggressive integration can erase a real difference.
- **Always inspect Branch A first.** If conditions already mix, you may not need integration.
- **Compare methods.** CCA, RPCA and Harmony differ; RPCA/Harmony are more conservative and faster.

---

## Checkpoint — save for Part 3

```r
saveRDS(ifnb, "results/02_ifnb_clustered.rds")
```

**Next:** [Part 3 — Differential expression & annotation](practical_part3_differential_expression.md), where we name the clusters and measure the stimulation response the right way.

---

### Sources
- [Seurat v5 — Integrative analysis (introduction)](https://satijalab.org/seurat/articles/integration_introduction)
- [Seurat v5 — Data integration methods](https://satijalab.org/seurat/articles/seurat5_integration)
- [Harmony (Korsunsky et al., 2019)](https://github.com/immunogenomics/harmony)
