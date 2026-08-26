# Practical · Part 3 — Differential Expression & Annotation

*Part of the [LISA scRNA-seq practical series](integration_comparison.md) · **Part 3 of 3** · prev → [Part 2: Analysis](practical_part2_analysis.md)*

With integrated clusters from Part 2, this final part turns clusters into **named cell types** and measures the **stimulation response** — first the quick per-cell way, then the statistically correct **pseudobulk** way.

```r
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
set.seed(1234)

ifnb <- readRDS("results/02_ifnb_clustered.rds")
Idents(ifnb) <- "integrated_clusters"
```

---

## 1. Re-join the layers

Differential expression reads from a single expression matrix, so re-join the per-condition layers that were split in Part 2.

```r
ifnb <- JoinLayers(ifnb)
```

---

## 2. Find cluster markers and annotate cell types

First, the genes that mark each cluster (the "who is this cluster?" question):

```r
markers <- FindAllMarkers(
  ifnb,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
top5
```

Check canonical PBMC markers against your clusters before naming anything:

```r
FeaturePlot(ifnb, reduction = "umap.integrated", features = c(
  "CD3D",   # T cells
  "CD8A",   # CD8 T
  "NKG7",   # NK / cytotoxic
  "MS4A1",  # B cells
  "LYZ",    # monocytes
  "FCGR3A", # CD16 monocytes
  "PPBP"    # platelets
), ncol = 3)

DotPlot(ifnb, features = c(
  "IL7R","CCR7","CD3D","CD8A","NKG7","GNLY",
  "MS4A1","CD79A","LYZ","S100A8","FCGR3A","PPBP"
)) + RotatedAxis()
```

Assign labels **only after** the marker evidence agrees. Edit the vector to match *your* cluster order (never paste labels blindly — see the caution in the [main workshop](README.md)):

```r
new_ids <- c("CD14 Mono","CD4 Naive T","CD4 Memory T","CD16 Mono",
             "B","CD8 T","NK","T activated","DC","B activated",
             "Mk","pDC","Eryth")           # adjust to your levels()!
# names(new_ids) <- levels(ifnb)           # uncomment when lengths match
# ifnb <- RenameIdents(ifnb, new_ids)
ifnb$celltype <- Idents(ifnb)              # store the labels in metadata

DimPlot(ifnb, reduction = "umap.integrated", label = TRUE, repel = TRUE) + NoLegend()
```

> `ifnb` also ships with expert annotations in `ifnb$seurat_annotations` — you can compare your labels against them: `table(Idents(ifnb), ifnb$seurat_annotations)`.

---

## 3. The stimulation response — quick per-cell test

Within one cell type, which genes change between STIM and CTRL? The fast (but statistically optimistic) approach compares cells directly:

```r
# pick a cell type present in both conditions
Idents(ifnb) <- "celltype"
cd14 <- subset(ifnb, idents = "CD14 Mono")

de_percell <- FindMarkers(
  cd14,
  group.by = "stim",
  ident.1  = "STIM",
  ident.2  = "CTRL"
)
head(de_percell, 15)
```

Interferon-response genes such as `ISG15`, `IFI6`, `IFIT1` and `ISG20` typically top the list — the biological signal, now measured *within* a cell type instead of confounding the clustering.

Visualize a couple of them:

```r
FeaturePlot(ifnb, reduction = "umap.integrated",
            features = c("ISG15", "IFI6"), split.by = "stim")

VlnPlot(cd14, features = c("ISG15", "IFI6"), group.by = "stim")
```

> **Caveat — this test treats every cell as an independent replicate.** With thousands of cells, p-values become extremely small and even tiny differences look "significant." Good for a quick look; **not** the basis for a condition-level claim. That is what Section 4 fixes.

---

## 4. The correct way — pseudobulk differential expression

For condition-level inference, aggregate counts **per biological replicate** (donor/sample) within each cell type, then run DE on those few aggregated profiles. This respects that **samples/patients — not cells — are the replicates**.

```r
# Aggregate to pseudobulk. Group by the replicate label AND cell type AND condition.
# ifnb pools several donors: use the donor/sample column as the replicate.
# Substitute your own replicate column name here (e.g. "sample" or "orig.ident").
pseudo <- AggregateExpression(
  ifnb,
  assays     = "RNA",
  return.seurat = TRUE,
  group.by   = c("stim", "donor", "celltype")   # <- replicate = donor
)

# Build a condition label on the pseudobulk object and test within one cell type
pseudo$celltype.stim <- paste(pseudo$celltype, pseudo$stim, sep = "_")
Idents(pseudo) <- "celltype.stim"

de_pseudobulk <- FindMarkers(
  pseudo,
  ident.1 = "CD14 Mono_STIM",
  ident.2 = "CD14 Mono_CTRL",
  test.use = "DESeq2"        # count-based test suited to pseudobulk
)
head(de_pseudobulk, 15)
```

> **If your data has no biological replicates** (e.g. one CTRL and one STIM sample, as in a minimal teaching subset), pseudobulk cannot estimate biological variance and a formal condition-level p-value is not meaningful. Report effect sizes and treat conclusions as hypothesis-generating until replicates exist. If `ifnb` in your install lacks a `donor` column, run `colnames(ifnb[[]])` to find the replicate column, or add one.

Compare the two approaches — the per-cell test almost always reports far more "significant" genes than pseudobulk. The pseudobulk list is the trustworthy one:

```r
nrow(de_percell[de_percell$p_val_adj < 0.05, ])
nrow(de_pseudobulk[de_pseudobulk$p_val_adj < 0.05, ])
```

---

## 5. Save results

```r
write.csv(markers,        "results/03_cluster_markers.csv", row.names = FALSE)
write.csv(de_pseudobulk,  "results/03_CD14Mono_STIM_vs_CTRL_pseudobulk.csv")
saveRDS(ifnb,             "results/03_ifnb_annotated.rds")
```

---

## Take-home messages

1. Annotate from **marker evidence**, not from cluster numbers — and keep uncertain labels broad.
2. Measure condition effects **within** aligned cell types, not across the whole dataset.
3. Per-cell DE is a quick look but inflates significance; **pseudobulk by replicate** is the defensible test.
4. **Cells are observations; samples/patients are the replicates.** No amount of cells substitutes for biological replication.

That completes the series: [Part 1 — Setup](practical_part1_setup.md) → [Part 2 — Analysis](practical_part2_analysis.md) → **Part 3 — Differential expression**.

---

### Sources
- [Seurat v5 — Differential expression testing](https://satijalab.org/seurat/articles/de_vignette)
- [Seurat v5 — Integrative analysis (perturbation response)](https://satijalab.org/seurat/articles/integration_introduction)
- [DESeq2 (Love et al., 2014)](https://bioconductor.org/packages/DESeq2/)
- Squair et al. (2021), *Confronting false discoveries in single-cell differential expression*, Nat. Commun. — why pseudobulk matters.
