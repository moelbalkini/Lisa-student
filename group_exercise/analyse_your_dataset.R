# =============================================================================
# LISA scRNA-seq -- Group exercise analysis template
# -----------------------------------------------------------------------------
# Each group: change GROUP below to your assigned folder (group_A / _B / _C / _D)
# then run this script top to bottom in RStudio, reading the comments as you go.
# This is a TEMPLATE, not a black box -- stop at each plot and interpret it.
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
set.seed(1234)

# ---- 0. Point to YOUR group's dataset --------------------------------------
GROUP    <- "group_A"                                   # <-- CHANGE THIS
data_dir <- file.path("datasets", GROUP)                # run from group_exercise/
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 1. Load the 10x matrix ------------------------------------------------
counts <- Read10X(data.dir = data_dir)
dim(counts)                                             # genes x cells

pbmc <- CreateSeuratObject(counts = counts, project = GROUP,
                           min.cells = 3, min.features = 100)
pbmc

# ---- 2. Quality control ----------------------------------------------------
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# LOOK at the distributions for THIS dataset before choosing thresholds:
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(pbmc, "nCount_RNA", "percent.mt")
FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA")

# Sensible starting thresholds for these teaching datasets
# (a batch of damaged high-mito cells and empty-ish low-count barcodes were
#  seeded in on purpose -- your filter should remove them):
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
pbmc                                                    # how many cells remain?

# ---- 3. Normalize, variable features, scale, PCA ---------------------------
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 30)
ElbowPlot(pbmc, ndims = 30)                             # how many PCs carry signal?

# ---- 4. Cluster and UMAP ---------------------------------------------------
dims_use <- 1:15
pbmc <- FindNeighbors(pbmc, dims = dims_use)
pbmc <- FindClusters(pbmc, resolution = 0.5)           # try 0.3 / 0.5 / 0.8 too
pbmc <- RunUMAP(pbmc, dims = dims_use)
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()

# ---- 5. Cluster markers ----------------------------------------------------
markers <- FindAllMarkers(pbmc, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
print(top5, n = 100)

# Canonical PBMC markers -- which cluster lights up for each?
canon <- c("CD3D","IL7R","CCR7",        # T (CD4 naive)
           "CD8A","GZMK",               # CD8 T
           "NKG7","GNLY","KLRD1",       # NK
           "MS4A1","CD79A",             # B
           "LYZ","S100A8","CD14",       # CD14 monocytes
           "FCGR3A","MS4A7",            # CD16 monocytes
           "FCER1A","CLEC10A",          # dendritic cells
           "PPBP","PF4",                # platelets
           "ISG15","IFIT1","IFI6")      # interferon-response program
canon <- canon[canon %in% rownames(pbmc)]
DotPlot(pbmc, features = canon) + RotatedAxis()
FeaturePlot(pbmc, features = c("CD3D","LYZ","MS4A1","NKG7"), ncol = 2)

# ---- 6. Name your clusters (ONLY after the markers agree) ------------------
# Example -- edit the vector to match YOUR clusters and marker evidence:
# new_ids <- c("CD4 T", "CD14 Mono", "CD8 T", "B", "NK")
# names(new_ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new_ids)
# pbmc$celltype <- Idents(pbmc)
# DimPlot(pbmc, label = TRUE, repel = TRUE) + NoLegend()

# ---- 7. Save your results for the presentation ----------------------------
ggsave(file.path("figures", paste0(GROUP, "_umap.png")),
       DimPlot(pbmc, label = TRUE) + NoLegend(), width = 6, height = 5, dpi = 150)
write.csv(top5, file.path("results", paste0(GROUP, "_top_markers.csv")), row.names = FALSE)
saveRDS(pbmc, file.path("results", paste0(GROUP, "_seurat.rds")))
