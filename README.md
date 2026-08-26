# LISA scRNA-seq — Student Practical Guide

### From Wet Lab to Biological Insight

A beginner-friendly introduction to single-cell RNA sequencing with **Cell Ranger**, **R** and **Seurat**.

| | |
|---|---|
| **Audience** | Life-science students with little or no bioinformatics experience |
| **Format** | Guided demonstration + hands-on R analysis (~2 hours) |
| **Instructors** | Mohamed Elbalkini · Imke Hinrichs · Abdus Salam — MHH / TWINCORE |

> This is the student handout. It contains the concepts, commands, exercises and troubleshooting notes you need to follow along. Code blocks can be copied straight into the terminal or RStudio. The full instructor script is in [`LISA_scRNA_Workshop_WetLab_to_Seurat.docx`](LISA_scRNA_Workshop_WetLab_to_Seurat.docx).

> **Companion practical series** (public data, hands-on, Seurat v5): [**overview**](integration_comparison.md) → [Part 1 — Setup](practical_part1_setup.md) → [Part 2 — Analysis (with vs. without integration)](practical_part2_analysis.md) → [Part 3 — Differential expression](practical_part3_differential_expression.md).

---

## Learning outcomes

By the end of this workshop you should be able to:

- Explain the journey from tissue or blood to a gene-by-cell count matrix.
- Describe cell barcodes, UMIs, reads, FASTQ files and alignment in biological language.
- Recognize the four common 10x FASTQ files and understand why R1 and R2 have different jobs.
- Run or explain a basic `cellranger count` command and identify its key outputs.
- Import a filtered 10x matrix into Seurat, and perform QC, normalization, PCA, clustering and UMAP.
- Use marker genes to propose cell identities while recognizing uncertainty and common artifacts.

---

## Part 1 — Start with the biological question

Begin with a concrete example: a PBMC sample from a patient *before* treatment and another *after* treatment. Questions you might ask include: Which immune cell types are present? Do their proportions change? Which cell state responds? Which pathways are activated?

> **Key message.** The analysis method must follow the biological question. Clustering is not the question; it is a tool for organizing cells so that biological patterns can be examined.

**Activity.** With a partner, write **one** biological question and identify the *unit of comparison*: cells, cell types, samples, patients or conditions. Remember: biological replicates are patients or samples — **not** thousands of cells from one patient.

### What scRNA-seq can and cannot tell you

Single-cell RNA sequencing samples the RNA molecules captured from individual cells. It gives a sparse, noisy snapshot of transcriptional state.

| It can support | It cannot prove alone |
|---|---|
| Cell-type composition and transcriptional states | Causal mechanism |
| Differential expression between well-designed groups | Protein abundance or cytokine secretion |
| Candidate pathways and cellular interactions | Physical cell migration without spatial/paired evidence |
| Rare populations when sampling is adequate | Absence of a population when capture depth is low |

---

## Part 2 — Wet lab to barcoded molecules

Before opening R, it helps to know what information was physically introduced in the wet lab — most downstream assumptions begin here.

### The physical workflow (and what can bias it)

| Stage | What happens | What can bias the result |
|---|---|---|
| Sample collection | Blood or tissue is collected and transported | Delay, temperature, anticoagulant, ischemia |
| Cell preparation | A single-cell suspension is produced | Stress response, selective cell loss, doublets |
| Counting and viability | Cell concentration and viability are measured | Incorrect loading concentration, debris |
| 10x partitioning | Cells and gel beads enter droplets/GEMs | Empty droplets, two cells in one droplet |
| Reverse transcription | RNA is copied to barcoded cDNA | Capture efficiency and transcript length |
| Library preparation | Amplifiable fragments are generated | PCR bias and low library complexity |
| Sequencing | Libraries become reads with quality scores | Insufficient depth or poor base quality |

### The three labels that make single-cell counting possible

Imagine every RNA molecule receives two stickers: one says *which cell* it came from, the other is a *serial number* for that original molecule. PCR then makes many photocopies — but copies that share the same gene, cell barcode and UMI should usually count as **one** captured molecule, not many.

| Label | Student-friendly analogy | Function |
|---|---|---|
| **Sample index** | Building address | Separates libraries pooled on the same sequencing run |
| **Cell barcode** | Apartment number | Groups reads that came from the same droplet/cell |
| **UMI** | Serial number on one RNA molecule | Collapses PCR copies and estimates original captured molecules |

> **Key message.** UMI counts are not absolute molecule counts inside the original cell. They are counts of *captured* molecules after imperfect sampling and processing.

---

## Part 3 — FASTQ files: what the sequencer gives us

A FASTQ file is simply a text file containing sequences and a confidence score for every base. The sequencer does **not** output a Seurat object and does not know which cluster is a T cell.

### One FASTQ record has four lines

```
@A00488:123:H7J2MDSX5:1:1101:1000:1000 1:N:0:ATCACG
ACGTGCTAGCTAGCTACGTA
+
FFFFFFFFFFFFFFFFFFFF
```

| Line | Meaning |
|---|---|
| 1 | Read identifier: instrument/run/position and read information |
| 2 | Nucleotide sequence called by the sequencer |
| 3 | Separator; may repeat the identifier |
| 4 | ASCII-encoded Phred quality scores, one character per base |

### Typical 10x FASTQ filenames

```
sample01_S1_L001_R1_001.fastq.gz
sample01_S1_L001_R2_001.fastq.gz
sample01_S1_L001_I1_001.fastq.gz
sample01_S1_L001_I2_001.fastq.gz
```

| Token | Meaning |
|---|---|
| `sample01` | Sample name from the sample sheet |
| `S1` | Sample number on the sequencing run |
| `L001` | Lane 1 |
| `R1` | Read 1: usually cell barcode + UMI for 10x GEX |
| `R2` | Read 2: cDNA insert used to identify the gene |
| `I1`/`I2` | Index reads used for demultiplexing pooled libraries |
| `001` | File chunk number |

### BCL versus FASTQ

Illumina instruments initially produce base-call data (BCL). Demultiplexing uses sample indices to separate pooled libraries and writes FASTQ files. If FASTQs were delivered by a sequencing core, this step is already done. If starting from BCL, current workflows commonly use **Illumina BCL Convert**; older Cell Ranger documentation also describes `cellranger mkfastq`.

```bash
# Safe inspection commands — these do not change the data
ls -lh /path/to/fastqs
zcat sample01_S1_L001_R1_001.fastq.gz | head -n 8
zcat sample01_S1_L001_R2_001.fastq.gz | head -n 8
```

> **Key message.** Do not open a full FASTQ in a normal text editor. Files may contain hundreds of millions of records and can be tens of gigabytes.

---

## Part 4 — Cell Ranger: FASTQ to gene-by-cell matrix

Cell Ranger is the primary processing pipeline for 10x data. It connects sequence-level information to count-level information. It does **not** convert FASTQ into another FASTQ — it transforms read-level data into barcode-level and gene-level counts.

### What `cellranger count` does conceptually

| Step | Question answered |
|---|---|
| Check inputs and chemistry | Are the files and expected read structures compatible? |
| Correct and assign barcodes | Which reads belong to the same candidate droplet? |
| Align or map transcript reads | Which gene is compatible with each biological read? |
| Process UMIs | Which reads are PCR copies of the same captured molecule? |
| Call cells | Which barcodes represent cell-containing droplets rather than ambient RNA? |
| Build matrices | How many UMIs were assigned to every gene in every barcode? |
| Secondary analysis | What broad PCA/clustering/UMAP structure is present for QC? |

### Compute reality

Cell Ranger runs on 64-bit Linux. Current 10x guidance lists at least **8 CPU cores, 64 GB RAM** and substantial free disk space; real requirements grow with dataset size. Institutional HPC or server infrastructure is usually more appropriate than a laptop.

> You will **not** run Cell Ranger from scratch during the workshop — a typical dataset needs Linux, lots of RAM and hours of runtime. We demonstrate the command, then inspect a completed `outs/` folder and use a small public matrix for the R practical.

### A practical command

```bash
cellranger count \
  --id=sample01_count \
  --transcriptome=/references/refdata-gex-GRCh38-2024-A \
  --fastqs=/data/project/fastqs \
  --sample=sample01 \
  --localcores=8 \
  --localmem=64 \
  --create-bam=true
```

| Argument | Meaning | Common mistake |
|---|---|---|
| `--id` | Name of the new output directory | Using a path or an existing pipestance unintentionally |
| `--transcriptome` | Cell Ranger-compatible reference directory | Passing a FASTA file rather than a built reference |
| `--fastqs` | Directory containing FASTQs | Pointing to one FASTQ file instead of the directory |
| `--sample` | Prefix used to select this sample's FASTQs | Mismatch with the filename prefix |
| `--localcores` | Maximum local CPU cores | Requesting more than the scheduler allocation |
| `--localmem` | Maximum RAM in GB | Confusing GB with MB |
| `--create-bam` | Whether to create the position-sorted BAM | Producing a large BAM when it is not required |

### HPC submission example (adapt locally)

```bash
#!/bin/bash
#SBATCH --job-name=cr_sample01
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=cr_sample01_%j.log

module load cellranger/10.1.0   # example; use the version available locally
cellranger count --id=sample01_count \
  --transcriptome=/references/refdata-gex-GRCh38-2024-A \
  --fastqs=/data/project/fastqs --sample=sample01 \
  --localcores=$SLURM_CPUS_PER_TASK --localmem=64
```

### Outputs you should recognize

```
sample01_count/outs/
├── web_summary.html
├── metrics_summary.csv
├── filtered_feature_bc_matrix/
├── filtered_feature_bc_matrix.h5
├── raw_feature_bc_matrix/
├── raw_feature_bc_matrix.h5
├── molecule_info.h5
├── possorted_genome_bam.bam
├── possorted_genome_bam.bam.bai
├── cloupe.cloupe
└── analysis/
```

| Output | Use |
|---|---|
| `web_summary.html` | Interactive overview of run quality, cells, reads and mapping |
| `metrics_summary.csv` | Machine-readable summary for multi-sample QC tables |
| `filtered_feature_bc_matrix` | Genes × barcodes called as cells; common Seurat starting point |
| `raw_feature_bc_matrix` | Genes × all detected barcodes; useful for ambient RNA/cell-calling methods |
| `molecule_info.h5` | Molecule-level information used by downstream 10x tools |
| `BAM + BAI` | Aligned reads and index; large files used for specialized inspection |
| `cloupe.cloupe` | Interactive exploration in Loupe Browser |
| `analysis/` | Cell Ranger secondary analysis such as dimensional reductions and clusters |

### Filtered versus raw matrix

The **raw** matrix includes barcodes across the whole barcode space, most representing empty droplets with ambient RNA. The **filtered** matrix includes barcodes Cell Ranger classified as cells. Filtered does not mean biologically perfect: doublets, stressed cells and low-quality cells may remain, while some real cells can be missed.

### Web summary challenge

Given a completed `web_summary.html`, report: estimated cells, mean reads per cell, median genes per cell, fraction of reads in cells, and sequencing saturation. Then ask: does any *single* number prove the experiment is good?

| Metric | Interpretation prompt |
|---|---|
| Estimated number of cells | Is it plausible given loading and sample type? |
| Median genes per cell | Does complexity match the expected cell population and chemistry? |
| Reads mapped confidently | Are read quality, reference and library type plausible? |
| Fraction reads in cells | Could ambient RNA or damaged cells be high? |
| Sequencing saturation | Would deeper sequencing likely discover many additional UMIs? |
| Barcode-rank plot | Is there a clear transition between cell-associated and empty droplets? |

> **Key message.** QC is contextual. Compare against the expected sample biology, loading target, chemistry, sequencing design and the other samples in the experiment.

---

## Part 5 — From Cell Ranger into Seurat

Cell Ranger has built the matrix; Seurat stores the matrix with metadata and derived analyses.

### What is inside the matrix?

| Dimension | Represents |
|---|---|
| Rows | Genes or other measured features |
| Columns | Cell-associated barcodes |
| Values | Sparse UMI counts; most values are zero |

> Why is the matrix sparse? A zero may mean the gene was not expressed — but it may also mean the molecule was not captured or sequenced. This is why *absence* at single-cell level must be interpreted carefully.

### R Markdown setup

````r
---
title: "LISA scRNA-seq practical"
author: "Student name"
output: html_document
---

```{r setup, message=FALSE, warning=FALSE}
library(Seurat)
library(dplyr)
library(ggplot2)
set.seed(1234)
```
````

### 1. Import the filtered matrix

```r
data_dir <- "data/filtered_feature_bc_matrix"
counts <- Read10X(data.dir = data_dir)
dim(counts)
counts[1:5, 1:5]

pbmc <- CreateSeuratObject(
  counts = counts,
  project = "LISA_PBMC",
  min.cells = 3,
  min.features = 200
)
pbmc
```

**Alternative: read the HDF5 matrix**

```r
counts_h5 <- Read10X_h5(
  filename = "data/filtered_feature_bc_matrix.h5"
)
pbmc_h5 <- CreateSeuratObject(counts = counts_h5, project = "LISA_PBMC")
```

### 2. Calculate mitochondrial percentage

```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc[[]])
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

| Metric | Biological / technical interpretation |
|---|---|
| `nFeature_RNA` | Number of detected genes; very low may indicate empty/poor cells, very high may indicate doublets |
| `nCount_RNA` | Total UMIs; low can indicate poor capture, high can indicate large cells or doublets |
| `percent.mt` | Fraction assigned to mitochondrial genes; high values can accompany damaged or stressed cells |

### 3. Explore relationships and filter

```r
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Example thresholds for the PBMC teaching dataset only
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 2500 &
           percent.mt < 5
)
```

> **Key message.** QC thresholds are not universal. Inspect distributions by sample and cell type, and document the rationale. A threshold copied from a tutorial can remove real biology.

### 4. Normalize and identify variable features

```r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(pbmc)
head(VariableFeatures(pbmc), 10)
```

Log normalization divides each cell by its total counts, multiplies by a scale factor and log-transforms the values. Highly variable features carry informative variation for dimensional reduction.

### 5. Scale and run PCA

```r
all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
ElbowPlot(pbmc, ndims = 30)
DimPlot(pbmc, reduction = "pca")
```

> PCA creates axes that summarize coordinated variation across many genes. We keep a subset of PCs that captures useful biological structure without carrying all gene-level noise.

### 6. Build a neighbor graph, cluster and calculate UMAP

```r
dims_use <- 1:10
pbmc <- FindNeighbors(pbmc, dims = dims_use)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = dims_use)
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE)
```

> **Key message.** UMAP axes have no direct biological units. Distance can help visualize local similarity, but cluster position and island size should not be over-interpreted.

---

## Part 6 — Markers and biological annotation

### Find cluster markers

```r
markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
```

### Inspect canonical PBMC markers

```r
FeaturePlot(pbmc, features = c(
  "MS4A1",  # B cells
  "CD3D",   # T cells
  "NKG7",   # NK/cytotoxic lymphocytes
  "LYZ",    # monocytes
  "FCGR3A", # CD16 monocytes
  "PPBP"    # platelets
), ncol = 3)

DotPlot(pbmc, features = c(
  "IL7R", "CCR7", "CD3D", "CD8A", "NKG7", "GNLY",
  "MS4A1", "CD79A", "LYZ", "S100A8", "FCGR3A", "PPBP"
)) + RotatedAxis()
```

| Candidate lineage | Useful markers | Interpretation caution |
|---|---|---|
| T cells | CD3D, CD3E, TRBC1/2 | CD4 and CD8 transcripts can be sparse |
| NK cells | NKG7, GNLY, PRF1 | Cytotoxic T cells can share these genes |
| B cells | MS4A1, CD79A, CD37 | Plasma cells use different markers |
| CD14 monocytes | LYZ, S100A8, S100A9, CTSD | Inflammatory programs can blur subtypes |
| CD16 monocytes | FCGR3A, LST1, IFITM3 | May form a continuum with CD14 monocytes |
| Platelets | PPBP, PF4 | Platelet RNA can contaminate other droplets |

### Rename clusters only after checking evidence

```r
new_ids <- c(
  "Naive CD4 T", "CD14 Mono", "Memory CD4 T", "B",
  "CD8 T", "FCGR3A Mono", "NK", "DC", "Platelet"
)
names(new_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new_ids)
DimPlot(pbmc, label = TRUE, repel = TRUE) + NoLegend()
```

> **Key message.** The mapping above only works if its order matches *your* actual cluster levels. Annotation is an evidence argument, not a renaming exercise — use several positive markers, absence of incompatible markers, cluster context and sample consistency. Keep uncertain labels broad, and never paste cluster labels from another analysis without checking markers.

---

## Part 7 — Common mistakes and how to explain them

| Mistake | Why it matters | Better practice |
|---|---|---|
| Treating every cell as an independent biological replicate | Inflates significance and ignores patient-level variation | Use sample-aware/pseudobulk methods for group inference |
| Applying one QC threshold to every sample | Can selectively remove one group or cell type | Inspect distributions per sample; report decisions |
| Calling clusters solely from one marker | Markers are shared and dropout is common | Use marker panels and incompatible-lineage checks |
| Reading UMAP distance as a trajectory | UMAP is a visualization, not time | Use dedicated trajectory methods and biological evidence |
| Comparing raw cell proportions without design checks | Recovery, viability and compositional effects confound results | Model sample-level proportions and batch/covariates |
| Removing batch before checking whether batch equals biology | Can erase the true condition signal | Inspect design; integrate only when justified |
| Using filtered matrix as proof all cells are high quality | Cell calling and cell QC solve different problems | Perform independent QC, doublet and ambient RNA assessment |

### Troubleshooting terminal errors

| Symptom | Likely cause | First check |
|---|---|---|
| `cellranger: command not found` | Software is not on PATH / module not loaded | `cellranger --version`; `module avail` |
| No input FASTQs found | Wrong directory or sample prefix | `ls` FASTQ directory; compare `--sample` with filenames |
| Reference error | Incorrect/incomplete reference directory | Check reference path and Cell Ranger compatibility |
| Out of memory / killed | Scheduler allocation too small | Job log and requested memory |
| Output ID already exists | Previous pipestance has same `--id` | Inspect existing run; choose a new ID or resume appropriately |
| Very low cells / genes | Sample, library, read structure or chemistry problem | web_summary, barcode rank, mapping and chemistry |

### Optional extension topics

- **Multiple samples:** `cellranger multi`, or one count per GEM well followed by sample-aware R workflows.
- **Ambient RNA correction:** tools such as SoupX or CellBender, using raw and filtered matrices.
- **Doublet detection:** combine computational scores with expected loading and marker incompatibility.
- **Batch-aware integration:** Harmony, Seurat integration or other methods — only after examining study design. → worked example: [`integration_comparison.md`](integration_comparison.md).
- **Differential expression:** pseudobulk aggregation by sample and cell type for group-level inference.
- **TCR/BCR analysis:** use the appropriate V(D)J pipeline and connect clonotypes to transcriptomic cells.

---

## Five take-home messages

1. FASTQ is read-level data; the count matrix is cell-by-gene data.
2. Cell barcodes identify droplets; UMIs help count captured molecules without counting PCR copies repeatedly.
3. Cell Ranger performs primary processing; Seurat supports downstream analysis and interpretation.
4. QC thresholds and cell labels require biological context — not copying from a tutorial.
5. Cells are observations, but samples or patients are the biological replicates for condition-level inference.

### Exit ticket

Answer these three on a card or in the chat:

1. What is the difference between a cell barcode and a UMI?
2. Which Cell Ranger output enters Seurat?
3. Name one reason a zero count does not necessarily mean a gene was not expressed.

---

## Appendix A — Copy-ready R Markdown practical

Save this content as `LISA_scRNA_practical.Rmd`:

````r
---
title: "LISA scRNA-seq practical"
output: html_document
---

```{r setup, message=FALSE, warning=FALSE}
library(Seurat)
library(dplyr)
library(ggplot2)
set.seed(1234)
```

## Import and create object
```{r}
counts <- Read10X("data/filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts, project="LISA_PBMC", min.cells=3, min.features=200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern="^MT-")
pbmc
```

## Quality control
```{r}
VlnPlot(pbmc, c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)
pbmc <- subset(pbmc, subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)
```

## Standard workflow
```{r}
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, nfeatures=2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
ElbowPlot(pbmc, ndims=30)
pbmc <- FindNeighbors(pbmc, dims=1:10)
pbmc <- FindClusters(pbmc, resolution=.5)
pbmc <- RunUMAP(pbmc, dims=1:10)
DimPlot(pbmc, label=TRUE)
```

## Marker evidence
```{r}
markers <- FindAllMarkers(pbmc, only.pos=TRUE, min.pct=.25, logfc.threshold=.25)
markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n=5)
DotPlot(pbmc, c("CD3D","IL7R","CD8A","NKG7","GNLY","MS4A1","CD79A","LYZ","S100A8","FCGR3A","PPBP")) + RotatedAxis()
```

## Save the result
```{r}
saveRDS(pbmc, "results/LISA_pbmc_seurat.rds")
write.csv(markers, "results/LISA_cluster_markers.csv", row.names=FALSE)
```
````

---

## Appendix B — Suggested workshop folder structure

```bash
# Create a clear workshop folder structure
mkdir -p LISA_scRNA_workshop/{data,cellranger_example,results,scripts}

# Confirm downloaded FASTQs
find /data/project/fastqs -maxdepth 1 -name 'sample01*.fastq.gz' -type f

# Confirm Cell Ranger and reference
cellranger --version
ls /references/refdata-gex-GRCh38-2024-A
```

---

## Appendix C — Sources and version note

**Version note:** this guide was prepared on 3 August 2026. Cell Ranger changes over time; the current 10x documentation identified Cell Ranger v10.1 as the latest release when this guide was prepared. Confirm local software, chemistry compatibility and command-line arguments before teaching or processing study data.

- [10x Genomics — Running Cell Ranger count](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct)
- [10x Genomics — Cell Ranger outputs overview](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview)
- [10x Genomics — System requirements](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-system-requirements)
- [10x Genomics — Cell Ranger command-line arguments](https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments)
- [Seurat — PBMC 3K guided clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
- [Seurat — Read10X reference](https://satijalab.org/seurat/reference/read10x)
