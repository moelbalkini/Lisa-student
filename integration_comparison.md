# scRNA-seq Practical Series — With vs. Without Integration

A hands-on, three-part companion to the [LISA scRNA-seq workshop](README.md). You download a **public** two-condition PBMC dataset (`ifnb`: control vs. IFN-β–stimulated), analyse it **without** and **with** batch integration, and measure the stimulation response the statistically correct way.

Work through the parts in order — each saves a checkpoint the next one loads.

| Part | File | What you do |
|---|---|---|
| **1 — Setup & installation** | [`practical_part1_setup.md`](practical_part1_setup.md) | Install every package, verify Seurat v5, download the public data, set up the project. |
| **2 — Loading & downstream analysis** | [`practical_part2_analysis.md`](practical_part2_analysis.md) | QC → normalize → PCA → cluster → UMAP, run it **without** and **with** integration (CCA + Harmony), and compare side by side. |
| **3 — Differential expression & annotation** | [`practical_part3_differential_expression.md`](practical_part3_differential_expression.md) | Name the clusters, find the stimulation response **within** cell types, and do it properly with pseudobulk. |

---

## The dataset

**`ifnb`** (Kang et al., 2018): ~14,000 human PBMCs in two conditions — resting **control (CTRL)** and **interferon-β–stimulated (STIM)**. It is public and ships through the `SeuratData` package, so no manual download or Cell Ranger run is needed. Because the two conditions were processed separately and stimulation changes many genes, cells cluster **by condition** unless you integrate — a clear, reproducible example of the effect integration addresses.

> **Swap in your own data.** Part 1 shows where to replace the `ifnb` load with your own merged samples; the rest of the workflow is identical as long as a metadata column labels the batch/sample.

---

## The big picture

1. **Setup (Part 1)** — get a clean, loaded Seurat object.
2. **Analyse (Part 2)** — the standard workflow, run twice; *see* the batch effect and correct it. **Integration aligns cell types across samples so clustering reflects biology, not batch.**
3. **Interpret (Part 3)** — recover the treatment signal *within* aligned cell types. **Cells are observations; samples/patients are the replicates** — so condition-level DE is done with pseudobulk.

Everything uses **Seurat v5** syntax. Start with **[Part 1 →](practical_part1_setup.md)**.

---

### Sources
- [Seurat v5 — Integrative analysis](https://satijalab.org/seurat/articles/integration_introduction)
- [Seurat v5 — Differential expression](https://satijalab.org/seurat/articles/de_vignette)
- [SeuratData](https://github.com/satijalab/seurat-data) · [Harmony](https://github.com/immunogenomics/harmony)
- Kang et al. (2018), Nat. Biotechnol. — the `ifnb` dataset.
