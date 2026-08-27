# Instructor Key — Group Exercise

**Do not hand this to students.** Ground truth, expected findings, and facilitation notes for the four PBMC group datasets.

---

## At a glance

| Dataset | Theme | Cell types (true) | Cells (total → after QC*) | Expected clusters |
|---|---|---|---|---|
| **group_A** | Standard PBMC | CD4 T, CD14 Mono, CD8 T, B, NK | 740 → ~628 | 5 |
| **group_B** | Rare populations | CD4 T, CD14 Mono, B, NK, **DC**, **Platelet** | 615 → ~519 | 6 |
| **group_C** | Two monocyte subsets | **CD14 Mono**, **CD16 Mono**, CD4 T, CD8 T, NK, B | 751 → ~642 | 6 |
| **group_D** | Activation states (IFN) | CD4 T, **CD14 Mono (resting)**, **CD14 Mono (IFN)**, NK, B | 597 → ~507 | 5 |

\* After `nFeature_RNA > 200 & percent.mt < 5`. All datasets have ~2500 genes; good cells sit at median ~340 genes, ~650 UMIs, ~1% mito.

Each dataset was seeded with **damaged high-mito cells** (~8%) and **empty-ish low-count barcodes** (~5%). A correct QC filter removes essentially all of them; students who skip QC will see a messy UMAP with a low-quality smear.

---

## group_A — Standard PBMC *(easiest — good warm-up group)*

Ground-truth composition (good cells): CD4 T (naive) 210 · CD14 Mono 150 · CD8 T 120 · B 95 · NK 80.

- **Expected result:** five clean, well-separated clusters, one per lineage.
- **Marker evidence:** CD4 T → CD3D/IL7R/CCR7; CD8 T → CD3D/CD8A/GZMK; NK → NKG7/GNLY/KLRD1 (CD3D-negative); B → MS4A1/CD79A; CD14 Mono → LYZ/S100A8/CD14.
- **Common trap:** calling NK and CD8 T the same cluster. Point them to **CD3D** — present in CD8 T, absent in NK.

## group_B — Rare populations

Composition: CD4 T 190 · CD14 Mono 150 · B 85 · NK 75 · **DC 26** · **Platelet 18**.

- **Expected result:** the four common clusters plus **two tiny clusters**.
- **Marker evidence:** DC → FCER1A/CLEC10A/CD1C; Platelet → PPBP/PF4 (extremely distinct — usually the cleanest marker in the whole exercise).
- **Teaching point:** small ≠ unimportant. At low resolution the platelet/DC clusters may sit as a handful of cells; the marker signal is unambiguous. Watch for groups that over-merge and report only 4 types — ask them to look for PPBP.

## group_C — Two monocyte subsets

Composition: **CD14 Mono 165** · **CD16 Mono 110** · CD4 T 150 · CD8 T 100 · NK 70 · B 70.

- **Expected result:** six clusters, but the two monocyte clusters sit next to each other on the UMAP.
- **Marker evidence:** CD14 Mono → CD14/S100A8/LYZ; CD16 Mono → FCGR3A/MS4A7 with **low CD14**.
- **Teaching point:** resolution matters. At `resolution = 0.3` the two monocyte subsets may **merge**; at `0.5–0.8` they split. This is the intended lesson — clustering granularity is a choice, justified by markers. (Note FCGR3A also marks NK, so require MS4A7/LST1 to call CD16 monocytes, not FCGR3A alone.)

## group_D — Activation states (IFN) *(ties directly to the ifnb practical)*

Composition: CD4 T 150 · **CD14 Mono resting 120** · **CD14 Mono IFN-activated 115** · NK 72 · B 72.

- **Expected result:** five clusters — but **two of them are both monocytes**.
- **Marker evidence:** both monocyte clusters share LYZ/S100A8/CD14; the activated one *additionally* has high **ISG15/IFIT1/IFI6/MX1**.
- **Teaching point:** a cluster can split by **state/condition**, not identity. This is the whole motivation for the `ifnb` integration story (Parts 1–3): if these were two samples (CTRL vs STIM), you would integrate so the two monocyte states align, then recover the interferon response *within* the monocyte cluster by differential expression. Great bridge into the integration practical.
- **Extension:** `FindMarkers` between the two monocyte clusters returns the interferon-response gene set — the same biology as the real Kang et al. data.

---

## Running the groups

- **Group size:** 3–4 students per dataset works well. With more than four groups, assign the same dataset to two groups — comparing two independent analyses of the *same* data is itself a good discussion (did they filter the same, find the same clusters?).
- **Suggested pairing to outcomes:** give group_A to a group that needs an easier start; save group_D for a confident group, as it previews integration.
- **Environment:** groups need the Part-1 packages installed (Seurat v5, dplyr, ggplot2). No internet is required for this exercise — the data is local, so it works offline in the room.
- **Deliverable:** the 5-slide format in the student README. Budget ~5 min per group + 2 min questions.

## What "good" looks like in a presentation

- QC justified by looking at *this* dataset's plots, not a copied threshold.
- Every cell-type call backed by **≥2 markers** and an absent incompatible marker.
- Honest uncertainty (e.g. CD8 T vs NK, or CD14 vs CD16 boundary).
- The twist identified and explained.

## Regenerating or editing the datasets

`generate_datasets.py` rebuilds all four folders deterministically and rewrites `datasets/ground_truth.json`. To change composition, edit the `GROUPS` dictionary (cell counts / types) and re-run. Requires `numpy` and `scipy`. To verify separability after editing, `check_datasets.py` reruns a QC→PCA→clustering pass and reports the adjusted Rand index against ground truth (expect ~1.0).
