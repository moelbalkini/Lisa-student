#!/usr/bin/env python3
"""
LISA scRNA-seq workshop -- group exercise test-data generator.

Creates small, SIMULATED PBMC-like single-cell datasets in 10x Genomics format
(matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz) so student groups can each
load their own dataset with Seurat's Read10X(), run the standard workflow, and
present what cell types they find.

The data is SIMULATED (not real patient data). Cell types carry canonical human
PBMC marker genes so the standard Seurat pipeline recovers them, and each dataset
contains deliberate QC problems (high-mito "damaged" cells and empty-ish low-count
barcodes) so the QC-filtering step is meaningful.

Ground truth for every dataset is written to ground_truth.json (instructor key).
"""

import os, gzip, json, struct
import numpy as np

RNG_BASE = 20260827
OUT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "datasets")

# ----------------------------------------------------------------------------
# Canonical human PBMC marker panels (signature genes per cell type)
# ----------------------------------------------------------------------------
SIGNATURES = {
    "CD4 T (naive)": ["CD3D","CD3E","CD3G","IL7R","CCR7","TCF7","LEF1","CD27","LTB","MAL","NOSIP"],
    "CD8 T":         ["CD3D","CD3E","CD8A","CD8B","GZMK","CCL5","GZMA","NKG7","KLRG1","CD7"],
    "NK":            ["NKG7","GNLY","KLRD1","KLRF1","NCAM1","GZMB","PRF1","FCGR3A","KLRB1","TYROBP"],
    "B":             ["MS4A1","CD79A","CD79B","CD37","CD19","TCL1A","BANK1","IGHM","HLA-DRA","VPREB3"],
    "CD14 Mono":     ["LYZ","S100A8","S100A9","CD14","FCN1","VCAN","S100A12","CTSD","MNDA"],
    "CD16 Mono":     ["FCGR3A","LST1","IFITM3","MS4A7","CDKN1C","LYPD2","AIF1","COTL1","SERPINA1"],
    "DC":            ["FCER1A","CLEC10A","CD1C","HLA-DQA1","FCGR2B","CST3","GPR183"],
    "Platelet":      ["PPBP","PF4","GP9","ITGA2B","TUBB1","NRGN","CAVIN2","GNG11"],
    # activation program overlaid on a monocyte identity (IFN-beta response, ifnb theme)
    "CD14 Mono (IFN-activated)": ["LYZ","S100A8","CD14","FCN1",  # still a monocyte
                                  "ISG15","IFI6","IFIT1","IFIT3","ISG20","MX1","OAS1","IRF7","STAT1"],
}

# genes broadly expressed in essentially all cells (housekeeping / ribosomal)
HOUSEKEEPING = ["ACTB","GAPDH","B2M","TMSB4X","MALAT1","FTL","FTH1","TPT1","EEF1A1",
                "RPL13","RPS27","RPS18","RPL10","RPLP1","VIM","RPS6","RPL3","JUNB"]

# mitochondrial genes -- drive percent.mt for the QC step
MT_GENES = ["MT-CO1","MT-CO2","MT-CO3","MT-ND1","MT-ND2","MT-ND4",
            "MT-ATP6","MT-CYB","MT-ND3","MT-ND5","MT-ND1","MT-ATP8"]
MT_GENES = list(dict.fromkeys(MT_GENES))  # dedupe, keep order

TOTAL_GENES = 2500  # pad with filler genes to a realistic panel size

# ----------------------------------------------------------------------------
# Per-group dataset composition.  Each entry: celltype -> n good cells.
# "theme" and "teaching_point" go into the instructor key.
# ----------------------------------------------------------------------------
GROUPS = {
    "group_A": {
        "title": "Standard PBMC",
        "teaching_point": "A clean textbook mix of the five commonest PBMC types. "
                          "Every group should be able to name all five from markers.",
        "compose": {"CD4 T (naive)":210, "CD14 Mono":150, "CD8 T":120, "B":95, "NK":80},
    },
    "group_B": {
        "title": "Rare populations hiding in the mix",
        "teaching_point": "Two rare populations (dendritic cells and platelets) sit among "
                          "the common types. Tests whether students notice small clusters "
                          "and resist over-merging them.",
        "compose": {"CD4 T (naive)":190, "CD14 Mono":150, "B":85, "NK":75, "DC":26, "Platelet":18},
    },
    "group_C": {
        "title": "Two monocyte subsets",
        "teaching_point": "CD14 and CD16 monocytes are a continuum. Tests distinguishing "
                          "closely related subsets with FCGR3A vs CD14/S100A8, and the idea "
                          "that clustering resolution changes how finely you split them.",
        "compose": {"CD14 Mono":165, "CD16 Mono":110, "CD4 T (naive)":150, "CD8 T":100, "NK":70, "B":70},
    },
    "group_D": {
        "title": "Same cell type, two activation states (IFN)",
        "teaching_point": "Monocytes appear as resting AND interferon-activated states "
                          "(ISG15/IFIT1 high). Ties to the ifnb integration story: a cluster "
                          "can split by CONDITION/STATE, not identity.",
        "compose": {"CD4 T (naive)":150, "CD14 Mono":120, "CD14 Mono (IFN-activated)":115, "NK":72, "B":72},
    },
}

# QC-failing cells added to every dataset (fractions of the good-cell count)
FRAC_DAMAGED = 0.08   # high mitochondrial %, filtered by percent.mt
FRAC_EMPTYISH = 0.05  # very low counts, filtered by nFeature_RNA > 200


def build_gene_list(rng):
    named = []
    for genes in SIGNATURES.values():
        named += genes
    named += HOUSEKEEPING + MT_GENES
    named = list(dict.fromkeys(named))  # unique, keep order
    n_filler = TOTAL_GENES - len(named)
    # filler genes get realistic-looking clone/accession-style symbols
    filler = [f"AC{rng.integers(100000,999999)}.{rng.integers(1,4)}" for _ in range(n_filler)]
    genes = named + filler
    # dedupe filler collisions
    seen, final = set(), []
    for g in genes:
        gg, k = g, 1
        while gg in seen:
            k += 1; gg = f"{g}_{k}"
        seen.add(gg); final.append(gg)
    return final


def simulate_group(name, cfg, out_dir):
    rng = np.random.default_rng(RNG_BASE + abs(hash(name)) % 100000)
    genes = build_gene_list(rng)
    gidx = {g: i for i, g in enumerate(genes)}
    G = len(genes)
    mt_rows = [gidx[g] for g in MT_GENES if g in gidx]
    hk_rows = [gidx[g] for g in HOUSEKEEPING if g in gidx]

    # per-gene background rate (most genes low; a few naturally variable)
    base = rng.lognormal(mean=np.log(0.12), sigma=0.5, size=G)
    for r in hk_rows:
        base[r] = rng.uniform(5, 11)          # housekeeping expressed everywhere
    for r in mt_rows:
        base[r] = rng.uniform(0.4, 0.8)       # modest MT in healthy cells

    cells = []          # list of (barcode, celltype, quality)
    col_counts = []     # list of dense count vectors (we keep sparse at write time)

    def emit(celltype, quality, size_factor, mt_mult=1.0, mix_with=None):
        rate = base.copy()
        # apply this cell's identity signature
        for ct in ([celltype] if mix_with is None else [celltype, mix_with]):
            sig = SIGNATURES.get(ct, [])
            for g in sig:
                if g in gidx:
                    rate[gidx[g]] = max(rate[gidx[g]], rng.uniform(9, 26))
        for r in mt_rows:
            rate[r] *= mt_mult
        lam = rate * size_factor
        counts = rng.poisson(lam).astype(np.int64)
        col_counts.append(counts)
        bc = f"{name.upper()}-{len(cells)+1:04d}-1"
        cells.append((bc, celltype if mix_with is None else f"{celltype}+{mix_with}", quality))

    # ---- good cells ----
    for ct, n in cfg["compose"].items():
        for _ in range(n):
            s = rng.lognormal(mean=0.0, sigma=0.35)     # library-size variation
            emit(ct, "good", s, mt_mult=1.0)

    n_good = len(cells)
    present_types = list(cfg["compose"].keys())

    # ---- damaged cells: high mito, moderate-low counts (filter on percent.mt) ----
    for _ in range(int(round(n_good * FRAC_DAMAGED))):
        ct = present_types[rng.integers(len(present_types))]
        s = rng.uniform(0.35, 0.6)
        emit(ct, "damaged_high_mt", s, mt_mult=rng.uniform(9, 16))

    # ---- empty-ish / low-quality: tiny counts (filter on nFeature_RNA > 200) ----
    for _ in range(int(round(n_good * FRAC_EMPTYISH))):
        ct = present_types[rng.integers(len(present_types))]
        s = rng.uniform(0.02, 0.05)
        emit(ct, "low_count", s, mt_mult=1.0)

    # assemble sparse matrix (genes x cells), MatrixMarket coordinate, integer
    N = len(cells)
    mat = np.array(col_counts).T  # G x N
    # drop all-zero genes? No -- keep panel fixed so all datasets share gene space feel.
    rows, cols = np.nonzero(mat)
    vals = mat[rows, cols]

    os.makedirs(out_dir, exist_ok=True)

    # matrix.mtx.gz  (1-based indices, header: genes cells nnz)
    with gzip.open(os.path.join(out_dir, "matrix.mtx.gz"), "wt") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write("%metadata_json: {\"software_version\": \"LISA-sim-1.0\", \"format_version\": 2}\n")
        fh.write(f"{G} {N} {len(vals)}\n")
        # column-major order is conventional; sort by column then row
        order = np.lexsort((rows, cols))
        r = rows[order] + 1
        c = cols[order] + 1
        v = vals[order]
        fh.write("\n".join(f"{a} {b} {d}" for a, b, d in zip(r, c, v)))
        fh.write("\n")

    # features.tsv.gz  (Ensembl-like id, symbol, type)  -- 3 columns like Cell Ranger
    with gzip.open(os.path.join(out_dir, "features.tsv.gz"), "wt") as fh:
        for i, g in enumerate(genes):
            ens = f"ENSGSIM{i:08d}"
            fh.write(f"{ens}\t{g}\tGene Expression\n")

    # barcodes.tsv.gz
    with gzip.open(os.path.join(out_dir, "barcodes.tsv.gz"), "wt") as fh:
        for bc, _, _ in cells:
            fh.write(bc + "\n")

    # per-cell truth table (for the checker + instructor key)
    total = mat.sum(0)
    nfeat = (mat > 0).sum(0)
    mt_counts = mat[mt_rows, :].sum(0) if mt_rows else np.zeros(N)
    pmt = np.where(total > 0, 100.0 * mt_counts / np.maximum(total, 1), 0.0)

    truth = {
        "dataset": name,
        "title": cfg["title"],
        "teaching_point": cfg["teaching_point"],
        "n_cells_total": int(N),
        "n_good_cells": int(n_good),
        "n_genes": int(G),
        "composition_good": {k: int(v) for k, v in cfg["compose"].items()},
        "qc_cells": {
            "damaged_high_mt": int(sum(1 for _, _, q in cells if q == "damaged_high_mt")),
            "low_count": int(sum(1 for _, _, q in cells if q == "low_count")),
        },
        "expected_clusters": len(cfg["compose"]),
        "markers_present": {ct: SIGNATURES[ct][:6] for ct in cfg["compose"]},
        "qc_summary_good_cells": {
            "median_nCount": float(np.median(total[:n_good])),
            "median_nFeature": float(np.median(nfeat[:n_good])),
            "median_percent_mt": round(float(np.median(pmt[:n_good])), 2),
        },
    }
    # also stash per-cell labels for the checker (not shipped to students)
    labels = [q if q != "good" else ct for (bc, ct, q) in cells]
    np.save(os.path.join(out_dir, "_truth_labels.npy"),
            np.array(labels, dtype=object), allow_pickle=True)
    return truth


def main():
    all_truth = {}
    for name, cfg in GROUPS.items():
        out_dir = os.path.join(OUT_ROOT, name)
        t = simulate_group(name, cfg, out_dir)
        all_truth[name] = t
        print(f"[{name}] {t['title']}: {t['n_cells_total']} cells "
              f"({t['n_good_cells']} good), {t['n_genes']} genes | "
              f"median nFeature {t['qc_summary_good_cells']['median_nFeature']:.0f}, "
              f"median %mt {t['qc_summary_good_cells']['median_percent_mt']}")
    with open(os.path.join(OUT_ROOT, "ground_truth.json"), "w") as fh:
        json.dump(all_truth, fh, indent=2)
    print("\nWrote ground_truth.json")


if __name__ == "__main__":
    main()
