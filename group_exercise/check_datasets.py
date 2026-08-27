#!/usr/bin/env python3
"""
Sanity check: mimic the Seurat student pipeline in Python and confirm the
intended cell types are recoverable on each simulated dataset.

Pipeline: QC filter (nFeature>200, percent.mt<5) -> log-normalize -> top HVGs
-> scale -> PCA -> Leiden-ish (KMeans at k = #expected types) -> ARI vs truth.
Also reports, per true cell type, its top marker by mean scaled expression to
confirm canonical markers are discriminative.
"""
import os, gzip, json
import numpy as np
from scipy.io import mmread
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "datasets")

CANON = {  # a couple of canonical markers we expect to light up per lineage
    "CD4 T (naive)": ["IL7R","CCR7","CD3D"],
    "CD8 T": ["CD8A","GZMK","CD3D"],
    "NK": ["GNLY","NKG7","KLRD1"],
    "B": ["MS4A1","CD79A"],
    "CD14 Mono": ["LYZ","S100A8","CD14"],
    "CD16 Mono": ["FCGR3A","MS4A7"],
    "DC": ["FCER1A","CLEC10A"],
    "Platelet": ["PPBP","PF4"],
    "CD14 Mono (IFN-activated)": ["ISG15","IFIT1","IFI6"],
}

def read_10x(d):
    genes = []
    with gzip.open(os.path.join(d, "features.tsv.gz"), "rt") as fh:
        for line in fh:
            genes.append(line.rstrip("\n").split("\t")[1])
    bcs = [l.strip() for l in gzip.open(os.path.join(d,"barcodes.tsv.gz"),"rt")]
    m = mmread(os.path.join(d, "matrix.mtx.gz")).tocsr()  # genes x cells
    return np.array(genes), np.array(bcs), m

def main():
    truth = json.load(open(os.path.join(ROOT, "ground_truth.json")))
    ok = True
    for name in sorted(truth):
        d = os.path.join(ROOT, name)
        genes, bcs, m = read_10x(d)
        labels = np.load(os.path.join(d,"_truth_labels.npy"), allow_pickle=True)
        X = m.toarray().astype(float).T  # cells x genes
        total = X.sum(1)
        nfeat = (X > 0).sum(1)
        mt = np.array([i for i,g in enumerate(genes) if g.startswith("MT-")])
        pmt = 100.0 * X[:, mt].sum(1) / np.maximum(total,1)

        keep = (nfeat > 200) & (pmt < 5) & (total > 0)
        removed = (~keep).sum()
        # how many of the removed were actually QC-junk vs real cells?
        junk = np.array([str(l) in ("damaged_high_mt","low_count") for l in labels])
        tp = int((junk & ~keep).sum()); fn = int((junk & keep).sum())
        fp = int((~junk & ~keep).sum())

        Xf = X[keep]; lab = labels[keep]
        # log-normalize (CPM-ish to 1e4) + log1p
        sf = Xf.sum(1, keepdims=True); sf[sf==0]=1
        Xln = np.log1p(Xf / sf * 1e4)
        # HVGs: top 800 by variance
        var = Xln.var(0); hv = np.argsort(var)[::-1][:800]
        Z = Xln[:, hv]
        Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-8)
        pcs = PCA(n_components=20, random_state=0).fit_transform(Z)
        k = truth[name]["expected_clusters"]
        km = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(pcs)
        # ARI only over "good" cells (ignore any residual junk that passed QC)
        good = np.array([str(l) in truth[name]["composition_good"] for l in lab])
        ari = adjusted_rand_score(lab[good], km[good])

        # marker discrimination: for each true type, is a canonical marker enriched?
        marker_ok = []
        for ct in truth[name]["composition_good"]:
            sel = np.array([str(l)==ct for l in lab])
            best = None
            for g in CANON.get(ct, []):
                if g in set(genes):
                    gi = list(genes).index(g)
                    hi = Xln[sel, gi].mean(); lo = Xln[~sel, gi].mean()
                    if best is None or (hi-lo) > best[1]:
                        best = (g, hi-lo)
            marker_ok.append(f"{ct.split(' (')[0]}:{best[0]}(+{best[1]:.1f})" if best else f"{ct}:?")

        status = "OK" if ari >= 0.80 else "LOW"
        if ari < 0.80: ok = False
        print(f"[{name}] {status}  ARI={ari:.3f}  k={k}  "
              f"kept={keep.sum()}/{len(keep)}  QCremoved={removed} "
              f"(junk caught {tp}, junk missed {fn}, real removed {fp})")
        print("        markers:", " ".join(marker_ok))
    print("\nRESULT:", "ALL DATASETS SEPARABLE" if ok else "SOME DATASETS NEED TUNING")

if __name__ == "__main__":
    main()
