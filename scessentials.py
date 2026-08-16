#!/usr/bin/env python3
"""
scEssentials -- implementation for detecting / scoring essential genes
in single-cell RNA-seq data.

Background
----------
scEssentials (https://github.com/huiwenzh/scEssentials) is a reference gene set of
"essential genes" that are robustly, stably and highly expressed across >60 cell
types and 10+ sequencing platforms.  The authors publish two precomputed lists:

    Data/scEssential_hsa.csv   (human, 1969 genes, with ES_score)
    Data/scEssential_mmu.csv   (mouse, 733  genes, with ES_score)

The published R scripts reproduce the *discovery* of those lists from huge
atlases (Tabula Sapiens / Muris, benchmarking datasets).  For a normal user the
practical, implementable workflow is:

  (A) Score every gene in your own scRNA-seq data by the SAME logic the authors
      used to justify scEssentials --- a gene is "essential" if it is
      *detected widely* (high % of cells) and *highly expressed* and *stable*
      (low variability).  This gives a per-gene essentiality score you can rank.

  (B) Intersect / annotate your genes with the published scEssential reference
      lists so you can flag the known-stable essential genes present in your data.

This module implements both.  It needs only numpy/scanpy/anndata (or a plain
genes x cells count matrix).

ESS score (this implementation)
-------------------------------
    detection  = fraction of cells with count > 0                  >=0..1
    mean_expr  = mean of log1p(cpm-normalised) counts              > 0
    stability  = 1 / (1 + cv_of_log_expression)                    in (0..1]

    ESS = detection * z-score_rank(mean_expr) * stability

A gene is called "essential in this dataset" when it is in the top-k by ESS
across the data OR it is a member of the published scEssential reference with a
non-trivial detection rate.  Thresholds are tunable.
"""

import os
import sys

import numpy as np
import pandas as pd

__all__ = [
    "load_reference",
    "essentiality_score",
    "detect_essential_genes",
    "annotate_with_reference",
    "save_report",
]

DEFAULT_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "scEssentials", "Data")


# --------------------------------------------------------------------------- #
# (A) Reference lists
# --------------------------------------------------------------------------- #
def load_reference(species="human", data_dir=DEFAULT_DATA_DIR):
    """Load the published scEssential gene list.

    species : 'human' (hsa) or 'mouse' (mmu)
    Returns a DataFrame with columns ['Gene','ES_score'] (ES_score is the
    author's original essentiality score; higher = more essential).
    """
    fname = {"human": "scEssential_hsa.csv", "mouse": "scEssential_mmu.csv"}
    if species not in fname:
        raise ValueError("species must be 'human' or 'mouse'")
    path = os.path.join(data_dir, fname[species])
    if not os.path.exists(path):
        if species == "human":
            alt = os.path.join(data_dir, "scEssential_hsa.csv")
            path = alt if os.path.exists(alt) else path
        raise FileNotFoundError(f"Reference file not found: {path}\n"
                                "Clone it first:  git clone "
                                "https://github.com/huiwenzh/scEssentials.git")
    df = pd.read_csv(path)
    # normalise column names / strip surrounding quotes
    df.columns = [c.strip().strip('"').strip() for c in df.columns]
    map_ = {c.lower(): c for c in df.columns}
    gene_col = map_.get("gene", df.columns[-1])
    score_col = map_.get("es_score", "ES_score")
    out = pd.DataFrame({
        "Gene": df[gene_col].astype(str).str.strip(),
        "ES_score": pd.to_numeric(pd.to_numeric(df[score_col].astype(str)
                                  .str.replace('"', ''), errors="coerce"),
                                  errors="coerce"),
    })
    out = out.dropna(subset=["Gene"]).drop_duplicates("Gene")
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# (B) Data-driven essentiality scoring
# --------------------------------------------------------------------------- #
def _as_count_matrix(X, feature_names, cell_names=None):
    """Accept either an anndata.AnnData, a genes x cells numpy matrix, or a
    pd.DataFrame (genes x cells) and return (counts, feature_names)."""
    # AnnData -- scanpy convention is cells x genes; return genes x cells.
    if X.__class__.__name__ == "AnnData":
        counts = X.raw.X if X.raw is not None else X.X
        if hasattr(counts, "toarray"):
            counts = counts.toarray()
        counts = np.asarray(counts, dtype=float)
        if counts.ndim != 2:
            raise ValueError("AnnData .X must be 2-D")
        feats = np.asarray(X.raw.var_names) if X.raw is not None \
            else np.asarray(X.var_names)
        if len(feats) == counts.shape[1]:
            counts = counts.T          # cells x genes -> genes x cells
        return counts, feats

    # pandas DataFrame (genes x cells)
    if isinstance(X, pd.DataFrame):
        return X.values.astype(float), np.asarray(X.index)

    # ndarray / sparse: assume genes x cells; feature_names required
    if isinstance(X, np.ndarray) or hasattr(X, "toarray"):
        if hasattr(X, "toarray"):
            X = X.toarray()
        arr = np.asarray(X, dtype=float)
        if feature_names is None:
            feature_names = np.arange(arr.shape[0])
        return arr, np.asarray(feature_names)

    raise TypeError("X must be AnnData, DataFrame, ndarray or scipy sparse "
                    "(genes x cells)")


def essentiality_score(X, feature_names=None, min_detection=0.01,
                       cpm=True):
    """Compute a per-gene essentiality score from a single-cell count matrix.

    Parameters
    ----------
    X : AnnData | DataFrame(genes x cells) | ndarray(genes x cells)
    feature_names : feature (gene) names aligned to rows
    min_detection  : exclude genes detected in fewer than this fraction of cells
    cpm            : normalise counts to counts-per-10k (CP10k, scRNA convention)
                     before computing mean expression

    Returns
    -------
    DataFrame with columns:
        Gene, detection, mean_expr, cv, stability, ESS, essential(bool)
    sorted by ESS descending.
    """
    counts, feats = _as_count_matrix(X, feature_names)
    counts = np.asarray(counts, dtype=float)
    n_cells = counts.shape[1]

    if cpm:
        total = counts.sum(axis=1, keepdims=True)
        total[total == 0] = np.nan
        norm = counts / total * 1e4
        norm = np.nan_to_num(norm)
    else:
        norm = counts

    logx = np.log1p(norm)

    detection = (counts > 0).sum(axis=1) / n_cells          # fraction of cells
    mean_expr = logx.mean(axis=1)                            # mean log expr
    std_expr = logx.std(axis=1, ddof=0)
    cv = np.where(mean_expr > 0, std_expr / np.maximum(mean_expr, 1e-9), np.inf)
    stability = 1.0 / (1.0 + cv)                             # 0..1

    # guard against genes with zero variance
    stability[np.isnan(cv)] = 0.0

    # z-score(mean_expr) -> 0..1 rank via logistic; keep interpretable scale
    m = mean_expr[~np.isnan(mean_expr)]
    mu, sd = m.mean(), (m.std() + 1e-12)
    mean_rank = 1.0 / (1.0 + np.exp(-(mean_expr - mu) / sd))  # 0..1

    idx = np.arange(len(feats))
    keep = (detection >= min_detection) & np.isfinite(mean_expr)

    df = pd.DataFrame({
        "Gene": feats,
        "detection": detection,
        "mean_expr": mean_expr,
        "cv": cv,
        "stability": stability,
        "mean_rank": mean_rank,
    })

    score = np.zeros(len(df), dtype=float)
    score[keep] = (df.loc[keep, "detection"]
                   * df.loc[keep, "mean_rank"]
                   * df.loc[keep, "stability"])
    df["ESS"] = score

    # default call: gene is "essential in this data" if it meets detection
    # threshold AND is in the top quartile of non-trivial scores.
    pos = df.loc[keep, "ESS"]
    thr = pos.quantile(0.75) if len(pos) else 0.0
    df["essential"] = (df["ESS"] >= thr) & keep

    df = df.sort_values("ESS", ascending=False).reset_index(drop=True)
    return df


def detect_essential_genes(X, feature_names=None, top_k=None, min_detection=0.01,
                           cpm=True):
    """Return only the genes called essential in the dataset, ranked by ESS.

    X can be a raw count matrix (AnnData/DataFrame/ndarray) OR a DataFrame
    already produced by essentiality_score().

    top_k : if given, return top-k by ESS instead of the default threshold.
    """
    if isinstance(X, pd.DataFrame) and "ESS" in X.columns:
        df = X
    else:
        df = essentiality_score(X, feature_names=feature_names,
                                min_detection=min_detection, cpm=cpm)
    if top_k is not None:
        pos = df[df["detection"] >= min_detection]
        return pos.head(top_k).reset_index(drop=True)
    return df[df["essential"]].reset_index(drop=True)


# --------------------------------------------------------------------------- #
# (C) Annotation with published scEssential reference
# --------------------------------------------------------------------------- #
def annotate_with_reference(scores, species="human", data_dir=DEFAULT_DATA_DIR):
    """Join a scored gene table (must have a 'Gene' column) with the published
    scEssential reference and add columns:
        ref_ES_score  : author's reference score (NaN if not a scEssential)
        in_scEssential: bool flag
    """
    ref = load_reference(species=species, data_dir=data_dir)
    out = scores.merge(
        ref.rename(columns={"ES_score": "ref_ES_score"}),
        left_on="Gene", right_on="Gene", how="left",
    )
    out["in_scEssential"] = out["ref_ES_score"].notna()
    return out


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def save_report(df, out_path, species="human"):
    """Write the scored table to CSV (and mark which genes are scEssentials)."""
    cols = [c for c in ["Gene", "ESS", "detection", "mean_expr", "stability",
                        "essential", "in_scEssential", "ref_ES_score"]
            if c in df.columns]
    df[cols].to_csv(out_path, index=False)
    return out_path


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    import argparse

    ap = argparse.ArgumentParser(
        prog="scessentials",
        description="Detect / score essential genes in scRNA-seq data "
                    "(implementation of huiwenzh/scEssentials).",
    )
    ap.add_argument("input", help="Count matrix: .h5ad / .txt/.tsv/.csv "
                                  "(genes x cells) / .mtx(+genes)")
    ap.add_argument("--species", choices=["human", "mouse"], default="human")
    ap.add_argument("--genes-file", help="gene names file for .mtx input "
                                         "(one per line)")
    ap.add_argument("--top-k", type=int, default=None,
                    help="return top-K essentials instead of default threshold")
    ap.add_argument("--min-detection", type=float, default=0.01)
    ap.add_argument("--out", default="essential_genes.csv")
    args = ap.parse_args(argv)

    f = args.input
    if f.endswith(".h5ad"):
        import anndata
        ad = anndata.read_h5ad(f)
        X, feats = _as_count_matrix(ad.X if ad.raw is None else ad.raw,
                                    np.asarray(ad.raw.var_names) if ad.raw else
                                    np.asarray(ad.var_names))
        # _as_count_matrix expects AnnData; handle orientation
        scores = essentiality_score(ad, min_detection=args.min_detection)
    else:
        sep = "\t" if f.endswith((".txt", ".tsv")) else ","
        df = pd.read_csv(f, sep=sep, index_col=0)
        scores = essentiality_score(df, min_detection=args.min_detection)

    if args.top_k:
        det = detect_essential_genes(scores, top_k=args.top_k,
                                     min_detection=args.min_detection)
    else:
        det = scores

    det = annotate_with_reference(det, species=args.species)
    save_report(det, args.out, species=args.species)
    ref_n = int(det["in_scEssential"].sum()) if "in_scEssential" in det else 0
    print(f"[scessentials] scored {len(det)} genes; "
          f"{int(det['essential'].sum()) if 'essential' in det else ''} "
          f"called essential; {ref_n} are published scEssentials.")
    print(f"[scessentials] report written to {args.out}")


if __name__ == "__main__":
    sys.exit(main())
