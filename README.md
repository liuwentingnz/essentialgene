# scEssentials — essential-gene detection for scRNA-seq

Implementation of the **[scEssentials](https://github.com/huiwenzh/scEssentials)**
approach for detecting / scoring *essential genes* (stable, housekeeping-like
genes robustly expressed across cell types and platforms) in single-cell RNA-seq
data.

## What's here

```
essentialgene/
├── scEssentials/            # upstream repo (cloned, read-only reference)
│   ├── Data/
│   │   ├── scEssential_hsa.csv     # human essential genes (1969) + ES_score
│   │   └── scEssential_mmu.csv     # mouse essential genes (733)  + ES_score
│   └── scripts/                    # R scripts reproducing the manuscript analysis
├── scessentials.py          # << IMPLEMENTATION (Python, runnable) >>
├── demo_make_data.py        # builds a synthetic scRNA-seq demo dataset
├── data_demo.h5ad           # demo data (5000 genes x 300 cells)
└── demo_out.csv             # example scoring output
```

## Approach

The upstream `scEssentials` repo is primarily a *reference gene set*: the
authors pre-computed which genes are "essential" (stably and highly expressed
across >60 cell types and 10+ sequencing platforms) and published two ranked
lists with an essentiality score (`ES_score`). Their R scripts only *reproduce
the discovery* from huge atlases (Tabula Sapiens/Muris), so they aren't
practical to run on a normal dataset.

`scEssentials.py` turns that into a usable tool that works on **your** counts:

1. **Score every gene in your data** using the same logic that defines an
   essential gene:

       detection = fraction of cells where the gene is detected (count > 0)
       mean_expr = mean of log1p(CP10k-normalised) counts
       stability = 1 / (1 + CV of log expression)     # low variability = stable

       ESS = detection * logit-ranking(mean_expr) * stability

2. **Call essential genes** as the top of the ESS ranking (default: top
   quartile of non-trivial scores, or `--top-k N`).

3. **Annotate** your results against the published scEssential reference
   (`in_scEssential` / `ref_ES_score`), so you can see which of your top genes
   are independently confirmed as stable essential genes.

## Requirements

Python 3.8+ with `numpy`, `pandas`. For `.h5ad` input, also `scanpy`/`anndata`.

## Usage

```bash
# plain counts matrix, genes in rows (txt/tsv/csv)
python3 scessentials.py counts.tsv --species human --out results.csv

# 10x-style anndata object (cells x genes)
python3 scessentials.py data.h5ad --species human --top-k 100 --out results.csv

# species: 'human' (default) or 'mouse'
# top-k:   return the K most essential genes (default: automatic threshold)
```

Output columns:

| column         | meaning                                                    |
|----------------|------------------------------------------------------------|
| `Gene`         | gene name                                                  |
| `ESS`          | essentiality score (higher = more essential)               |
| `detection`    | fraction of cells with counts > 0                          |
| `mean_expr`    | mean log1p(CP10k) expression                               |
| `stability`    | 1/(1+CV) — low variability, closer to 1                     |
| `essential`    | whether the gene is called essential in this dataset       |
| `in_scEssential` | True if also in the published scEssential reference      |
| `ref_ES_score` | the author's reference score (if a scEssential)             |

## Run the demo

```bash
python3 demo_make_data.py          # builds data_demo.h5ad
python3 scessentials.py data_demo.h5ad --species human --top-k 30 --out demo_out.csv
```

The 30 planted housekeeping genes rank at the top by ESS, and 16 of them match
the published human scEssential reference.
