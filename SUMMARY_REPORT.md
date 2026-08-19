---
title: "Deliverables Summary Report"
subtitle: "Essential-Gene Detection: Implementation, Research Survey, and Review Manuscript"
author: "liuwenting"
date: "2026-08-18"
---

# Summary Report — Essential-Gene Detection Workstream

This report consolidates everything produced in this session into one document:
the `scEssentials` implementation, the research/model/benchmark survey, the
comparison of methods, and the review manuscript. All deliverables live under
`/data/liuwenting/tools/essentialgene/` and are pushed to
`github.com/liuwentingnz/essentialgene` (branch `main`).

---

## 1. Scope

The work had two threads that converged:

1. **A working tool.** Implement an "essential-gene" detector for single-cell
   RNA-seq data, based on the published scEssentials reference.
2. **A research review.** Survey what data and models the field uses, compare
   them, identify the best current approaches, and propose future directions —
   delivered as a review-style manuscript plus supporting reports.

---

## 2. Deliverables (this session)

| File | Contents |
|------|----------|
| `scessentials.py` | Runnable Python implementation: scores essential genes in scRNA-seq data and flags them against the published human/mouse scEssential reference |
| `demo_make_data.py`, `data_demo.h5ad`, `demo_out.csv` | Synthetic demo dataset and its scoring output (validates the pipeline end-to-end) |
| `README.md` | Usage and design documentation |
| `research_ideas_report.md` | Initial English research-ideas report (models, benchmarks, v2 roadmap) |
| `research_report_zh.md` | Chinese research report: required data, best-model survey, comparison, recommended directions |
| `manuscript_review.md` / `.docx` | **Review-style manuscript** (the main deliverable): methods, benchmarks, comparison, future directions, references |
| `scEssentials/` | Upstream reference clone (gene lists + R discovery scripts) |
| `.gitignore` | Excludes upstream clone + generated data from the repo |

---

## 3. The Implementation

**scEssentials** (upstream) is a *reference gene set*, not a runnable detector: the
authors pre-computed which genes are stably/highly expressed across >60 cell types
and 10+ platforms, and published ranked lists with an essentiality score:

- `scEssential_hsa.csv` — 1969 human essential genes
- `scEssential_mmu.csv` — 733 mouse essential genes
- R scripts that only *reproduce the discovery* from huge atlases (impractical to run)

`scEssentials.py` makes this usable on any dataset:

```
ESS = detection × expression-rank × stability
   detection  = fraction of cells where the gene is detected (count > 0)
   mean_expr  = mean log1p(CP10k-normalised) expression
   stability  = 1 / (1 + CV of log expression)     # low variability = stable
```

It also annotates results against the published reference (`in_scEssential`,
`ref_ES_score`). **Verified end-to-end**: on a synthetic demo, the 30 planted
housekeeping genes ranked top by ESS, and 16 matched the published human reference.
Both human and mouse reference lists load correctly.

Command-line usage:
```
python3 scessentials.py data.h5ad --species human --top-k 100 --out results.csv
```

---

## 4. What Data Are Needed to Identify Essential Genes

**Features (model input):**
- scRNA-seq expression — Tabula Sapiens, Tabula Muris Senis, 10x PBMC
- Cross-platform expression — Ziegenhain 2017 (6 platforms), Chen 2020 (7 platforms)
- Chromatin accessibility (scATAC), PPI/metabolic networks (STRING, BioGRID)
- Sequence/genomic features + DNA language models (for unexpressed genes)
- Multi-omics cancer covariates — CCLE

**Labels (ground truth):**
- CRISPR essentiality — **DepMap CERES** (~1100 lines × 18k genes), **Project Score** (Sanger)
- Core/true-essential sets — **CoRe**, **Hart/BEAN**
- Context-specific — **ADAM/ADAM2**
- RNAi cross-validation — DRIVE, Achilles, GenomeRNAi
- Microbial — Keio (*E. coli*), DEG, OGEE
- Housekeeping catalogs — Hounkpe, Eisenberg–Levanon, PanglaoDB

**Key insight:** scEssentials currently uses only expression features (no labels);
the most direct route to a better model is to supervise the expression score with
functional labels (DepMap/Project Score).

---

## 5. State of the Art (Surveyed)

### 5.1 Expression-robustness
**scEssentials**, **scRoGG** (co-expression). Expression-only, interpretable.

### 5.2 Functional essentiality (CRISPR gold standard)
**DepMap/CERES**, **Project Score**, **CoRe**, **ADAM/ADAM2**, **MAGeCK**,
**Chronos**, **ACE**.

### 5.3 Deep-learning predictors
DeepHE/DeepSF, **EPGAT**, **FluxGAT**, **DeepHEN**, **LeAP**, topology-ML,
multiplex-network DL, explainable target DL, **EssentialGIN**, **DeepHEM**,
**DeEPsnap**, **DeepVul**, **XGEP**, DeeplyEssential, spiking NN, **MetaGEM**,
**dnaHNet**.

### 5.4 Generative / scCRISPR (frontier)
**DEGAS** (2023), **scLAMBDA** (2024), **PerturbNet** (2025), scCRISPR tiling
(2021).

---

## 6. Comparison (High Level)

| Dimension | Expression-robust (scEssentials) | Functional (DepMap/CRISPR) | Deep learning |
|-----------|--------------------------------|----------------------------|---------------|
| Data need | scRNA-seq only | pooled knockout screen | features + labels |
| Interpretability | High | High | Low–medium |
| Genomic coverage | expressed only | genome-wide | feature-dependent |
| Cross-species | recompute | hard | transferable |
| Runs on new data | Yes | No | needs trained model |
| Meaning | housekeeping/reference | knockout lethality | fits functional labels |

Main trade-offs: expression-robust is cheap but not "functional"; functional is
authoritative but expensive and cancer-biased; deep learning is powerful but less
interpretable and label-dependent. Graph/topology and generative scCRISPR methods
add network semantics and single-cell heterogeneity respectively.

---

## 7. Recommended Future Directions (ranked)

1. **Learn the expression score against functional labels** (GBM → DNN/GNN on
   DepMap/Project Score). *Highest biological value-to-effort.*
2. **Context-/cell-type-conditional essentiality** (core vs. context; CCLE
   covariates). *Biggest novelty.*
3. **Network centrality/topology as an essentiality signal** (EPGAT, FluxGAT,
   topology-ML); scores unexpressed genes.
4. **Sequence / DNA-LLM features** for unexpressed genes and cross-species
   transfer.
5. **Generative scCRISPR alignment** (DEGAS/scLAMBDA-style) where scCRISPR data
   exist.
6. **Co-essentiality modules** with GO/KEGG enrichment instead of flat lists.
7. **Built-in benchmark suite** (platform-consistency, DepMap AUC, Hart/BEAN
   overlap, top-k enrichment).

### Roadmap (v2.0)
Benchmark suite → supervised DepMap score → context-conditional → network/module
delivery → DNA-LLM / cross-species (stretch).

---

## 8. Key Sources & Repositories (verified 2026-08-18)

**Repositories:** huiwenzh/scEssentials, huiwenzh/scRoGG, DepMap-Analytics/CoRe,
DepMap-Analytics/ADAM2, francescojm/ADAM, broadinstitute/depmap-portal,
liulab-dfci/MAGeCK, WubingZhang/MAGeCKFlute, xzhang2016/DeepHE,
QIANJINYDX/EssTFNet, ucrbioinfo/DeeplyEssential, alaa27/DeepVul,
BioDataLearning/XGEP, Bitbol-Lab/ProteomeLM.

**Datasets:** DepMap, Project Score, CCLE, DRIVE, Keio, Hart/BEAN, Perturb-seq,
GenomeRNAi, L1000/CMAP.

**Note on verification:** arXiv IDs for EPGAT/FluxGAT/DeepHEN/LeAP etc. are solid.
DEGAS and scEVE/CHNet lack confirmed canonical public repositories; DEGAS is
journal-only (Nat Commun 2023) and its exact DOI could not be auto-verified from
title search — flagged for confirmation before submission.

---

## 9. Outstanding Items / Blockers

- **NASH dataset access** (requested earlier): `/media/nfs/nfs02/huyang/NASH/
  Monkey_NASH_model/raw_data` on `192.168.66.231` is an empty root-owned skeleton —
  the NFS share is not mounted and SMB/NFS ports time out from this host. Not
  resolved.
- **DEGAS / scEVE / CHNet DOIs** to confirm before manuscript submission.
- Optionally: add a comparison figure/table to the manuscript; add a Chinese
  version of the expanded manuscript.
