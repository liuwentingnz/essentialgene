---
title: "Computational Identification of Essential Genes in Single-Cell and CRISPR Data: A Review of Methods, Benchmarks, and Future Directions"
author:
  - "L. Wenting"
date: "2026-08-16"
abstract: |
  Essential genes — genes whose disruption is incompatible with cell viability — are
  central to understanding cellular fitness, cancer dependencies, and housekeeping
  functions. Recent advances in single-cell RNA sequencing (scRNA-seq) and pooled
  CRISPR-Cas9 screening have produced two largely separate literatures: (i)
  expression-robustness methods such as scEssentials that nominate "essential"
  genes from stable, broad, high expression across cell types and platforms, and
  (ii) functional-essentiality methods (DepMap, CoRe, ADAM) that measure the
  phenotypic consequence of gene knockout. This review formalises the data
  requirements of each approach, surveys and benchmarks the current state of the
  art — including deep-learning predictors (DeepHE, EssTFNet, DeepVul, XGEP) — and
  critically compares the two paradigms. We argue that the most promising direction
  is a supervised, context-aware synthesis that uses functional CRISPR labels
  (DepMap) to learn and reweight the expression-derived essentiality score, extended
  by sequence / DNA-language-model features for unexpressed genes and cross-species
  transfer, and delivered as co-essentiality modules rather than flat gene lists.
keywords: [essential genes, scRNA-seq, CRISPR, DepMap, deep learning, gene dependency, review]
---

# 1. Introduction

Genes whose loss is lethal or strongly impairs cell fitness are termed *essential*
genes. They underpin a cell's baseline viability, serve as reference points for
normalisation in omics assays, and, in cancer, define the concept of "dependencies"
that are promising therapeutic targets. Historically, essentiality was defined by
conservation and single-gene deletion studies in microbes (DEG, OGEE). The
advent of single-cell transcriptomics and pooled CRISPR screens has transformed how
essentiality is measured and predicted in mammalian and complex systems.

Two communities have evolved somewhat independently:

- **Expression-robustness approach.** Methods such as **scEssentials** identify
  genes that are *stably, broadly, and highly expressed* across cell types and
  sequencing platforms. The implicit assumption is that a gene transcribed at high,
  low-variance levels everywhere is essential for general cellular function (a
  refined "housekeeping" definition tailored to single-cell data).
- **Functional-essentiality approach.** Pooled **CRISPR-Cas9** screens directly
  measure the effect of gene knockout on cell fitness, yielding per-gene
  "essentiality" or "dependency" scores (e.g. DepMap CERES). This is the
  operational gold standard for "essential".

The two definitions are related but **not identical**: expression-robust genes are
often functionally essential, but the reverse is not guaranteed, and highly
expressed genes can be dispensable. The central open question — and the subject of
this review — is how to reconcile and improve these paradigms.

This review makes three contributions:

1. **Data requirements** — what inputs each approach needs (Section 2);
2. **State of the art and comparison** — a survey and head-to-head of current
   expression-based, functional, and deep-learning methods, with strengths and
   weaknesses (Sections 3–5);
3. **Future directions** — a rank-ordered agenda for the next generation of
   essential-gene tools (Section 6).

# 2. Data Requirements for Essential-Gene Identification

We separate inputs into *features* (used to score or predict) and *labels*
(the ground truth used to supervise or validate).

## 2.1 Feature data (model input)

| Modality | Representative resources | Purpose |
|----------|--------------------------|---------|
| scRNA-seq expression | Tabula Sapiens, Tabula Muris Senis, 10x PBMC | detection rate, mean expression, CV, cell-type breadth |
| Cross-platform expression | Ziegenhain 2017 (6 platforms), Chen 2020 (7 platforms) | cross-platform robustness |
| Chromatin accessibility | scATAC (Cusanovich et al.) | link to open chromatin |
| Protein–protein interaction | STRING, BioGRID | network centrality as importance feature |
| Sequence / genomic | conservation, allele annotation, raw DNA | cover unexpressed genes |
| Protein language model | ProteomeLM, ESM | transferable, cross-species embeddings |

## 2.2 Label data (ground truth)

| Label type | Representative resources | Purpose |
|------------|--------------------------|---------|
| CRISPR essentiality | **DepMap CERES** (~1000 lines × 18k genes) | functional gold-standard score |
| Core-fitness sets | DepMap-Analytics/CoRe | core/common-essential positives |
| Context-specific sets | ADAM / ADAM2 | core vs. context distinction |
| Essential-gene catalogs | DEG, OGEE | cross-species validation |
| Housekeeping catalogs | Hounkpe, Eisenberg–Levanon, PanglaoDB | overlap / discordance analysis |

> **Key point.** scEssentials currently uses only expression features (no labels).
> The most direct route to a *better* model is to introduce functional labels
> (DepMap) as a supervision signal to *learn* the expression-derived score.

# 3. State of the Art: Method Survey

## 3.1 Expression-robustness methods

**scEssentials** (huiwenzh/scEssentials). Defines a reference set of genes stably
and highly expressed across >60 cell types and 10+ sequencing platforms. Score is
`ESS ≈ detection · expression-rank · stability`. Directly runnable on new data (the
accompanying implementation `scessentials.py`). **scRoGG** (huiwenzh/scRoGG) reuses
this set as an anchor for robust co-expression networks.

## 3.2 Functional-essentiality methods

**DepMap / CERES.** Genome-wide CRISPR-Cas9 dependency scores across ~1000 cancer
cell lines. The reference label set. **CoRe** (DepMap-Analytics/CoRe) integrates
multiple screens to define core-fitness and common-essential genes. **ADAM /
ADAM2** (francescojm/ADAM) classify genes as core-fitness (essential everywhere)
versus context-specific (essential in a subset of lines). **MAGeCK**
(liulab-dfci/MAGeCK) is the upstream tool that turns raw screen FASTQs into
essentiality calls.

## 3.3 Deep-learning predictors

| Tool | Repository | Approach |
|------|-----------|----------|
| DeepHE / DeepSF | xzhang2016/DeepHE | DNN over sequence/expression/network features |
| EssTFNet | QIANJINYDX/EssTFNet | temporal-frequency DNN + DNA language model |
| DeeplyEssential | ucrbioinfo/DeeplyEssential | microbial essentiality DNN |
| DeepVul | alaa27/DeepVul | multi-task transformer: essentiality + drug response |
| XGEP | BioDataLearning/XGEP | expression-based essentials + lncRNAs |

**Methods not yet located in a canonical public repository.** DEGAS (deep
generative model of essentiality from single-cell CRISPR screens), scEVE, and
CHNet appear in the literature but their canonical GitHub repositories could not
be confirmed in this survey; they should be sourced from their papers directly.

# 4. Comparison and Strengths / Weaknesses

## 4.1 Cross-cutting comparison

| Dimension | Expression-robust (scEssentials) | Functional (DepMap/CRISPR) | Deep learning |
|-----------|--------------------------------|----------------------------|---------------|
| Data need | scRNA-seq only | pooled knockout screen | features + labels |
| Interpretability | High (transparent formula) | High (experimental) | Low–medium |
| Genomic coverage | expressed genes only | genome-wide | feature-dependent |
| Cross-species | recompute needed | hard (per-species screens) | transferable |
| Runs on new data | Yes | No (needs screen data) | needs trained model |
| Biological meaning | housekeeping / reference | direct knock-out lethality | fits functional labels |

## 4.2 Detailed strengths and weaknesses

**scEssentials (baseline).**
*Strengths:* lightweight, expression-only, interpretable, runs on any new
scRNA-seq data.
*Weaknesses:* (1) "expression-robust" ≠ "functionally essential"; (2) cannot
score unexpressed genes; (3) hand-tuned, not supervised; (4) fixed, non-contextual.

**DepMap/CRISPR (functional gold standard).**
*Strengths:* most direct definition of essentiality; benchmark for everything.
*Weaknesses:* requires large screens; cell lines ≠ primary tissue; cancer bias.

**Deep learning.**
*Strengths:* multi-modal features; covers unexpressed genes; cross-species; can
align to functional labels.
*Weaknesses:* poor interpretability; label-quality dependent; training bias
towards cancer cell lines.

# 5. Discussion

The field is converging on a synthesis. **Expression-robustness and functional
essentiality are complementary, not competing**: expression features are cheap,
available, and single-cell-resolved, while functional labels are expensive but
authoritative. Neither alone is sufficient. Expression-only scores mis-classify
dispensable highly-expressed genes and miss essential genes that are conditionally
expressed; function-only labels are unavailable for most tissues and species.

Deep learning offers the machinery to merge these: learn a mapping from
single-cell expression and sequence features to functional essentiality labels,
thereby transferring the authority of DepMap to any scRNA-seq dataset. The
principal risks are interpretability and label bias, which argue for
interpretable (GBM) baselines plus model-agnostic explanation, and for
explicit context modelling rather than a single global score.

# 6. Future Directions (ranked)

1. **Learn the expression score against functional labels.** Replace the
   hand-tuned `detection × expr × stability` formula with a supervised model
   (GBM → small DNN) trained on DepMap CERES labels from scRNA-derived features.
   *Highest biological-value-to-effort.*
2. **Context- / cell-type-conditional essentiality.** Output both pan-essential
   and context-specific essentials (core-vs-context framing of ADAM/CoRe).
   *Biggest novelty increment.*
3. **Sequence / DNA-LLM features.** Cover unexpressed genes and enable
   cross-species transfer (EssTFNet, ProteomeLM).
4. **Network centrality as an essentiality signal.** Replace or augment
   stability with co-expression-network degree/centrality (scRoGG, STRING).
5. **Co-essentiality modules.** Deliver functionally-enriched modules instead of
   flat lists.
6. **Built-in benchmark suite.** Reproducibly report platform-consistency,
   housekeeping/Panglao overlap, DepMap AUC, and top-k enrichment — the
   measurement stick for every other improvement.

# 7. Recommended Roadmap (v2.0)

1. Benchmark suite (enables measurement).
2. Supervised score vs. DepMap.
3. Context-conditional essentiality.
4. Network- and module-based delivery.
5. DNA-LLM / cross-species (stretch).

# 8. Verified Repository Annex (2026-08-16)

| Repository | Category |
|-----------|----------|
| huiwenzh/scEssentials, huiwenzh/scRoGG | expression-robust / co-expression |
| DepMap-Analytics/CoRe, DepMap-Analytics/ADAM2 | functional essentiality (R) |
| francescojm/ADAM, broadinstitute/depmap-portal | functional / portal |
| liulab-dfci/MAGeCK, WubingZhang/MAGeCKFlute | CRISPR upstream analysis |
| xzhang2016/DeepHE, xzhang2016/DeepSF | deep-learning essentiality |
| QIANJINYDX/EssTFNet, ucrbioinfo/DeeplyEssential | DNA-LLM / microbial DNN |
| alaa27/DeepVul, BioDataLearning/XGEP | transformer multi-task / lncRNA |
| Bitbol-Lab/ProteomeLM | protein language-model embeddings |

# References

*Placeholder — resolve before submission. Candidates:*
- Ziegenhain et al., 2017 (Mol. Cell) — single-cell platform benchmarking.
- Chen et al., 2020 — multi-platform human B-cell line data.
- DepMap / Tsherniak et al., 2017 (Nat. Cancer / Cell) — dependency mapping.
- The Tabula Sapiens Consortium, 2022 (Science); Tabula Muris Consortium, 2018.
- Tan et al. (DEGAS), Nat. Commun. 2023 — deep generative essentiality.
- Hounkpe et al., 2021 — housekeeping gene list.
- Original scEssentials manuscript (huiwenzh).
