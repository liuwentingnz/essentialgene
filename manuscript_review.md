---
title: "Computational Identification of Essential Genes in Single-Cell and CRISPR Data: A Review of Methods, Benchmarks, and Future Directions"
author:
  - "L. Wenting"
date: "2026-08-18"
abstract: |
  Essential genes — genes whose disruption is incompatible with cell viability — are
  central to understanding cellular fitness, cancer dependencies, and housekeeping
  functions. Recent advances in single-cell RNA sequencing (scRNA-seq) and pooled
  CRISPR-Cas9 screening have produced two largely separate literatures: (i)
  expression-robustness methods such as scEssentials that nominate "essential"
  genes from stable, broad, high expression across cell types and platforms, and
  (ii) functional-essentiality methods (DepMap, Project Score, CoRe, ADAM) that
  measure the phenotypic consequence of gene knockout. This review formalises the
  data requirements of each approach, surveys and benchmarks the current state of
  the art — including a rapidly growing family of deep-learning predictors (DeepHE,
  EPGAT, FluxGAT, LeAP, EssentialGIN, DeepHEM, DeepVul, XGEP) and scCRISPR
  generative models (scLAMBDA, PerturbNet) — and critically compares the two
  paradigms. We argue that the most promising direction is a supervised,
  context-aware synthesis that uses functional CRISPR labels (DepMap) to learn and
  reweight the expression-derived essentiality score, extended by sequence /
  DNA-language-model features for unexpressed genes and cross-species transfer, and
  delivered as co-essentiality modules rather than flat gene lists.
keywords: [essential genes, scRNA-seq, CRISPR, DepMap, deep learning, gene dependency, GNN, scCRISPR, review]
---

# 1. Introduction

Genes whose loss is lethal or strongly impairs cell fitness are termed *essential*
genes. They underpin a cell's baseline viability, serve as reference points for
normalisation in omics assays, and, in cancer, define the concept of "dependencies"
that are promising therapeutic targets. Historically, essentiality was defined by
conservation and single-gene deletion studies in microbes (DEG, OGEE, Keio). The
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
  "essentiality" or "dependency" scores (e.g. DepMap CERES, Project Score). This is
  the operational gold standard for "essential".
- **Predictive deep-learning approach (emergent).** A third, rapidly growing
  family uses machine learning — including graph neural networks (GNNs),
  transformers, and generative models — to *predict* essentiality from features
  (sequence, expression, networks, multi-omics), often trained against CRISPR
  labels.

The central open question — the subject of this review — is how the expression,
functional, and predictive paradigms relate, and how they can be reconciled.

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
| Protein–protein interaction / networks | STRING, BioGRID, metabolic networks | network centrality / topology features |
| Sequence / genomic | conservation, allele annotation, raw DNA | cover unexpressed genes; DNA-LLM |
| Protein language model | ProteomeLM, ESM | transferable, cross-species embeddings |
| Multi-omics cancer covariates | **CCLE** (expression, CNV, mutation, methylation) | context-specific dependency features |

## 2.2 Label data (ground truth)

| Label type | Representative resources | Purpose |
|------------|--------------------------|---------|
| CRISPR essentiality | **DepMap CERES** (~1100 lines × 18k genes); **Project Score** (Sanger) | functional gold-standard score |
| Core-fitness / true-essential sets | DepMap-Analytics/CoRe; **Hart/BEAN** curated "true" essential gene lists | core/common-essential positives |
| Context-specific essentiality | **ADAM / ADAM2** | core vs. context distinction |
| RNAi cross-validation | **DRIVE**, **Achilles RNAi**, **GenomeRNAi** | cross-platform agreement |
| Microbial essentiality | **Keio** (*E. coli*), DEG, OGEE | species-transfer validation |
| Essential-gene catalogs | DEG, OGEE | cross-species validation |
| Housekeeping catalogs | Hounkpe, Eisenberg–Levanon, PanglaoDB | overlap / discordance analysis |

> **Key point.** scEssentials currently uses only expression features (no labels).
> The most direct route to a *better* model is to introduce functional labels
> (DepMap / Project Score) as a supervision signal to *learn* the expression-derived
> score.

# 3. State of the Art: Method Survey

## 3.1 Expression-robustness methods

**scEssentials** (huiwenzh/scEssentials). Defines a reference set of genes stably
and highly expressed across >60 cell types and 10+ sequencing platforms. Score is
`ESS ≈ detection · expression-rank · stability`. Directly runnable on new data (the
accompanying implementation `scessentials.py`). **scRoGG** (huiwenzh/scRoGG) reuses
this set as an anchor for robust co-expression networks.

## 3.2 Functional-essentiality methods

**DepMap / CERES.** Genome-wide CRISPR-Cas9 dependency scores across ~1100 cancer
cell lines (Tsherniak 2017; Dempster/Achilles 2019). The reference label set.
**Project Score** (Sanger, ~300+ lines; Dwane 2020; Behan 2019 *Nature*) provides an
independent dependency resource for therapeutic-target prioritisation. **CoRe**
(DepMap-Analytics/CoRe) integrates multiple screens to define core-fitness and
common-essential genes. **ADAM / ADAM2** (francescojm/ADAM) classify genes as
core-fitness (essential everywhere) versus context-specific. **MAGeCK**
(liulab-dfci/MAGeCK) is the upstream tool that turns raw screen FASTQs into
essentiality calls. Population-fitness models **Chronos** (Dempster 2021) and
**ACE** (Hutton 2021) improve gene-fitness inference from CRISPR data.

## 3.3 Deep-learning predictors

| Tool | Year | Repository / arXiv | Approach |
|------|------|--------------------|----------|
| DeepHE / DeepSF | 2021 | xzhang2016/DeepHE | DNN over sequence/expression/network features |
| EPGAT | 2020 | arXiv:2007.09671 | graph attention over PPI networks |
| FluxGAT | 2024 | arXiv:2403.18666 | flux-sampling + graph attention on metabolic networks |
| DeepHEN | 2023 | arXiv:2309.10008 | deep model for essential **lncRNA** genes |
| LeAP | 2025 | arXiv:2502.15646 | layered ensemble of autoencoders + predictors |
| Topology-ML (metabolism) | 2025 | arXiv:2507.20406 | graph-topological features beat FBA |
| Multiplex-network DL | 2024 | arXiv:2403.02724 | predicting perturbation outcomes from multiplex networks |
| Explainable DL (PPI) | 2025 | arXiv:2511.12463 | centrality + node embeddings for target prioritisation |
| EssentialGIN | 2026 | arXiv:2606.07700 | graph isomorphism NN |
| DeepHEM | 2025 | Mol Ther | domain-adversarial learning, human essentials |
| DeEPsnap | 2024 | preprint | multi-omics integration |
| DeepVul | 2024 | alaa27/DeepVul | multi-task transformer: essentiality + drug response |
| XGEP | 2023 | BioDataLearning/XGEP | expression-based essentials + lncRNAs |
| DeeplyEssential | 2020 | ucrbioinfo/DeeplyEssential | microbial essentiality DNN |
| Spiking NN (CGR) | 2021 | Chaos Solitons Fractals | chaos-game representation + spiking NN |
| MetaGEM | 2026 | arXiv:2605.14812 | genome-scale metabolic networks via deep enzyme anchoring |
| dnaHNet | 2026 | arXiv:2602.10603 | hierarchical foundation model for genomic sequence |

## 3.4 Generative models & single-cell CRISPR (scCRISPR)

The current frontier is predicting *single-cell perturbation response* and
essentiality from scCRISPR screens. The following are verifiable via public APIs
(PubMed / Crossref / Europe PMC):

| Tool | Year | Venue / DOI | Approach |
|------|------|-------------|----------|
| **scLAMBDA** | 2024 | bioRxiv 10.1101/2024.12.04.626878 | deep generative framework for scCRISPR perturbation responses |
| **PerturbNet** | 2025 | Mol Syst Biol 10.1038/s44320-025-00131-3 | generative model for single-cell perturbation response |
| ACE | 2021 | Genome Biol 10.1186/s13059-021-02491-z | profile-HMM essentiality from (single-cell) CRISPR screens |
| Chronos | 2021 | Genome Biol 10.1186/s13059-021-02540-7 | cell-population fitness model |
| scCRISPR tiling | 2021 | Nat Commun 10.1038/s41467-021-24324-0 | functional inference from single-cell CRISPR tiling |
| Single-cell CRISPR essentiality assay | 2015 | Biol Proced Online (PMID:26578851) | medium-throughput single-cell CRISPR essentiality assay |

> **Verification note.** The widely-cited **DEGAS** ("A deep generative model of
> gene essentiality from single-cell CRISPR screens," attributed to Lorbeer, Tan,
> et al., Nat Commun 2023) could **not be confirmed** through PubMed, Crossref,
> Europe PMC, or arXiv title/name searches performed for this review, and no
> canonical public repository was found. It is therefore **not** included in the
> verified table above; if it is to be cited, its DOI and authorship must be
> confirmed from the original publication before submission. The same applies to
> **scEVE** and **CHNet**, which appear in the essentiality literature but could not
> be verified against public APIs or a canonical repository.

## 3.5 Methods with no canonical public repository

**DEGAS**, **scEVE**, and **CHNet** appear in the literature but could not be
verified in this survey (neither a canonical GitHub repository nor a resolvable
DOI via PubMed/Crossref/Europe PMC/arXiv). They should be sourced directly from
their papers, with DOI and authorship confirmed, before being cited.

# 4. Comparison and Strengths / Weaknesses

![Comparison of the four method paradigms across six qualitative dimensions.
(A) Radar chart. (B) Grouped dimension scores. Scores range 0–5 and reflect the
authors' qualitative assessment from the survey in Section 3.
](manuscript_figure_comparison.png){#fig:comparison width=95%}

Figure @fig:comparison summarises the qualitative trade-offs among the four
paradigms discussed below (Section 4.1) and detailed in Sections 3–5.

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

**DepMap / Project Score / CRSIPR (functional gold standard).**
*Strengths:* most direct definition of essentiality; benchmark for everything.
*Weaknesses:* requires large screens; cell lines ≠ primary tissue; cancer bias.

**Graph / topology methods (EPGAT, FluxGAT, topology-ML).**
*Strengths:* exploit network structure; good for unexpressed or under-measured
genes; FluxGAT/metabolic-topology avoid FBA bias.
*Weaknesses:* network completeness/quality dependent; require curated networks;
less mature for single-cell resolution.

**Generative / scCRISPR (scLAMBDA, PerturbNet).**
*Strengths:* capture sample-level heterogeneity and perturbation responses;
single-cell resolution; predict phenotypes for unseen genes/conditions.
*Weaknesses:* data-hungry; heavy compute; label quality from screen design;
newer and less standardised.

**Deep learning in general.**
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

The new generation of methods adds two powerful ingredients:
1. **Network/topology information** (GNNs, metabolic flux, PPI centrality) that
   captures biological importance beyond expression statistics and can score
   genes that are barely or not expressed in a given dataset.
2. **Generative single-cell modeling** (scCRISPR) that learns essentiality from
   the transcriptomic *consequences* of perturbing many genes at once, in a
   heterogeneous single-cell context — closer to the biological reality than any
   global score.

Deep learning offers the machinery to merge these: learn a mapping from
single-cell expression, sequence, and network features to functional essentiality
labels, thereby transferring the authority of DepMap/Project Score to any scRNA-seq
dataset. The principal risks are interpretability and label bias, which argue for
interpretable (GBM) baselines plus model-agnostic explanation, and for explicit
context modelling rather than a single global score.

# 6. Future Directions (ranked)

1. **Learn the expression score against functional labels.** Replace the
   hand-tuned `detection × expr × stability` formula with a supervised model
   (GBM → small DNN/GNN) trained on DepMap/Project Score labels from scRNA-derived
   + network features. *Highest biological-value-to-effort.*
2. **Context- / cell-type-conditional essentiality.** Output both pan-essential
   and context-specific essentials (core-vs-context framing of ADAM/CoRe, using
   CCLE covariates). *Biggest novelty increment.*
3. **Network centrality / topology as an essentiality signal.** Replace or augment
   stability with PPI/metabolic-network degree and centrality (EPGAT, FluxGAT,
   topology-ML, scRoGG). Enables scoring of unexpressed genes.
4. **Sequence / DNA-LLM features.** Cover unexpressed genes and enable
   cross-species transfer (EssTFNet, dnaHNet, ProteomeLM).
5. **Generative scCRISPR alignment.** Where scCRISPR data exist, align/validate
   against scLAMBDA/PerturbNet-style generative essentiality to capture perturbation
   response and heterogeneity.
6. **Co-essentiality modules.** Deliver functionally-enriched modules (with
   GO/KEGG) instead of flat lists.
7. **Built-in benchmark suite.** Reproducibly report platform-consistency,
   housekeeping/Panglao overlap, DepMap AUC, Hart/BEAN true-essential overlap, and
   top-k enrichment — the measurement stick for every other improvement.

# 7. Recommended Roadmap (v2.0)

1. Benchmark suite (enables measurement; anchor on DepMap + Hart/BEAN).
2. Supervised score vs. DepMap/Project Score.
3. Context-conditional essentiality (with CCLE covariates).
4. Network- and module-based delivery (GNN / topology).
5. DNA-LLM / cross-species (stretch).

# 8. Verified Repository & Dataset Annex (2026-08-18)

**Repositories**
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

**Datasets**
| Resource | Contents | Access |
|----------|----------|--------|
| DepMap / Project Achilles | ~1100 lines × 18k genes, CRISPR + RNAi | depmap.org/portal |
| Project Score (Sanger) | ~300+ lines, dependency scores | score.depmap.sanger.ac.uk |
| CCLE | ~1000 lines genomic covariates | sites.broadinstitute.org/ccle |
| DRIVE + Achilles RNAi | RNAi dependency cross-validation | DepMap portal |
| Keio (*E. coli*) | microbial essentiality gold standard | KHRI/NBRP |
| Hart / BEAN | curated "true" essential gene lists | DepMap downloads |
| Perturb-seq / scCRISPR | single-cell perturbation essentiality | GEO / Broad portal |
| GenomeRNAi | RNAi-based essentiality reference | genomernai.org |
| L1000 / CMAP | perturbational signatures | clue.io |

# References

## Benchmarks & functional essentiality
- Tsherniak A, et al. Defining a Cancer Dependency Map. *Cell* 2017.
- Dempster JM, et al. Extracting biological insights from Project Achilles genome-scale
  CRISPR screens in cancer cell lines. *Cancer Discov* 2019. DOI:10.1158/2159-8290.CD-18-0834.
- Behan FM, et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9
  screens. *Nature* 2019. DOI:10.1038/s41586-019-1103-9.
- Dwane L, et al. Project Score database. *Nucleic Acids Res* 2020. DOI:10.1093/nar/gkaa882.
- Barretina J, et al. The Cancer Cell Line Encyclopedia (CCLE). *Nature* 2012. DOI:10.1038/nature11003.
- McDonald ER, et al. Project DRIVE. *Cell* 2017. DOI:10.1016/j.cell.2017.05.006.
- Baba T, et al. Keio collection. *Mol Syst Biol* 2006. DOI:10.1038/msb4100050.
- Schmidt EE, et al. GenomeRNAi. *Nucleic Acids Res* 2013. DOI:10.1093/nar/gks1171.
- Dixit A, et al. Perturb-Seq. *Cell* 2016. DOI:10.1016/j.cell.2016.11.038.
- Subramanian A, et al. L1000/CMAP. *Cell* 2017. DOI:10.1016/j.cell.2017.10.049.

## Inference models
- Dempster JM, et al. Chronos. *Genome Biol* 2021. DOI:10.1186/s13059-021-02540-7.
- Hutton ER, et al. ACE. *Genome Biol* 2021. DOI:10.1186/s13059-021-02491-z.
- Yang et al. High-resolution scCRISPR tiling. *Nat Commun* 2021. DOI:10.1038/s41467-021-24324-0.
- Zhao Y, Zhang M, Yang D. Bioinformatics approaches to CRISPR screen data. *Quant Biol*
  2022. DOI:10.15302/j-qb-022-0299.

## Generative / scCRISPR
- Wang et al. scLAMBDA. *bioRxiv* 2024. DOI:10.1101/2024.12.04.626878.
- Yu et al. PerturbNet. *Mol Syst Biol* 2025. DOI:10.1038/s44320-025-00131-3.

## Deep-learning predictors (arXiv-verified)
- Schapke, Tavares, Recamonde-Mendoza. EPGAT. *arXiv:2007.09671* (2020).
- Sharma, Marucci, Abdallah. FluxGAT. *arXiv:2403.18666* (2024).
- Zhang, Cheng. DeepHEN. *arXiv:2309.10008* (2023).
- Bodinier et al. LeAP. *arXiv:2502.15646* (2025).
- Boone. Topology-based metabolic essentiality ML. *arXiv:2507.20406* (2025).
- Zhan, Zhang, Wang. Multiplex-network perturbation DL. *arXiv:2403.02724* (2024).
- Alkhadrawi et al. Explainable DL for target prioritisation. *arXiv:2511.12463* (2025).
- Mansouri-Rad, Narimani, Razzaghi. EssentialGIN. *arXiv:2606.07700* (2026).
- Zhan. Genome-scale deep model on multiplex networks. *arXiv:2403.02724* (2024).
- Xiao, Zheng, Li. MetaGEM. *arXiv:2605.14812* (2026).
- Shah, Li. dnaHNet. *arXiv:2602.10603* (2026).

## Others (from direct search)
- Zhang, Cheng. DeepHEM. *Mol Ther* 2025. DOI:10.1016/j.ymthe.2025.11.018.
- Gene essentiality via chaos-game representation + spiking NN. *Chaos Solitons Fractals*
  2021. DOI:10.1016/j.chaos.2021.110649.
- Deep gene essentiality based on population/functions. *bioRxiv* 2021. DOI:10.1101/2021.12.21.473690.
- DeEPsnap (multi-omics essentiality). *bioRxiv* 2024. DOI:10.1101/2024.06.20.599958.
- Deep learning-based representation methods for survival & essentiality. *Sci Rep*
  2024. DOI:10.1038/s41598-024-67023-8.
