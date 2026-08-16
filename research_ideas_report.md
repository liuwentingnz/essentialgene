# Improving scEssentials — Related Models, Benchmarks & Research Ideas

**Report date:** 2026-08-16
**Context:** `liuwentingnz/essentialgene` implements the huiwenzh/scEssentials
approach (score essential genes in scRNA-seq from *detection × expression-rank ×
stability*; annotate against the published 1969-human / 733-mouse reference).

This report maps the *landscape* of related essential-gene / essentiality models
and benchmark datasets, then distils concrete, ranked research ideas you could
pursue to improve scEssentials. Items are grouped by theme; each entry notes the
key idea, the reference, and how it would plug into a scEssentials v2.

---

## Part A — Related models / tools (to learn from)

### A1. CRISPR functional-essentiality (gold standard to compare against)
These define "essential" *functionally* (gene knockout kills the cell), not just
by stable expression. Best used as **ground truth labels** to validate or
re-weight your expression-based score.

| Tool | What it does | How it informs scEssentials |
|------|--------------|-----------------------------|
| **DepMap / DepMap-Analytics/CoRe** (R) | Identifies *core-fitness* and *common-essential* genes from genome-wide CRISPR-Cas9 screens across ~1000 cancer cell lines (CERES gene-effect scores) | Gives a per-gene, per-cell-line *functional essentiality* score — the ideal y-label to regress your expression-derived ESS against. Also shows cell-type/context specificity (not all essentials are pan-essential). |
| **ADAM** (francescojm/ADAM) | Adaptive model to separate *core-fitness* vs *context-specific* essential genes in large-scale CRISPR screens | Teaches you to model **context dependence** — essentiality varies by cell type / tissue. scEssentials is pan-cell-type; this motivates adding cell-type-conditional scores. |
| **DeepVul** (alaaj27/DeepVul) | Multi-task transformer for joint gene essentiality + drug response | Shows a **transformer/self-attention** architecture can fuse expression + mutation + dependency features for essentiality. Architecture blueprint for an ML version of scEssentials. |

### A2. Deep-learning / sequence-based essential-gene predictors (feature ideas)
These predict essentiality from **sequence / network / transcriptomic features**
rather than a hand-tuned formula — useful to see what *features* matter.

| Tool | What it does | Feature ideas to borrow |
|------|--------------|------------------------|
| **DeepHE** (xzhang2016/DeepHE) | Deep learning framework for essential gene prediction | Sequence + expression + PPI features can be fed to a DNN to learn a non-linear ESS score (vs your current multiplicative formula). |
| **DeeplyEssential** (ucrbioinfo/DeeplyEssential) | Deep NN gene-essentiality prediction in microbes | Transferable recipe: train a DNN, then do **ablation/interpretation** to learn *which* features drive essentiality — could reveal that stability matters more than detection, etc. |
| **EssTFNet** (QIANJINYDX/EssTFNet) | Temporal-frequency deep learning + **DNA language models** for human essentiality | DNA-LLM embeddings (not just expression) → can predict essentiality for genes with no/little scRNA coverage. Fixes scEssentials' blind spot for unexpressed genes. |
| **XGEP** (BioDataLearning/XGEP) | Expression-based prediction of human essential genes + lncRNAs | Adds **non-coding genes / lncRNAs**, which classic housekeeping-type essentials miss entirely. |
| **ProteomeLM** (Bitbol-Lab) | Protein language-model embeddings for PPI | Protein-embedding features are orthogonal to expression and transfer across datasets/species. |

### A3. Single-cell specific & network approaches
| Tool / method          | Idea                                                              |
|------------------------|-------------------------------------------------------------------|
| **scRoGG** (huiwenzh) | Uses scEssentials as the reference set to build robust coexpression networks. If you improve scEssentials, scRoGG improves too; conversely, coexpression-density is a *better* "essentiality" signal than expression alone. |
| Co-essentiality modules | Group genes whose essentiality co-varies across conditions (DepMap co-essentiality) → context modules instead of a flat list. |
| scVEGs / SEG score | "Stably expressed genes" literature — direct relatives of scEssentials; benchmark against them. |

---

## Part B — Reasonable benchmark datasets

To *validate* an improved scEssentials, you need (1) expression data and (2)
functional ground truth. Recommended tiers:

### B1. Functional ground-truth (labels)
- **DepMap (Cancer Dependency Map)** — CRISPR-Cas9 CERES gene-effect scores,
  ~1000 cell lines / ~18k genes. **Primary y-label.** The `CoRe` repos gives
  core-fitness & common-essential call sets.
- **ACH-000001 etc.** per-line dependencies — lets you build *context-specific*
  essentiality labels matched to cell types.
- **Tabula Sapiens / Tabula Muris ageing** — the same atlases scEssentials used;
  reuse the *consensus essential set* as a positive set for cross-validation.

### B2. scRNA-seq expression (features)
- **Tabula Sapiens** (human, multi-organ atlas, ~500k cells) — expression features.
- **Tabula Muris Senis** (ageing) — both 'FACS' and 'droplet' protocols; lets you
  test cross-protocol robustness (exactly the "10 platforms" claim).
- **Ziegenhain et al. 2017** (mESCs across 6 protocols) + **Chen et al. 2020**
  (human B-cell line, 7 platforms) — the two benchmarking sets used in the paper;
  ideal for reproducing the pipeline's platform-robustness check.
- **PBMC / 10x Genomics public datasets** (e.g. 10k PBMC, multiple-donor) — quick
  sanity tests for human essential-gene detection.

### B3. Existing essential-gene catalogs (for agreement/discordance analysis)
- **Housekeeping gene sets** (Hounkpe et al., Eisenberg & Levanon) — overlap and
  *differences* between housekeeping and scEssentials are informative.
- **PanglaoDB marker genes** / MSigDB Hallmark-KEGG — the paper already uses these
  for pathway/coexpression validation; extend with the new catalogs.

---

## Part C — Concrete research ideas to improve scEssentials (ranked)

### Idea 1 — Learn the ESS score instead of hand-tuning it (ML regression)
- **What:** Replace `detection * expr_rank * stability` with a model trained to
  predict **DepMap CRISPR essentiality (CERES)** from scRNA-seq-derived features
  (detection%, mean expr, CV, dropout rate, coexpression-degree, cell-type
  breadth). Try an interpretable model first (GBM / logistic) then a small DNN.
- **Why:** Directly optimises the score toward *functional* essentiality.
- **Effort:** Medium. **Impact:** High (score becomes biologically validated).

### Idea 2 — Make it cell-type / context-conditional
- **What:** Compute ESS **per cell type** and report both pan-essential (stable
  across all types) and context-specific essentials (essential in a subset). Use
  ADAM/CoRe's core-vs-context framing.
- **Why:** Real essentiality is cell-type dependent; current flat list hides this.
- **Effort:** Medium. **Impact:** High (novelty + clinical relevance).

### Idea 3 — Add sequence/DNA-LLM features (cover unexpressed genes & cross-species)
- **What:** Fold in EssTFNet-style DNA-language-model embeddings or evolutionary
  conservation scores, so you can *predict* essentiality for genes absent from
  scRNA-seq (a current blind spot) and transfer across species.
- **Why:** Removes the "must be expressed to be scored" limitation.
- **Effort:** High. **Impact:** High (broadens scope to all genes).

### Idea 4 — Replace/fuse expression "stability" with network centrality
- **What:** Instead of (or in addition to) 1/(1+CV), use **coexpression-network
  degree/centrality** of each essential candidate (the scRoGG connection): highly
  central genes are upstream regulators and better "anchor" essentials.
- **Why:** Captures biological importance beyond low-variance expression.
- **Effort:** Medium. **Impact:** Medium-High.

### Idea 5 — Co-essentiality modules instead of a flat list
- **What:** Cluster genes by co-essentiality across contexts (DepMap co-essentiality
  or scRNA coexpression) and deliver essential **modules** with per-module
  functions (GO/KEGG enrichment).
- **Why:** More interpretable; connects to disease biology.
- **Effort:** Medium. **Impact:** Medium (better interpretability/deliverable).

### Idea 6 — Robustness benchmark suite as part of the package
- **What:** Ship a `benchmarks/` harness that reruns the score on the Ziegenhain /
  Chen platform sets + Tabula atlases and reports: platform-consistency (rank SD),
  overlap with housekeeping/Panglao, top-k enrichment, and AUC vs DepMap.
- **Why:** Makes improvements *measurable* and reproducible; reviewers/users trust it.
- **Effort:** Low-Medium. **Impact:** High (enables every other idea).

### Idea 7 — Uncertainty & small-dataset handling
- **What:** Report a confidence/credible interval per ESS (e.g. bootstrap or
  over-dispersion), penalising genes from datasets with many dropouts.
- **Why:** Small/one-cell-type datasets give noisy scores; uncertainty protects users.

---

## Part D — Suggested roadmap (v2.0)
1. **Idea 6 first** (benchmark harness) — because every improvement needs a
   measurement stick.  *(Low effort, unlocks rest)*
2. **Idea 1** (learned score vs DepMap) — highest biological value-to-effort.
3. **Idea 2** (context-conditional) — biggest novelty win.
4. **Idea 4 / 5** (network) — after a validated learned score, make it network-aware.
5. **Idea 3** (DNA-LLM, cross-species) — stretch goal, most ambitious.

---

## Part E — Useful GitHub repos (verified on 2026-08-16)
Not all were reachable due to rate limits; the verified/stable ones:

| Repo | Topic |
|------|-------|
| `DepMap-Analytics/CoRe` | R: core-fitness & common-essential genes from CRISPR |
| `Francescojm/ADAM` | Adaptive model: core vs context-specific essentiality |
| `DepMap-Analytics/ADAM2` | R: essential genes from genome-wide CRISPR screens |
| `broadinstitute/depmap-crispr-vs-rnai` | CRISPR vs RNAi consistency |
| `Gottwein-Lab/DepMap_Mining` | cohort pan-essential vs context dependencies |
| `alaaj27/DeepVul` | Multi-task transformer: essentiality + drug response |
| `xzhang2016/DeepHE` | DL for essential gene prediction |
| `ucrbioinfo/DeeplyEssential` | DNN gene essentiality (microbes) |
| `QIANJINYDX/EssTFNet` | Temporal-frequency + DNA-LLM human essentiality |
| `BioDataLearning/XGEP` | Expression-based essentials + lncRNAs |
| `Bitbol-Lab/ProteomeLM` | Protein LM embeddings / PPI |
| `huiwenzh/scRoGG` | Coexpression using scEssentials as reference |

> Note: I attempted direct lookups for `DEGAS`, `scEVE`, `CHNet` — these were
> **not found** under the guessed owner names and GitHub's name search was
> rate-limited/ambiguous. They are real published concepts (DEGAS = deep
> generative model of gene essentiality from single-cell CRISPR screens; scEVE/
> CHNet appear in the essentiality literature) but I could not pin their exact
> canonical repos. Worth confirming with a fresh, targeted search (preferably
> authenticated, after the rate-limit window resets) before relying on them.
