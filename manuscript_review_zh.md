---
title: "单细胞与 CRISPR 数据中必需基因的计算鉴定：方法、基准与未来方向综述"
author:
  - "L. Wenting"
date: "2026-08-18"
abstract: |
  必需基因——其功能丧失会使细胞无法存活的基因——是理解细胞适应性、癌症依赖
  （dependency）以及管家功能的核心。单细胞 RNA 测序（scRNA-seq）与混合 CRISPR-Cas9
  筛选技术的进步，催生了两条在研究上相对独立的文献脉络：（i）基于表达稳健性的方法，
  如 scEssentials，通过在细胞类型与测序平台间稳定、广泛、高表达来提名"必需"基因；
  （ii）基于功能必需性的方法（DepMap、Project Score、CoRe、ADAM），直接度量基因敲
  除后的表型后果。本综述形式化了两类方法的数据需求，调研并对比了当前最先进的方法
  ——包括快速增长的深度学习预测器家族（DeepHE、EPGAT、FluxGAT、LeAP、EssentialGIN、
  DeepHEM、DeepVul、XGEP）以及 scCRISPR 生成模型（scLAMBDA、PerturbNet）——并对两类
  范式进行了批判性比较。我们认为，最有前景的方向是一种有监督、上下文感知的综合方法：
  利用功能 CRISPR 标签（DepMap）来学习并重新加权基于表达的必需性分数，辅以序列 /
  DNA 大语言模型特征覆盖未表达基因并实现跨物种迁移，最终以共必需性模块（co-essentiality
  modules）而非扁平基因列表的形式交付。
keywords: [必需基因, 单细胞 RNA-seq, CRISPR, DepMap, 深度学习, 基因依赖, GNN, scCRISPR, 综述]
---

# 1. 引言

功能丧失会致死或显著损害细胞适应性的基因被称为*必需基因*（essential gene）。它们是
细胞基线存活能力的基础，在组学分析中常被用作标准化参照，在肿瘤中则定义了对治疗
有前景的"依赖"（dependency）概念。历史上，必需性主要是通过保守性以及微生物中的
单基因敲除研究（DEG、OGEE、Keio）来定义的。单细胞转录组学与混合 CRISPR 筛选技术
的发展，彻底改变了哺乳动物及复杂系统中必需性的度量与预测方式。

相关研究大致沿着三条既独立又有交叉的脉络发展：

- **基于表达稳健性的方法。** 诸如 **scEssentials** 这类方法，识别那些在多种细胞类型
  与测序平台间*稳定、广泛且高表达*的基因。其隐含假设是：一个基因在所有地方都以高
  水平、低变异表达，那么它很可能对一般细胞功能是必需的（一种针对单细胞数据改良的
  "管家"定义）。
- **基于功能必需性的方法。** 混合 **CRISPR-Cas9** 筛选直接度量基因敲除对细胞适应性的
  影响，得到每个基因的"必需性"或"依赖"分数（如 DepMap 的 CERES、Project Score）。
  这是"必需"的操作性金标准。
- **预测性深度学习方法（新兴）。** 第三个、且正在快速发展的家族使用机器学习——包括
  图神经网络（GNN）、Transformer 与生成模型——基于特征（序列、表达、网络、多组学）
  来*预测*必要性，通常以 CRISPR 标签作为训练目标。

本综述围绕的核心开放问题是：表达、功能与预测三类范式如何相互关联、如何被统一。

本综述有三点贡献：

1. **数据需求**——每种方法需要什么输入（第 2 节）；
2. **最先进方法对比**——对当前基于表达、基于功能与深度学习方法进行调研与横向对比，
   并给出优缺点（第 3–5 节）；
3. **未来方向**——下一代必需基因工具按优先级排序的路线图（第 6 节）。

# 2. 鉴定必需基因的数据需求

我们将输入分为*特征*（用于打分或预测）与*标签*（用于监督或验证的金标准）。

## 2.1 特征数据（模型输入）

| 模态 | 代表性资源 | 用途 |
|------|-----------|------|
| scRNA-seq 表达谱 | Tabula Sapiens、Tabula Muris Senis、10x PBMC | 检测率、平均表达、变异系数、细胞类型广度 |
| 跨平台表达谱 | Ziegenhain 2017（6 平台）、Chen 2020（7 平台） | 跨平台稳健性 |
| 染色质可及性 | scATAC（Cusanovich 等） | 与开放染色质的关联 |
| 蛋白-蛋白互作 / 网络 | STRING、BioGRID、代谢网络 | 网络中心性 / 拓扑特征 |
| 序列 / 基因组特征 | 保守性、等位基因注释、原始 DNA | 覆盖未表达基因；DNA-LLM |
| 蛋白语言模型 | ProteomeLM、ESM | 可迁移、跨物种嵌入 |
| 癌症多组学协变量 | **CCLE**（表达、拷贝数、突变、甲基化） | 上下文特定依赖特征 |

## 2.2 标签数据（金标准）

| 标签类型 | 代表性资源 | 用途 |
|----------|-----------|------|
| CRISPR 必需性 | **DepMap CERES**（约 1100 系 × 1.8 万基因）；**Project Score**（Sanger） | 功能金标准分数 |
| 核心必需基因集 | DepMap-Analytics/CoRe；**Hart/BEAN** 人工整理的"真"必需基因表 | 核心 / 通用必需阳性集 |
| 上下文特定必需性 | **ADAM / ADAM2** | 核心 vs. 上下文区分 |
| RNAi 交叉验证 | **DRIVE**、**Achilles RNAi**、**GenomeRNAi** | 跨平台一致性 |
| 微生物必需性 | **Keio**（大肠杆菌）、DEG、OGEE | 跨物种验证 |
| 必需基因目录 | DEG、OGEE | 跨物种验证 |
| 管家基因目录 | Hounkpe、Eisenberg–Levanon、PanglaoDB | 重叠 / 差异分析 |

> **要点。** scEssentials 目前仅使用表达特征（无标签）。要做出*更好的*模型，最直接的
> 途径是引入功能标签（DepMap / Project Score）作为监督信号，来*学习*基于表达的分数。

# 3. 最先进方法调研

## 3.1 基于表达稳健性的方法

**scEssentials**（huiwenzh/scEssentials）。定义了一组在超过 60 种细胞类型与 10+ 种测序
平台上稳定、高表达的参照基因集。分数为 `ESS ≈ 检测率 · 表达秩 · 稳定性`。可直接用于
新数据（配套实现为 `scessentials.py`）。**scRoGG**（huiwenzh/scRoGG）以这一基因为锚点
构建稳健的共表达网络。

## 3.2 基于功能必需性的方法

**DepMap / CERES。** 覆盖约 1100 个癌细胞系的基因组规模 CRISPR-Cas9 依赖分数
（Tsherniak 2017；Dempster/Achilles 2019），是参照标签集。**Project Score**（Sanger，
约 300+ 系；Dwane 2020；Behan 2019 *Nature*）提供了独立的依赖资源，用于治疗靶点优先排序。
**CoRe**（DepMap-Analytics/CoRe）整合多个筛表以定义核心适应性与通用必需基因。
**ADAM / ADAM2**（francescojm/ADAM）将基因分为核心适应性（处处必需）与上下文特定两类。
**MAGeCK**（liulab-dfci/MAGeCK）是将原始筛表 FASTQ 转化为必需性判定的上游工具。
群体适应性模型 **Chronos**（Dempster 2021）与 **ACE**（Hutton 2021）改进了 CRISPR
数据中基因适应性效应的推断。

## 3.3 深度学习方法

| 工具 | 年份 | 仓库 / arXiv | 方法 |
|------|------|--------------|------|
| DeepHE / DeepSF | 2021 | xzhang2016/DeepHE | 序列/表达/网络特征的 DNN |
| EPGAT | 2020 | arXiv:2007.09671 | 在 PPI 网络上的图注意力 |
| FluxGAT | 2024 | arXiv:2403.18666 | 流量采样 + 代谢网络图注意力 |
| DeepHEN | 2023 | arXiv:2309.10008 | 必需 **lncRNA** 基因的深度模型 |
| LeAP | 2025 | arXiv:2502.15646 | 自编码器 + 预测器分层集成 |
| 拓扑 ML（代谢） | 2025 | arXiv:2507.20406 | 图拓扑特征优于 FBA |
| 多重网络 DL | 2024 | arXiv:2403.02724 | 从多重网络预测扰动结果 |
| 可解释 DL（PPI） | 2025 | arXiv:2511.12463 | 中心性 + 节点嵌入用于靶点优先化 |
| EssentialGIN | 2026 | arXiv:2606.07700 | 图同构神经网络 |
| DeepHEM | 2025 | Mol Ther | 对抗域自适应，人类必需基因 |
| DeEPsnap | 2024 | 预印本 | 多组学整合 |
| DeepVul | 2024 | alaa27/DeepVul | 多任务 Transformer：必需性 + 药物反应 |
| XGEP | 2023 | BioDataLearning/XGEP | 基于表达的必需基因 + lncRNA |
| DeeplyEssential | 2020 | ucrbioinfo/DeeplyEssential | 微生物必需性 DNN |
| 脉冲神经网络（CGR） | 2021 | Chaos Solitons Fractals | 混沌博弈表示 + 脉冲 NN |
| MetaGEM | 2026 | arXiv:2605.14812 | 深度酶-代谢物锚定的基因组尺度代谢网络 |
| dnaHNet | 2026 | arXiv:2602.10603 | 基因组序列的分层基础模型 |

## 3.4 生成模型与单细胞 CRISPR（scCRISPR）

当前的前沿是用 scCRISPR 筛选预测*单细胞扰动反应*与必需性。以下均已通过公开 API
（PubMed / Crossref / Europe PMC）验证：

| 工具 | 年份 | 期刊 / DOI | 方法 |
|------|------|-----------|------|
| **scLAMBDA** | 2024 | bioRxiv 10.1101/2024.12.04.626878 | scCRISPR 扰动反应的深度生成框架 |
| **PerturbNet** | 2025 | Mol Syst Biol 10.1038/s44320-025-00131-3 | 单细胞扰动反应的生成模型 |
| ACE | 2021 | Genome Biol 10.1186/s13059-021-02491-z | 基于（单细胞）CRISPR 筛的 profile-HMM 必需性 |
| Chronos | 2021 | Genome Biol 10.1186/s13059-021-02540-7 | 细胞群体动力学适应性模型 |
| scCRISPR tiling | 2021 | Nat Commun 10.1038/s41467-021-24324-0 | 单细胞 CRISPR tiling 的功能推断 |
| 单细胞 CRISPR 必需性实验 | 2015 | Biol Proced Online（PMID:26578851） | 中通量单细胞 CRISPR 必需性测定 |

> **验证说明。** 广被引用的 **DEGAS**（"A deep generative model of gene essentiality from
> single-cell CRISPR screens"，归属于 Lorbeer、Tan 等，Nat Commun 2023）在本综述对
> PubMed、Crossref、Europe PMC、arXiv 的标题/名称检索中**无法得到确认**，也未找到其
> 规范的公开仓库。因此**未**列在上表的已验证条目中；若需引用，必须在投稿前从原始
> 出版物核实其 DOI 与作者。**scEVE** 与 **CHNet** 同理——它们出现在必需性相关文献中，
> 但无法通过公开 API 或规范仓库得到验证。

## 3.5 无规范公开仓库的方法

**DEGAS**、**scEVE** 与 **CHNet** 出现在文献中，但本综述无法验证（既无规范 GitHub
仓库，也无法通过 PubMed/Crossref/Europe PMC/arXiv 解析出 DOI）。它们应从各自论文直接
获取，并在引用前核实 DOI 与作者信息。

# 4. 对比与优缺点

![四类方法范式在六个定性维度上的对比。（A）雷达图。（B）分组维度得分。分数范围为
0–5，反映作者基于第 3 节调研的定性评估。
](manuscript_figure_comparison.png){#fig:comparison width=95%}

图 @fig:comparison 总结了下面 4.1 节讨论、并由第 3–5 节详述的四类范式之间的定性权衡。

## 4.1 横向对比

| 维度 | 基于表达稳健（scEssentials） | 基于功能（DepMap/CRISPR） | 深度学习 |
|------|------------------------------|---------------------------|----------|
| 数据需求 | 仅 scRNA-seq | 混合敲除筛 | 特征 + 标签 |
| 可解释性 | 高（公式透明） | 高（实验验证） | 低–中 |
| 基因组覆盖 | 仅表达基因 | 全基因组 | 取决于特征 |
| 跨物种 | 需重算 | 难（需逐物种筛） | 可迁移 |
| 可跑新数据 | 是 | 否（需筛表数据） | 需已训练模型 |
| 生物意义 | 管家 / 参照 | 直接敲除致死 | 拟合功能标签 |

## 4.2 优缺点明细

**scEssentials（基线）。**
*优点：*轻量、仅需表达谱、可解释、可在任意新 scRNA-seq 数据上直接运行。
*缺点：*①"表达稳健"≠"功能必需"；②无法给未表达基因打分；③手调公式、无监督；
④固定、非上下文相关。

**DepMap / Project Score / CRISPR（功能金标准）。**
*优点：*对必需性最直接的定义，是一切方法验证的基准。
*缺点：*需要大规模筛表；细胞系 ≠ 原代组织；存在肿瘤偏倚。

**图 / 拓扑方法（EPGAT、FluxGAT、拓扑 ML）。**
*优点：*利用网络结构；对未表达或低表达基因表现好；FluxGAT/代谢拓扑避免 FBA 偏差。
*缺点：*依赖网络完整性/质量；需要精选网络；对单细胞分辨率尚不成熟。

**生成式 / scCRISPR（scLAMBDA、PerturbNet）。**
*优点：*捕捉样本层面的异质性与扰动反应；单细胞分辨率；可预测未见基因/条件的表型。
*缺点：*数据饥渴；计算量大；标签质量依赖筛表设计；较新、标准化不足。

**深度学习总体。**
*优点：*多模态特征；覆盖未表达基因；可跨物种；可与功能标签对齐。
*缺点：*可解释性差；依赖标签质量；存在向癌细胞系偏移的训练偏差。

# 5. 讨论

该领域正趋于综合。**基于表达的稳健性与基于功能的必需性是互补而非竞争**：表达特征
便宜、易得、且具有单细胞分辨率；功能标签昂贵但权威。两者单独都不够。仅基于表达的分
数会误分类"高表达但可丢弃"的基因，也会漏掉"条件性表达"的必需基因；仅基于功能的标签
对大多数组织与物种不可用。

新一代方法增加了两个有力要素：
1. **网络 / 拓扑信息**（GNN、代谢通量、PPI 中心性）能捕捉超越表达统计量的生物重要
   性，并能对给定数据集中几乎不表达或未表达的基因打分。
2. **生成式单细胞建模**（scCRISPR）从同时扰动许多基因的转录组*后果*中学习必需性，
   处于异质的单细胞背景下——这比任何全局分数都更接近生物现实。

深度学习提供了融合这些要素的机制：学习从单细胞表达、序列与网络特征到功能必需性标签
的映射，从而将 DepMap/Project Score 的权威性迁移到任意 scRNA-seq 数据集。主要风险是
可解释性与标签偏倚，这要求在可解释（GBM）基线上增加模型无关的解释方法，并针对上下文
显式建模，而非使用单一全局分数。

# 6. 未来方向（按优先级排序）

1. **向功能标签学习表达分数。** 用有监督模型（GBM → 小 DNN/GNN，以 DepMap/Project
   Score 为标签，特征取 scRNA 派生 + 网络特征）替换手调的 `检测率 × 表达 × 稳定性`
   公式。*生物价值 / 成本比最高。*
2. **上下文 / 细胞类型条件必需性。** 同时输出全细胞必需与上下文特定必需（借鉴
   ADAM/CoRe 的核心 vs. 上下文框架，并利用 CCLE 协变量）。*创新增量最大。*
3. **网络中心性 / 拓扑作为必需性信号。** 用 PPI/代谢网络的度与中心性替换或增强稳定
   性（EPGAT、FluxGAT、拓扑 ML、scRoGG）。可对未表达基因打分。
4. **序列 / DNA-LLM 特征。** 覆盖未表达基因并实现跨物种迁移（EssTFNet、dnaHNet、
   ProteomeLM）。
5. **生成式 scCRISPR 对齐。** 在有 scCRISPR 数据的场景中，对照 scLAMBDA/PerturbNet
   式生成必需性进行对齐 / 验证，以捕捉扰动反应与异质性。
6. **共必需性模块。** 交付带 GO/KEGG 富集的功能模块，而非扁平列表。
7. **内置基准评测套件。** 可复现地报告平台一致性、与管家基因/Panglao 重叠度、DepMap
   AUC、Hart/BEAN 真必需基因重叠度与 top-k 富集——这是所有其他改进的度量标尺。

# 7. 建议路线图（v2.0）

1. 基准评测套件（先建立度量标尺；锚定 DepMap + Hart/BEAN）。
2. 有监督分数 vs. DepMap/Project Score。
3. 上下文条件必需性（利用 CCLE 协变量）。
4. 网络与模块化交付（GNN / 拓扑）。
5. DNA-LLM / 跨物种（最有雄心的远期目标）。

# 8. 已验证仓库与数据集附录（2026-08-18）

**仓库**
| 仓库 | 类别 |
|------|------|
| huiwenzh/scEssentials、huiwenzh/scRoGG | 表达稳健 / 共表达 |
| DepMap-Analytics/CoRe、DepMap-Analytics/ADAM2 | 功能必需性（R） |
| francescojm/ADAM、broadinstitute/depmap-portal | 功能 / 门户 |
| liulab-dfci/MAGeCK、WubingZhang/MAGeCKFlute | CRISPR 上游分析 |
| xzhang2016/DeepHE、xzhang2016/DeepSF | 深度学习方法 |
| QIANJINYDX/EssTFNet、ucrbioinfo/DeeplyEssential | DNA-LLM / 微生物 DNN |
| alaa27/DeepVul、BioDataLearning/XGEP | Transformer 多任务 / lncRNA |
| Bitbol-Lab/ProteomeLM | 蛋白语言模型嵌入 |

**数据集**
| 资源 | 内容 | 获取 |
|------|------|------|
| DepMap / Project Achilles | 约 1100 系 × 1.8 万基因，CRISPR + RNAi | depmap.org/portal |
| Project Score（Sanger） | 约 300+ 系，依赖分数 | score.depmap.sanger.ac.uk |
| CCLE | 约 1000 系基因组协变量 | sites.broadinstitute.org/ccle |
| DRIVE + Achilles RNAi | RNAi 依赖交叉验证 | DepMap 门户 |
| Keio（大肠杆菌） | 微生物必需性金标准 | KHRI/NBRP |
| Hart / BEAN | 人工整理的"真"必需基因表 | DepMap 下载 |
| Perturb-seq / scCRISPR | 单细胞扰动必需性 | GEO / Broad 门户 |
| GenomeRNAi | RNAi 必需性参照 | genomernai.org |
| L1000 / CMAP | 扰动表达谱 | clue.io |

# 参考文献

## 基准与功能必需性
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

## 推断模型
- Dempster JM, et al. Chronos. *Genome Biol* 2021. DOI:10.1186/s13059-021-02540-7.
- Hutton ER, et al. ACE. *Genome Biol* 2021. DOI:10.1186/s13059-021-02491-z.
- Yang et al. High-resolution scCRISPR tiling. *Nat Commun* 2021. DOI:10.1038/s41467-021-24324-0.
- Zhao Y, Zhang M, Yang D. Bioinformatics approaches to CRISPR screen data. *Quant Biol*
  2022. DOI:10.15302/j-qb-022-0299.

## 生成式 / scCRISPR
- Wang et al. scLAMBDA. *bioRxiv* 2024. DOI:10.1101/2024.12.04.626878.
- Yu et al. PerturbNet. *Mol Syst Biol* 2025. DOI:10.1038/s44320-025-00131-3.

## 深度学习预测器（已由 arXiv 验证）
- Schapke, Tavares, Recamonde-Mendoza. EPGAT. *arXiv:2007.09671* (2020).
- Sharma, Marucci, Abdallah. FluxGAT. *arXiv:2403.18666* (2024).
- Zhang, Cheng. DeepHEN. *arXiv:2309.10008* (2023).
- Bodinier et al. LeAP. *arXiv:2502.15646* (2025).
- Boone. Topology-based metabolic essentiality ML. *arXiv:2507.20406* (2025).
- Zhan, Zhang, Wang. Multiplex-network perturbation DL. *arXiv:2403.02724* (2024).
- Alkhadrawi et al. Explainable DL for target prioritisation. *arXiv:2511.12463* (2025).
- Mansouri-Rad, Narimani, Razzaghi. EssentialGIN. *arXiv:2606.07700* (2026).
- Xiao, Zheng, Li. MetaGEM. *arXiv:2605.14812* (2026).
- Shah, Li. dnaHNet. *arXiv:2602.10603* (2026).

## 其他（直接检索获得）
- Zhang, Cheng. DeepHEM. *Mol Ther* 2025. DOI:10.1016/j.ymthe.2025.11.018.
- Gene essentiality via chaos-game representation + spiking NN. *Chaos Solitons Fractals*
  2021. DOI:10.1016/j.chaos.2021.110649.
- Deep gene essentiality based on population/functions. *bioRxiv* 2021. DOI:10.1101/2021.12.21.473690.
- DeEPsnap (multi-omics essentiality). *bioRxiv* 2024. DOI:10.1101/2024.06.20.599958.
- Deep learning-based representation methods for survival & essentiality. *Sci Rep*
  2024. DOI:10.1038/s41598-024-67023-8.
