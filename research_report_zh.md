# 必需基因（Essential Gene）鉴定 —— 数据、模型调研与方向建议

**报告日期：** 2026-08-16
**对象：** 针对 scEssentials 体系的改进（`liuwentingnz/essentialgene`）
**语言：** 中文

---

## 0. 一句话结论

"必需基因"（essential gene）有两类定义，对应两套完全不同的数据与模型栈：

1. **表达稳健型可比基因（scEssentials 的定义）** —— 只需 scRNA-seq 表达谱，
   用「检测率 × 表达量 × 稳定性」打分，本身**不需要**敲除实验数据。
2. **功能必需基因（DepMap/CRISPR 的定义）** —— 需要大规模 CRISPR-Cas9 敲除
   筛表型数据作为金标准标签，把"敲除后细胞死亡/增殖受损"的程度作为必需性。

真正能**验证并提升** scEssentials 的路线，是把第 2 类的功能标签（DepMap）
作为监督信号，来**学习**第 1 类的表达打分——这是本报告推荐的主攻方向。

---

## Part 1 —— 鉴定必需基因需要哪些数据

按用途分为 **特征数据（features）** 和 **标签数据（labels）** 两大类。

### 1.1 特征数据（模型输入，用于打分/预测）
| 数据类型 | 具体资源 | 用途 |
|----------|----------|------|
| scRNA-seq 表达谱 | Tabula Sapiens（人，多器官）、Tabula Muris Senis（鼠）、PBMC/10x 公共数据 | 计算检测率、平均表达、变异系数、细胞类型广度 —— scEssentials 的核心输入 |
| 跨平台表达谱 | Ziegenhain 2017（mESC，6 平台）、Chen 2020（人 B 细胞系，7 平台） | 检验"跨平台稳健性"（scEssentials 的核心卖点） |
| 染色质可及性 / 表观 | scATAC（Cusanovich et al.） | 关联必需基因与染色质开放程度 |
| 蛋白-蛋白互作（PPI） | STRING、BioGRID | 网络中心性可作为"重要性"特征 |
| 序列 / 基因组特征 | 序列保守性、等位基因注释、DNA 序列本身 | DeepHE/EssTFNet 用序列特征预测，可覆盖未表达基因 |
| 蛋白质语言模型嵌入 | ProteomeLM、ESM 系列 | 跨物种、跨数据集的迁移特征 |

### 1.2 标签数据（监督/验证，金标准）
| 数据类型 | 具体资源 | 用途 |
|----------|----------|------|
| **CRISPR 必需性打分** | **DepMap（Cancer Dependency Map）**：CERES 基因效应分，约 1000 个癌细胞系 × ~18k 基因 | 功能必需性的权威金标准标签 |
| 核心必需基因集合 | DepMap-Analytics/CoRe（core-fitness / common-essential 列表） | 明确的阳性集合 |
| 上下文特定必需基因 | ADAM / ADAM2 输出 | 区分"全细胞必需"与"特定细胞类型必需" |
| 传统必需基因目录 | DEG (Database of Essential Genes)、OGEE | 跨物种验证 |
| 管家基因目录 | Hounkpe / Eisenberg-Levanon 管家基因、PanglaoDB 标记基因 | 与 scEssentials 做重叠/差异分析 |

> **要点：** scEssentials 目前只用 1.1 的 表达特征（无标签）。要做"更好"的模型，
> 关键在于引入 1.2 的 **DepMap 功能标签** 来监督学习。

---

## Part 2 —— 目前最好的相关模型调研

按技术路线分类，每类给出代表工具与 GitHub 仓库（均已在本日验证）。

### 2.1 表达稳健型（与 scEssentials 同源）
| 工具 | 仓库 | 思路 |
|------|------|------|
| **scEssentials** | huiwenzh/scEssentials（本仓上游） | 检测率×表达量×稳定性，构建跨细胞类型/平台的稳健必需基因集 |
| **scRoGG**（同作者） | huiwenzh/scRoGG | 以 scEssentials 为参考集构建稳健共表达网络 |

### 2.2 功能必需性（CRISPR 金标准，最权威）
| 工具 | 仓库 | 思路 |
|------|------|------|
| **DepMap / CERES** | depmap.org（门户）、broadinstitute/depmap-portal | 大规模 CRISPR-Cas9 筛表型 → 基因效应分 |
| **CoRe** | DepMap-Analytics/CoRe | R 包：从多 CRISPR 筛中整合出 core-fitness 与 common-essential 基因 |
| **ADAM** | francescojm/ADAM | 自适应模型，区分核心必需 vs 上下文特定必需 |
| **ADAM2** | DepMap-Analytics/ADAM2 | 升级版，全基因组 CRISPR 筛必需基因鉴定 |
| **MAGeCK** | liulab-dfci/MAGeCK | 从原始 CRISPR 筛 FASTQ 识别必需基因（上游工具链） |

### 2.3 深度学习方法（用特征预测必需性）
| 工具 | 仓库 | 思路 |
|------|------|------|
| **DeepHE** | xzhang2016/DeepHE | 深度学习预测必需基因（融合序列/表达/网络特征） |
| **DeepSF** | xzhang2016/DeepSF | DeepHE 同作者 V2 |
| **EssTFNet** | QIANJINYDX/EssTFNet | 时频深度学习 + DNA 大语言模型预测人类必需基因 |
| **DeeplyEssential** | ucrbioinfo/DeeplyEssential | 微生物必需性 DNN 预测 |
| **DeepVul** | alaa27/DeepVul | 多任务 Transformer：联合必需性 + 药物响应 |
| **XGEP** | BioDataLearning/XGEP | 基于表达预测人必需基因 + lncRNA |

### 2.4 相关未能在 GitHub 直接定位的方法（需另行确认）
- **DEGAS**：深度生成模型，从单细胞 CRISPR 筛选推断基因必需性（概念真实，公开仓库归属待确认）。
- **scEVE** / **CHNet**：见于必需基因/染色质可及性文献，但 GitHub 检索未命中规范仓库。
> 提示：这些是文献中真实存在的方法，但若要用需从论文或作者主页另行获取，不能仅凭仓库名引用。

---

## Part 3 —— 对比与优缺点

### 3.1 三类路线的横向对比
| 维度 | 表达稳健型（scEssentials） | 功能必需性（DepMap/CRISPR） | 深度学习方法 |
|------|---------------------------|------------------------------|--------------|
| **数据需求** | 只要 scRNA-seq | 需要大规模敲除筛表型 | 需特征 + 标签（通常用 DepMap） |
| **可解释性** | 高（公式透明） | 高（实验验证） | 低-中（黑盒） |
| **覆盖范围** | 仅表达的基因 | 全基因组 | 取决于特征，可覆盖未表达基因 |
| **跨物种** | 需重新计算 | 难（每种都要筛） | 可迁移（序列/蛋白嵌入） |
| **能否跑在新数据上** | ✅ 是（本仓 scessentials.py 已实现） | ❌ 需已有筛表型数据 | ✅ 需训练/已有模型 |
| **生物意义** | 稳健内参/管家角色 | 直接"敲了会死" | 拟合功能标签 |

### 3.2 优缺点明细
**scEssentials（当前基线）**
- 优点：轻量、只需表达谱、可解释、能在任意新 scRNA 数据上直接跑。
- 缺点：① 非"功能必需"（表达稳健 ≠ 敲除致死）；② 无法覆盖未表达基因；
  ③ 是手调公式，未用标签学习；④ 固定打分方式，无法适配上下文/细胞类型。

**DepMap/CRISPR（功能金标准）**
- 优点：生物学上最直接的"必需"定义，是验证一切打分器的基准。
- 缺点：依赖大规模筛表型，多数研究拿不到；细胞系 ≠ 原代组织；肿瘤背景偏倚。

**深度学习方法（DeepHE/EssTFNet/DeepVul）**
- 优点：可用多维特征、可覆盖未表达基因、可跨物种迁移、可与功能标签对齐。
- 缺点：可解释性差；依赖高质量标签；数据偏差（多训练于癌细胞系）。

---

## Part 4 —— 值得做的方向（按性价比排序）

### 方向 1（首推）：向 DepMap 功能标签"学习"表达打分
- **做法**：把 scEssentials 当前的 `detection × expr_rank × stability` 换成有监督模型，
  用 DepMap CERES 分数做标签，特征取 scRNA 派生量（检测率、平均表达、变异、共表达度、
  细胞类型广度）。先用可解释的 GBM/逻辑回归，再上小 DNN。
- **理由**：直接把分数从"表达稳健"优化到"功能必需"，验证容易、价值最高。
- **参考**：DeepVul（多任务）、DeepHE（深度特征融合）。

### 方向 2：上下文 / 细胞类型条件必需性
- **做法**：按细胞类型分别打分，输出"全细胞必需"（pan-essential）+ "特定细胞必需"
  （context-specific）两级结果。
- **理由**：真实必需性高度依赖背景；这是 flat list 的天然缺陷。
- **参考**：ADAM / CoRe 的 core-vs-context 框架。

### 方向 3：加入序列 / DNA-LLM 特征，覆盖未表达基因并跨物种
- **做法**：融合 EssTFNet 式 DNA 语言模型嵌入，或进化保守性打分。
- **理由**：解决"没表达就打不了分"的盲区；跨物种迁移。
- **参考**：EssTFNet、ProteomeLM。

### 方向 4：用共表达网络中心性替代/融合"稳定性"
- **做法**：把 1/(1+CV) 换成（或叠加）共表达网络度/中心性。
- **理由**：上游调控基因比"低变异"更具生物学重要性。
- **参考**：scRoGG、STRING 拓扑特征。

### 方向 5：交付"必需基因模块"而非平铺列表
- **做法**：按跨上下文共必需性聚类成模块，配 GO/KEGG 富集。
- **理由**：更可解释、更贴近疾病机制。

### 方向 6（先修基建）：内建基准评测套件
- **做法**：在包里带 `benchmarks/`，自动在 Ziegenhain/Chen 平台集与 Tabula 图谱上
  重跑打分，汇报：平台一致性（秩 SD）、与管家/Panglao 重叠、DepMap AUC、top-k 富集。
- **理由**：没有评测尺子，所有改进都无法量化、无法发表。

---

## Part 5 —— 建议路线图（v2.0）
1. **方向 6** 基准套件（低成本，解锁一切测量）。
2. **方向 1** 学习式打分 vs DepMapAUC（生物价值/成本比最高）。
3. **方向 2** 上下文条件必需性（创新性最大的增量）。
4. **方向 4 / 5** 网络化与模块化交付（在有效分数基础上增强）。
5. **方向 3** DNA-LLM 跨物种（最雄心、最后做）。

---

## Part 6 —— 验证过的关键 GitHub 仓库清单（2026-08-16）
| 仓库 | 类别 |
|------|------|
| huiwenzh/scEssentials, huiwenzh/scRoGG | 表达稳健型 / 共表达 |
| DepMap-Analytics/CoRe, DepMap-Analytics/ADAM2 | 功能必需（R） |
| francescojm/ADAM, broadinstitute/depmap-portal | 功能必需 / 门户 |
| liulab-dfci/MAGeCK, WubingZhang/MAGeCKFlute | CRISPR 筛选上游分析 |
| xzhang2016/DeepHE, xzhang2016/DeepSF | 深度学习必需基因 |
| QIANJINYDX/EssTFNet, ucrbioinfo/DeeplyEssential | DNA-LLM / 微生物 DNN |
| alaa27/DeepVul, BioDataLearning/XGEP | Transformer 多任务 / lncRNA |
| Bitbol-Lab/ProteomeLM | 蛋白质语言模型嵌入 |
