# 经验网络约束的合成依赖优先排序

[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

一个系统性的生物信息学框架，通过整合统计关联与多层经验生物网络证据，从大规模突变和 CRISPR 依赖数据中优先排序候选**合成致死/合成依赖**基因对。

---

## 目录

1. [生物学背景](#生物学背景)
2. [方法概览](#方法概览)
3. [项目结构](#项目结构)
4. [环境要求](#环境要求)
5. [快速开始](#快速开始)
6. [输入数据格式](#输入数据格式)
7. [管线工作流](#管线工作流)
8. [命令行用法详解](#命令行用法详解)
9. [输出结果说明](#输出结果说明)
10. [配置文件参考](#配置文件参考)
11. [评分体系详解](#评分体系详解)
12. [证据分类](#证据分类)
13. [可视化](#可视化)
14. [模拟数据](#模拟数据)
15. [零模型与置换检验](#零模型与置换检验)
16. [生物学解读模板](#生物学解读模板)
17. [实验验证策略](#实验验证策略)
18. [常见问题排查](#常见问题排查)
19. [已知局限](#已知局限)
20. [创新点定位](#创新点定位)
21. [参考文献](#参考文献)

---

## 生物学背景

肿瘤基因组携带大量突变，但其中只有一部分能产生**条件性脆弱性**——即肿瘤细胞对特定基因的依赖增强。鉴定这些突变条件性依赖是精准肿瘤学的核心挑战。

本框架解决的核心生物学问题是：

> **如果基因 A 的突变与基因 B 的依赖性增强相关，那么经验生物网络能否帮助判断这种关联在机制上是否合理、在治疗上是否可行？**

纯统计关联可能产生误导，因为存在以下混杂因素：

- 肿瘤谱系效应和批次效应
- 肿瘤突变负荷 (TMB) 和基因组不稳定性
- 微卫星不稳定性 (MSI)
- 拷贝数变异 (CNV) 负荷
- 靶基因表达差异
- 低突变频率（样本量不足）
- Hub 基因假象
- 乘客突变效应

本框架通过在统计关联之上叠加**经验生物网络证据**，产生可解释的、有机制基础的候选排序结果。

---

## 方法概览

框架将原始突变和依赖数据转化为排序的、生物学上可解释的合成依赖候选对，通过九个顺序模块实现：

```
突变矩阵  +  CRISPR 依赖矩阵  +  样本元数据
         │                    │
         └────────────────────┼────────────────────────┘
                              │
                   回归关联分析（校正谱系/TMB/MSI/CNV）
                              │
                   候选基因对 (FDR < 0.05)
                              │
         ┌────────────────────┼────────────────────────────┐
         │                    │                            │
   STRING 网络             Reactome 通路              OmniPath 信号传导
   (邻近性)                (一致性)                   (方向性)
         │                    │                            │
         ├────────────────────┼────────────────────────────┤
         │                    │                            │
   CORUM 复合物          SynLethDB 先验              药物靶点证据
   (共复合体)            (已知 SL 关系)              (可药性)
         │                    │                            │
         └────────────────────┼────────────────────────────┘
                              │
                   加权多证据评分整合
                              │
                   排序的合成依赖候选对
                              │
         ┌────────────────────┼────────────────────────────┐
         │                    │                            │
   机制卡片生成          发表级图表                HTML/MD 报告
```

### 证据层

| 证据层 | 数据来源 | 测量维度 | 分值范围 |
|--------|----------|----------|----------|
| **统计** | 关联分析 | 突变→依赖效应量 × 显著性 | 0–1（缩放） |
| **网络邻近** | STRING v12 | 最短路径距离、置信度、Hub 惩罚 | 0–1（复合） |
| **通路一致性** | Reactome | 共享生物学通路 | 0 / 0.5 / 1 |
| **方向性信号** | OmniPath / SIGNOR | 有向调控路径 A→B 或 B→A | 0 / 0.5 / 0.7 / 1 |
| **蛋白复合物** | CORUM | 共享蛋白复合物成员 | 0 / 1 |
| **SL 先验** | SynLethDB | 已知或邻居支持的 SL 证据 | 0 / 0.3 / 0.5 / 1 |
| **共依赖** | DepMap | 依赖谱的 Pearson 相关性 | 0–1 (abs(r)) |
| **可药性** | DGIdb / Open Targets | 靶基因的治疗可行性 | 0 / 0.5 / 0.7 / 0.8 / 1 |

### 最终评分公式（平衡方案）

```
最终_SL_评分 = 0.35 × 统计分
            + 0.25 × 网络分
            + 0.15 × 通路分
            + 0.10 × 信号分
            + 0.05 × 复合物分
            + 0.05 × SL先验分
            + 0.05 × 可药性分
```

提供三种评分方案：**balanced**（平衡）、**mechanism**（机制导向，侧重网络/通路）、**therapeutic**（治疗导向，侧重可药性）。

### 回归模型

对每个突变基因 A × 依赖基因 B 对：

```
依赖_B ~ 突变_A + 表达_B + 谱系 + TMB + MSI状态 + CNV负荷
```

校正 CRISPR 依赖数据中已知的主要混杂因素。

---

## 项目结构

```
network_synthetic_dependency_project/
│
├── config/
│   └── config.yaml                    # 唯一配置来源：所有可调参数
│
├── R/
│   ├── 00_utils.R                     # 共享工具函数（缩放、安全路径、日志）
│   ├── 01_prepare_data.R              # 数据加载、验证、基因符号标准化
│   ├── 02_association_analysis.R      # LM 回归、FDR 校正、候选过滤
│   ├── 03_build_networks.R            # STRING/信号网络构建、查找索引
│   ├── 04_compute_network_features.R  # 最短路径、路径置信度、Hub 惩罚
│   ├── 05_compute_empirical_evidence.R # 通路、信号、复合物、SL先验、药物评分
│   ├── 06_score_candidates.R          # 加权整合、证据分类 (I–IV)
│   ├── 07_null_model_permutation.R    # 度匹配经验 p 值
│   ├── 08_visualization.R             # 点图、热图、网络路径图
│   ├── 09_generate_report.R           # 机制卡片、HTML/Markdown 报告
│   ├── simulate_data.R                # 真实合成测试数据生成器
│   └── run_pipeline.R                 # CLI 编排器，支持步骤控制和断点续跑
│
├── data/
│   ├── raw/                           # 输入数据（用户提供或 simulate_data.R 生成）
│   └── processed/                     # 管线输出（自动生成）
│
├── results/
│   ├── tables/                        # 最终排序表格
│   ├── figures/                       # 发表级 PDF/PNG 图表
│   └── reports/                       # 自动生成的报告
│
├── logs/
├── README.md                          # 英文 README
├── README_CN.md                       # 中文 README（本文件）
└── project.md                         # 完整规范文档
```

---

## 环境要求

### 系统要求

- **R** >= 4.0
- **操作系统**：Linux、macOS 或 Windows
- **内存**：全基因组分析建议 >= 8 GB
- **磁盘**：约 500 MB（取决于输入规模）

### R 包安装

```r
# 必须安装
install.packages(c("dplyr", "ggplot2", "tidyr", "purrr", "readr", "igraph", "yaml", "optparse"))

# 推荐安装
install.packages(c("pheatmap"))   # 增强热图（带行注释）
install.packages(c("furrr"))       # 并行处理（加速关联分析步骤）
```

各包用途：
- `igraph` — 网络分析（图构建、最短路径、中心性）
- `yaml` — 配置文件解析
- `optparse` — 命令行参数解析
- `dplyr`, `tidyr`, `purrr`, `readr` — 数据操作
- `ggplot2` — 数据可视化

---

## 快速开始

```bash
# 进入项目目录
cd network_synthetic_dependency_project

# 第一步：生成模拟测试数据
Rscript R/run_pipeline.R --simulate --stop-at 1

# 第二步：运行核心管线
Rscript R/run_pipeline.R --verbose

# 一键运行：生成数据 + 全管线
Rscript R/run_pipeline.R --simulate --verbose

# 查看 Top 10 候选
head data/processed/final_ranked_candidates.csv

# 打开自动生成的报告
open results/reports/synthetic_dependency_report.html

# 查看所有生成的图表
ls results/figures/
```

---

## 输入数据格式

### 必选文件

#### 1. 突变矩阵 (`mutation_matrix.csv`)

二值矩阵（样本 × 基因），1 = 存在非同义突变，0 = 无。

| sample_id | TP53 | KRAS | ATM | ARFGEF3 | ... |
|-----------|------|------|-----|---------|-----|
| S1        | 1    | 0    | 0   | 0       |     |
| S2        | 0    | 1    | 1   | 1       |     |

**建议**：包含非同义突变、移码插入/缺失、无义突变、剪接位点突变。避免混入同义突变。

#### 2. 依赖矩阵 (`dependency_matrix.csv`)

连续值矩阵（样本 × 基因），越负表示依赖越强（标准 CERES/Chronos 基因效应分数）。

| sample_id | EGFR  | ATM   | PARP1 | CDK6  | ... |
|-----------|-------|-------|-------|-------|-----|
| S1        | -0.45 | -0.12 | -0.78 | -0.20 |     |
| S2        | -0.10 | -0.35 | -0.15 | -0.60 |     |

#### 3. 样本元数据 (`sample_metadata.csv`)

| sample_id | lineage | TMB  | MSI_status | CNV_burden | batch |
|-----------|---------|------|------------|------------|-------|
| S1        | LUAD    | 5.2  | 0          | 0.12       | 1     |
| S2        | BRCA    | 3.1  | 1          | 0.08       | 2     |

**必含列**：`sample_id`, `lineage`, `TMB`, `MSI_status`, `CNV_burden`

**谱系命名建议**：使用 TCGA 肿瘤类型缩写（如 LUAD、BRCA、COAD、GBM 等）。

#### 4. STRING 边文件 (`STRING_edges.csv`)

| geneA | geneB | combined_score |
|-------|-------|----------------|
| TP53  | MDM2  | 999            |
| BRCA1 | BARD1 | 998            |

- `combined_score`：0–1000（STRING 综合评分）
- **重要**：基因符号必须与突变/依赖矩阵一致
- **注意**：如果使用 STRING 蛋白 ID（如 `9606.ENSP...`），须先映射为基因符号

#### 5. Reactome 通路文件 (`Reactome_gene_pathway.csv`)

| gene   | pathway                                |
|--------|----------------------------------------|
| BRCA1  | DNA Double-Strand Break Repair         |
| BARD1  | DNA Double-Strand Break Repair         |

#### 6. SynLethDB 合成致死对 (`SynLethDB_pairs.csv`)

| geneA | geneB  | evidence_type |
|-------|--------|---------------|
| BRCA1 | PARP1  | literature    |
| KRAS  | GATA2  | database      |

#### 7. 可药性靶点列表 (`druggable_targets.csv`)

| gene  | drug_category  | drug_name      |
|-------|----------------|-----------------|
| EGFR  | FDA_approved   | Osimertinib     |
| PARP1 | FDA_approved   | Olaparib        |
| CDK6  | clinical       | Palbociclib     |

**分类体系**：`FDA_approved` → `clinical` → `preclinical` → `druggable_class` → `unknown`

### 可选文件

#### 8. 表达矩阵 (`expression_matrix.csv`)

样本 × 基因。用作回归协变量，控制靶基因表达对依赖性的混杂影响。

#### 9. OmniPath 信号边 (`OmniPath_edges.csv`)

| source_gene | target_gene | effect      | confidence |
|-------------|-------------|-------------|------------|
| EGFR        | KRAS        | activation  | 0.95       |
| TP53        | CDKN1A      | activation  | 0.99       |

#### 10. CORUM 复合物 (`CORUM_complexes.csv`)

| gene  | complex_id |
|-------|------------|
| BRCA1 | COMPLEX_1  |
| BARD1 | COMPLEX_1  |
| SMARCA4 | COMPLEX_2 |

---

## 管线工作流

整个管线包含 9 个顺序步骤。每步将中间结果保存到 `data/processed/`，支持断点续跑。

| 步骤 | 模块 | 功能 | 输出 |
|------|------|------|------|
| **0** | `simulate_data.R` | 生成合成测试数据 | 10 个 CSV 文件到 `data/raw/` |
| **1** | `01_prepare_data.R` | 加载、验证、标准化基因符号 | `data_list`（内存对象） |
| **2** | `02_association_analysis.R` | 所有突变×依赖对的 LM 回归，BH-FDR 校正，候选过滤 | `association_results.csv`, `candidate_pairs.csv` |
| **3** | `03_build_networks.R` | 构建 STRING igraph、信号 igraph，通路/复合物/SL/药物查找索引 | `network_objects`（内存对象） |
| **4** | `04_compute_network_features.R` | STRING 最短路径、路径置信度、Hub 惩罚 | `network_features.csv` |
| **5** | `05_compute_empirical_evidence.R` | 通路一致性、方向信号、复合物、SL 先验、药物评分 | `empirical_features.csv` |
| **6** | `06_score_candidates.R` | 加权评分整合，证据分类 (I–IV)，排序 | `final_ranked_candidates.csv` |
| **7** | `07_null_model_permutation.R` | 度匹配置换检验，网络显著性 | `permutation_results.csv` |
| **8** | `08_visualization.R` | 点图、热图、网络路径图、分数分布 | PDF/PNG 图表 |
| **9** | `09_generate_report.R` | 候选机制卡片，HTML + Markdown 报告 | 报告文件 |

### 步骤依赖图

```
步骤 0 (模拟) ──> 步骤 1 (准备) ──> 步骤 2 (关联) ──> 步骤 3 (网络)
                     │                                      │
                     │              ┌────────────────────────┘
                     │              │
                     └──> 步骤 4 (网络特征) ──> 步骤 5 (经验证据) ──> 步骤 6 (评分)
                                                                       │
                                                            ┌──────────┼──────────┐
                                                            │          │          │
                                                       步骤 7      步骤 8      步骤 9
                                                      (置换检验)   (可视化)   (报告)
```

---

## 命令行用法详解

### 基本用法

```bash
# 一键运行：生成模拟数据 + 完整管线
Rscript R/run_pipeline.R --simulate --verbose

# 仅生成模拟数据，不运行分析
Rscript R/run_pipeline.R --simulate --stop-at 1

# 使用已有数据运行核心管线（跳过耗时的置换检验和图表）
Rscript R/run_pipeline.R --skip-permutation --skip-visualization --skip-report

# 从指定步骤恢复运行
Rscript R/run_pipeline.R --start-at 6              # 从评分步骤开始
Rscript R/run_pipeline.R --start-at 8 --stop-at 9  # 仅可视化 + 报告

# 使用自定义配置文件
Rscript R/run_pipeline.R --config my_custom_config.yaml
```

### 全部命令行选项

| 选项 | 说明 |
|------|------|
| `--config PATH` | YAML 配置文件路径（默认：`config/config.yaml`） |
| `--simulate` | 运行前先生成模拟数据 |
| `--skip-association` | 跳过步骤 2（使用已有关联结果） |
| `--skip-permutation` | 跳过步骤 7（置换检验耗时较长） |
| `--skip-visualization` | 跳过步骤 8 |
| `--skip-report` | 跳过步骤 9 |
| `--start-at N` | 从步骤 N 开始（1–9） |
| `--stop-at N` | 在步骤 N 停止（1–9） |
| `--verbose` | 详细日志输出 |

### 常用场景

```bash
# 场景 1：首次测试，快速验证管线是否能跑通
Rscript R/run_pipeline.R --simulate --skip-permutation --verbose

# 场景 2：已有真实数据，运行完整管线
Rscript R/run_pipeline.R --verbose

# 场景 3：修改了评分权重，从评分步骤恢复
Rscript R/run_pipeline.R --start-at 6

# 场景 4：仅重新生成图表和报告
Rscript R/run_pipeline.R --start-at 8

# 场景 5：关联分析太慢，跳过它（使用之前的缓存结果）
Rscript R/run_pipeline.R --skip-association

# 场景 6：调整了 config.yaml，重新运行全部
Rscript R/run_pipeline.R --simulate --verbose
```

---

## 输出结果说明

### 主要输出：`data/processed/final_ranked_candidates.csv`

按 `final_SL_score` 降序排列的候选对，包含完整证据注释。

#### 完整列参考

| 列名 | 类型 | 说明 |
|------|------|------|
| `rank` | int | 最终排名（1 = 最佳候选） |
| `mutation_gene` | chr | 突变基因 |
| `target_gene` | chr | 依赖靶基因 |
| `beta` | num | 突变回归系数 |
| `FDR` | num | Benjamini-Hochberg 校正 p 值 |
| `delta_dependency` | num | 平均依赖差值（突变组 − 野生型组） |
| `n_mut` | int | 突变样本数 |
| `stat_score` | num | 缩放统计分 (0–1) |
| `stat_score_raw` | num | 原始分：−log₁₀(FDR) × |beta| |
| `shortest_distance` | int | STRING 最短路径边数 |
| `path_genes` | chr | 路径基因（分号分隔） |
| `path_confidence` | num | 边置信度乘积 |
| `hub_penalty` | num | Hub 惩罚 = 1 / log₂(中间节点平均度 + 2) |
| `distance_score` | num | 距离映射分 (1→1.0, 2→0.7, ...) |
| `network_score` | num | distance_score × path_confidence × hub_penalty |
| `network_score_scaled` | num | Min-max 缩放至 0–1 |
| `pathway_score` | num | 0 / 0.5（父通路）/ 1（精确匹配） |
| `shared_pathways` | chr | 共享通路名（分号分隔） |
| `signaling_score` | num | 0 / 0.5（弱）/ 0.7（反向）/ 1（直接） |
| `signaling_direction` | chr | 调控方向说明 |
| `signaling_path` | chr | 有向路径基因 |
| `complex_score` | num | 1 = 相同复合物，否则 0 |
| `complex_evidence` | chr | 共享复合物 ID |
| `sl_prior_score` | num | 0 / 0.3（通路级）/ 0.5（邻居支持）/ 1（已知对） |
| `sl_prior_type` | chr | 证据层级说明 |
| `codep_r` | num | 依赖谱 Pearson 相关系数 |
| `codependency_score` | num | abs(codep_r)，0–1 |
| `drug_score` | num | 0 / 0.5 / 0.7 / 0.8 / 1 |
| `drug_category` | chr | 药物开发阶段 |
| `final_SL_score` | num | 平衡加权最终分 |
| `final_mechanism_score` | num | 机制导向加权分 |
| `final_therapeutic_score` | num | 治疗导向加权分 |
| `evidence_class` | chr | 证据类别 (I–IV 或 Low Priority) |
| `empirical_network_p` | num | 度匹配置换检验 p 值（步骤 7） |
| `observed_network_score` | num | 观察到的网络分 |
| `mean_random_score` | num | 随机对照平均网络分 |

### 中间文件

| 文件 | 说明 |
|------|------|
| `association_results.csv` | 所有测试对（过滤前），含 beta、p_value、FDR |
| `candidate_pairs.csv` | 通过统计阈值的候选对 |
| `network_features.csv` | 富集 STRING 网络特征的候选对 |
| `empirical_features.csv` | 富集全部经验证据的候选对 |
| `permutation_results.csv` | 含经验网络 p 值的候选对 |

### 报告和图表

| 路径 | 说明 |
|------|------|
| `results/figures/ranking_dotplot.pdf` | 候选排序点图（Top 30） |
| `results/figures/ranking_dotplot.png` | 同上，PNG 格式 |
| `results/figures/evidence_heatmap.pdf` | 证据热图 |
| `results/figures/score_distributions.pdf` | 各分数分布图 |
| `results/figures/network_path_*.pdf` | 网络路径子图（Top 候选） |
| `results/reports/synthetic_dependency_report.html` | HTML 格式完整报告 |
| `results/reports/synthetic_dependency_report.md` | Markdown 格式完整报告 |

---

## 配置文件参考

所有参数在 `config/config.yaml` 中统一管理。以下为关键配置段：

### 关联分析配置

```yaml
association:
  min_mut_count: 5          # 基因最少突变样本数
  fdr_method: "BH"          # 多重检验校正方法
  covariates:               # 回归协变量
    - "lineage"
    - "TMB"
    - "MSI_status"
    - "CNV_burden"
  use_expression_covariate: true
```

### 候选过滤

```yaml
candidate_filtering:
  fdr_threshold: 0.05       # FDR 阈值
  beta_direction: "negative" # 期望突变增强依赖 (beta < 0)
  min_n_mut: 5              # 最少突变样本数
  min_abs_delta_dependency: 0.2  # 最小效应量
```

### 网络构建

```yaml
networks:
  string:
    score_cutoff: 700       # STRING 综合评分阈值 (0–1000)
                            # 700 = 高置信度，900 = 极高置信度（更稀疏）
```

### 评分权重

```yaml
scoring:
  balanced:                 # 默认方案
    stat_weight: 0.35
    network_weight: 0.25
    pathway_weight: 0.15
    signaling_weight: 0.10
    complex_weight: 0.05
    sl_prior_weight: 0.05
    drug_weight: 0.05
  mechanism:                # 机制导向（侧重生物证据）
    stat_weight: 0.30
    network_weight: 0.30
    pathway_weight: 0.20
    signaling_weight: 0.10
    sl_prior_weight: 0.10
  therapeutic:              # 治疗导向（侧重可药性）
    stat_weight: 0.30
    drug_weight: 0.15
    ...
  default_scheme: "balanced"
```

### 证据分类阈值

```yaml
classification:
  class_I:
    fdr_max: 0.05
    network_distance_max: 2
    pathway_score_min: 0.5
    final_score_min: 0.75
```

---

## 评分体系详解

### 各组分计算方式

**统计分：**
```
stat_score_raw = −log₁₀(FDR + 10⁻³⁰⁰) × |beta|
stat_score = scale₀_₁(stat_score_raw)
```
（FDR=0 时用极小伪计数 10⁻³⁰⁰ 避免 Inf）

**网络分：**
```
distance_score ∈ {1.0, 0.7, 0.4, 0.1}  对应距离 {1, 2, 3, 4–5}
path_confidence = ∏(edge_confidences)  距离 ≤ 5 时
hub_penalty = 1 / log₂(mean_degree(中间基因) + 2)
network_score = distance_score × path_confidence × hub_penalty
```
Hub 惩罚防止通过 TP53、MYC、UBC 等高连接度基因的路径主导排名。

**通路分 (Reactome)：**
```
1.0 = 同一精确通路
0.5 = 同一父通路
0.0 = 无共享通路
```

**信号分 (OmniPath)：**
```
1.0 = 有向路径 A→B，长度 ≤ 3
0.7 = 有向路径 B→A，长度 ≤ 3
0.5 = 存在任意有向路径，长度 ≤ 5
0.0 = 无有向路径
```

**SL 先验分 (SynLethDB)：**
```
1.0 = 已知 SL 对（数据库中直接收录）
0.5 = 邻居支持（一个基因的 STRING 邻居与另一个基因是已知 SL）
0.3 = 通路级支持（两基因所在通路中有已知 SL 对）
0.0 = 无证据
```

**可药性分：**
```
1.0 = FDA 批准药物靶点
0.8 = 临床试验阶段靶点
0.7 = 临床前抑制剂靶点
0.5 = 可药类（激酶/酶/受体）
0.0 = 未知
```

---

## 证据分类

| 类别 | 标签 | 判定条件 | 解读建议 |
|------|------|----------|----------|
| **Class I** | 高置信度合成依赖候选 | FDR<0.05, beta<0, 距离≤2, 通路≥0.5, 总分≥0.75 | 最适合机制写作和实验验证 |
| **Class II** | 网络支持候选 | FDR<0.05, beta<0, 距离≤3, 总分≥0.60 | 统计+网络支持良好，适合假设生成 |
| **Class III** | SL 先验支持候选 | FDR<0.05, beta<0, SL先验分>0 | 有已知 SL 证据支持，适合验证/再发现 |
| **Class IV** | 纯统计候选 | FDR<0.05, beta<0, 网络分=0 | 统计显著但无网络证据——可能代表新生物学或假阳性 |
| **Low Priority** | 低优先级 | 不满足以上任何条件 | 需要更宽松阈值或额外证据 |

**重要**：Class IV 候选不会丢弃。它们可能代表尚未被现有数据库收录的真正新生物学（Bug #11）。

---

## 可视化

管线生成四类发表级图表：

### 1. 排序点图 (`ranking_dotplot.pdf`)
X 轴为 final_SL_score，Y 轴为 Top 30 候选对（mutation_gene → target_gene），颜色按证据类别着色，点大小按 −log₁₀(FDR) 缩放。

### 2. 证据热图 (`evidence_heatmap.pdf`)
行 = Top 候选对，列 = 各证据组分分。白（0）到深蓝（1）渐变。行注释标注证据类别，按分数相似性聚类。

### 3. 网络路径图 (`network_path_*.pdf`)
为每个有 STRING 路径的 Class I/II 候选对生成子图可视化：
- 突变基因（红色节点）
- 靶基因（蓝色节点）
- 中间基因（灰色节点）
- 边宽按 STRING 置信度缩放

### 4. 分数分布图 (`score_distributions.pdf`)
各组分分数和最终分数的密度分布图，分面排列检查分布形态。

---

## 模拟数据

`simulate_data.R` 生成含已知真实信号的合成测试数据：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `n_samples` | 500 | 细胞系/样本数 |
| `n_mutation_genes` | 200 | 突变基因数 |
| `n_dependency_genes` | 300 | 依赖靶基因数 |
| `n_true_sl_pairs` | 50 | 注入的真实阳性 SL 信号数 |
| `true_sl_effect_size` | −0.4 | 突变样本依赖值偏移量 |
| `string_coverage` | 0.85 | STRING 网络中的基因覆盖比例 |
| `n_pathway_terms` | 80 | Reactome 通路数 |
| `n_complexes` | 40 | CORUM 蛋白复合物数 |

**模拟数据关键特征**：
- 真实阳性 SL 对的依赖值在突变样本中向下偏移（增强依赖）
- 真实阳性对在 STRING 网络中通过 1 跳或 2 跳路径连接
- 真实阳性对共享 Reactome 通路，部分有向信号路径
- 真实阳性对的子集出现在 SynLethDB 对列表中
- 真实噪声：随机突变率、泛必需基因、谱系效应

真实阳性对已保存到 `data/raw/true_sl_pairs.csv`，可用来评估管线恢复率。

---

## 零模型与置换检验

网络邻近分对高连接度基因可能存在偏差（Hub 基因与一切都很接近）。为评估观察到的网络邻近是否显著强于随机期望，管线实现了度匹配置换检验。

**流程：**
1. 对每个候选对 (A, B)，计算观察到的 network_score
2. 保持突变基因 A 不变
3. 随机采样 1000 个与基因 B 的 STRING 度匹配的靶基因
4. 对每个置换计算 network_score(A, B_random)
5. 经验 p = (count(random ≥ observed) + 1) / (n_perm + 1)

**度匹配** 确保零分布考虑了高连接度基因普遍具有更短路径的事实。

---

## 生物学解读模板

对每个 Top 候选，管线生成结构化机制卡片：

```
候选: ARFGEF3 突变 → EGFR 依赖

1. 统计关联:
   ARFGEF3 突变样本对 EGFR 表现出更强的依赖性
   (beta = −0.38, FDR = 1e−3, delta_dependency = −0.34, n_mut = 10)

2. 网络支持:
   ARFGEF3 → 中间基因X → EGFR
   距离 = 2, 路径置信度 = 0.64, Hub 惩罚 = 0.82

3. 通路支持:
   共享通路: 膜转运 / 受体信号传导

4. SL 先验:
   已知 SL 证据为 无/邻居支持

5. 治疗相关性:
   EGFR 分类为 FDA_approved_target

6. 工作假设:
   ARFGEF3 突变可能重塑膜转运过程，增强细胞对 EGFR 信号的依赖。
   支持 EGFR 作为 ARFGEF3 突变背景下的候选合成依赖靶点。

7. 建议验证:
   - 比较 ARFGEF3 突变 vs 野生型模型中 EGFR 的依赖性差异
   - 在突变和野生型细胞中敲除或抑制 EGFR
   - 如有抑制剂，评估 ARFGEF3 突变细胞的 EGFR 抑制剂敏感性
```

### 用语规范

- 说 **"候选合成依赖"**，而非"已确认的合成致死"
- 说 **"处于已知信号级联中"**，而非"直接调控"（Bug #13）
- 说 **"网络邻近提示功能耦合"**，而非"物理相互作用"
- 说 **"可能增强细胞对...的依赖"**，而非"导致合成致死"

---

## 实验验证策略

### 计算验证

1. 在独立依赖数据集中检查关联一致性
2. 按肿瘤谱系分层分析（留一谱系分析）
3. 校正靶基因 CNV 和表达
4. 分别验证 CRISPR 和 RNAi 依赖数据
5. 与药物敏感性数据对比（CTRP、GDSC、PRISM）
6. 验证靶基因依赖是选择性的，非泛必需（Bug #18）
7. 与正交 SL 预测方法比较

### 湿实验验证设计

```
突变阳性细胞系             |     突变阴性匹配对照
         ↓                 |              ↓
   靶基因敲除/抑制         |    靶基因敲除/抑制
         ↓                 |              ↓
   细胞活力/克隆形成/凋亡  |    细胞活力/克隆形成/凋亡
         ↓                 |              ↓
         比较: IC50/AUC/活力差异
         ↓
   回复实验或通路生物标志物检测
```

**可药靶点**：比较抑制剂量效曲线 (IC50) 在突变阳性 vs 阴性细胞系中的差异。

**非可药靶点**：CRISPR 敲除或 siRNA 敲低后检测活力或增殖。

---

## 常见问题排查

| 症状 | 可能原因 | 解决方案 |
|------|----------|----------|
| 大多数基因不在网络中 | 基因符号不匹配（Bug #1） | 标准化为 HGNC 批准符号；检查大小写和别名 |
| STRING ID 类似 `9606.ENSP...` | 蛋白 ID 而非基因符号（Bug #2） | 通过 STRING 别名文件映射；使用 `map_protein_to_gene()` |
| Top 路径全经过 TP53/MYC/UBC | Hub 基因偏差（Bug #3） | Hub 惩罚默认启用；检查 hub_penalty 值 |
| 最短路径生物学意义不明 | 最短路径局限（Bug #4） | 提高 STRING 分数阈值至 900；人工检查路径基因 |
| 较长路径 path_confidence = 0 | 边分乘积趋近零（Bug #5） | 已限制距离≤5；长路径用 geometric_mean 方法 |
| −log₁₀(FDR) = Inf | FDR 数值为零（Bug #6） | safe_log10() 已使用伪计数 10⁻³⁰⁰ |
| scale_0_1 返回 NaN | 所有值相同（Bug #7） | 已处理：常数向量返回 0.5 |
| lm 系数为 NA | 模型矩阵奇异（Bug #8） | 提高 min_mut_count；检查协变量方差 |
| Top 对实际是谱系标记 | 谱系混杂（Bug #9） | 已在回归中包含 lineage；建议做谱系内分析 |
| 依赖被表达驱动 | 表达混杂（Bug #10） | 已将靶基因表达纳入协变量（默认启用） |
| 新对网络分始终为 0 | 网络过度惩罚新生物学（Bug #11） | Class IV 保护这些候选；人工检查 |
| 多处 NA 值 | 特征表合并后缺失（Bug #15） | 所有 NA 分数自动替换为 0 |
| 对数太多 | 组合爆炸（Bug #20） | 按突变频率和依赖方差预过滤基因 |
| 推送超时 | 网络限制 | 在本地终端上执行 git push |

---

## 已知局限

1. **网络证据是生物学合理性的参考，不是证明。** Class IV 候选可能代表尚未被数据库收录的新生物学。
2. **STRING 包含功能关联而不仅是物理 PPI。** STRING 边源自共表达、文本挖掘、基因组上下文等多种证据。
3. **最短路径未必是生物学最相关的路径。** 路径最小化算法可能选择 Hub 介导路径。Hub 惩罚可缓解但不能消除。
4. **泛癌分析可能混淆谱系特异性效应。** 一对可能在泛癌中显著只因突变集中于某癌种。务必进行谱系内敏感性分析。
5. **已知 SL 数据库不完整且有研究偏差。** 新的 SL 对可能即使真实也缺乏数据库支持。
6. **共依赖不等于合成致死。** 强相关依赖谱可能反映共享必需性、复合物共成员或技术假象。
7. **有向信号路径不证明因果性。** OmniPath 有向边表示文献证据支持调控关系，不一定是直接生化因果关系。
8. **可药性分数反映当前知识。** 今天被归类为 "unknown" 的靶点未来可能变为可药。治疗导向评分方案可上调可药性权重。

---

## 创新点定位

本框架与现有方法的关键区别：

### 1. 突变条件化的配对发现

输出不是"EGFR 是依赖基因"，而是 **"ARFGEF3 突变肿瘤可能对 EGFR 抑制更敏感"**。这直接对齐精准肿瘤学的生物标志物策略。

### 2. 多层经验证据整合

不依赖单一网络或数据库，而是整合 7 个独立证据层，每层衡量生物学合理性的不同维度。这形成了透明、可审计的排序框架。

### 3. 可解释的生物学路径

每个候选配有一个明确的机制路径（mutation_gene → 中间基因 → target_gene），可直接支持机制图和实验假设。

### 4. 混杂感知的统计框架

回归模型校正了谱系、TMB、MSI、CNV 负荷和靶基因表达——CRISPR 依赖分析中已知的主要混杂因素。

### 5. 度匹配零模型

置换检验控制 Hub 基因与一切接近的事实，提供经验显著性估计，而非依赖任意网络距离阈值。

### 6. 治疗优先排序

可药性评分支持按转化潜力分层，提供独立的机制导向和治疗导向评分方案。

---

## 参考文献

- **STRING:** Szklarczyk D, et al. The STRING database in 2023. *Nucleic Acids Research*, 2023.
- **Reactome:** Gillespie M, et al. The reactome pathway knowledgebase 2022. *Nucleic Acids Research*, 2022.
- **OmniPath:** Turei D, et al. Integrated intra- and intercellular signaling knowledge. *Molecular Systems Biology*, 2021.
- **CORUM:** Giurgiu M, et al. CORUM: the comprehensive resource of mammalian protein complexes. *Nucleic Acids Research*, 2019.
- **SynLethDB:** Guo J, et al. SynLethDB 2.0: a web-based knowledge graph database on synthetic lethality. *Nucleic Acids Research*, 2023.
- **DepMap:** Tsherniak A, et al. Defining a cancer dependency map. *Cell*, 2017.
- **DGIdb:** Freshour SL, et al. Integration of the Drug-Gene Interaction Database. *Nucleic Acids Research*, 2021.
- **Open Targets:** Ochoa D, et al. Open Targets Platform. *Nucleic Acids Research*, 2021.

---

## 许可证

本项目基于 MIT 许可证开源。

---

*完整规范、方法细节和设计理由见 [project.md](project.md)。英文 README 见 [README.md](README.md)。*
