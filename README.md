# Empirical Network-Constrained Synthetic Dependency Prioritization

[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A systematic bioinformatics framework for prioritizing candidate **synthetic lethal / synthetic dependency gene pairs** from large-scale mutation and CRISPR dependency data by integrating statistical association with multi-layer empirical biological network evidence.

---

## Table of Contents

1. [Biological Motivation](#biological-motivation)
2. [Methodology Overview](#methodology-overview)
3. [Project Structure](#project-structure)
4. [Quick Start](#quick-start)
5. [Installation & Dependencies](#installation--dependencies)
6. [Input Data](#input-data)
7. [Pipeline Workflow](#pipeline-workflow)
8. [Output & Results](#output--results)
9. [Configuration Reference](#configuration-reference)
10. [Scoring System](#scoring-system)
11. [Evidence Classification](#evidence-classification)
12. [Visualization](#visualization)
13. [Simulated Data](#simulated-data)
14. [Null Model & Permutation Testing](#null-model--permutation-testing)
15. [CLI Usage](#cli-usage)
16. [Biological Interpretation](#biological-interpretation)
17. [Experimental Validation Guide](#experimental-validation-guide)
18. [Troubleshooting](#troubleshooting)
19. [Known Limitations](#known-limitations)
20. [Novelty & Positioning](#novelty--positioning)
21. [References](#references)
22. [License](#license)

---

## Biological Motivation

Cancer genomes harbor thousands of mutations, but only a subset create **conditional vulnerabilities** — situations where tumor cells become more dependent on specific genes for survival. Identifying these mutation-conditioned dependencies is a central challenge in precision oncology.

The core biological question this framework addresses is:

> **If mutation of gene A is associated with increased dependency on gene B, can empirical biological networks help determine whether this association is mechanistically plausible and therapeutically actionable?**

Pure statistical associations between mutations and CRISPR dependency can be misleading due to:

- Cancer lineage effects and batch confounding
- Tumor mutation burden (TMB) and genomic instability
- Microsatellite instability (MSI) status
- Copy number burden artifacts
- Hub gene artifacts in biological networks
- Low mutation frequency (small sample sizes)
- Passenger mutation effects

This framework addresses these challenges by layering **empirical biological network evidence** on top of statistical association, producing interpretable, mechanistically-grounded candidate rankings.

---

## Methodology Overview

The framework transforms raw mutation and dependency data into ranked, biologically-interpretable synthetic dependency candidates through nine sequential modules:

```
Mutation Matrix  +  CRISPR Dependency Matrix  +  Sample Metadata
         │                    │                        │
         └────────────────────┼────────────────────────┘
                              │
                   Regression-Based Association
                   (lineage, TMB, MSI, CNV adjusted)
                              │
                   Candidate Pairs (FDR < 0.05)
                              │
         ┌────────────────────┼────────────────────────┐
         │                    │                        │
   STRING Network       Reactome Pathway        OmniPath Signaling
   (proximity)          (concordance)          (directionality)
         │                    │                        │
         ├────────────────────┼────────────────────────┤
         │                    │                        │
   CORUM Complex        SynLethDB Prior         Drug Target Evidence
   (co-complex)         (known SL)              (druggability)
         │                    │                        │
         └────────────────────┼────────────────────────┘
                              │
                   Weighted Multi-Evidence Scoring
                              │
                   Ranked Synthetic Dependency Candidates
                              │
         ┌────────────────────┼────────────────────────┐
         │                    │                        │
   Mechanism Cards      Publication Figures      HTML/MD Report
```

### Evidence Layers

| Layer | Data Source | What It Measures | Score Range |
|-------|-------------|-----------------|-------------|
| **Statistical** | Association analysis | Mutation→dependency effect size × significance | 0–1 (scaled) |
| **Network proximity** | STRING v12 | Shortest path distance, confidence, hub penalty | 0–1 (composite) |
| **Pathway concordance** | Reactome | Shared biological pathways | 0 / 0.5 / 1 |
| **Directional signaling** | OmniPath / SIGNOR | Directed regulatory path A→B or B→A | 0 / 0.5 / 0.7 / 1 |
| **Protein complex** | CORUM | Shared protein complex membership | 0 / 1 |
| **SL prior** | SynLethDB | Known or neighbor-supported SL evidence | 0 / 0.3 / 0.5 / 1 |
| **Co-dependency** | DepMap | Pearson correlation of dependency profiles | 0–1 (abs(r)) |
| **Druggability** | DGIdb / Open Targets | Target gene therapeutic tractability | 0 / 0.5 / 0.7 / 0.8 / 1 |

### Final Scoring Formula (Balanced)

```
Final_SL_Score = 0.35 × stat_score
               + 0.25 × network_score
               + 0.15 × pathway_score
               + 0.10 × signaling_score
               + 0.05 × complex_score
               + 0.05 × sl_prior_score
               + 0.05 × drug_score
```

Three scoring schemes are provided: **balanced**, **mechanism-focused** (emphasizes network/pathway), and **therapeutic-focused** (emphasizes druggability).

### Regression Model

For each mutation gene A × dependency gene B pair:

```
dependency_B ~ mutation_A + expression_B + lineage + TMB + MSI_status + CNV_burden
```

This controls for the major known confounders in CRISPR dependency data (lineage effects, target expression, genomic instability).

---

## Project Structure

```
network_synthetic_dependency_project/
│
├── config/
│   └── config.yaml                    # Single source of truth: all parameters
│
├── R/
│   ├── 00_utils.R                     # Shared utilities (scaling, safe paths, logging)
│   ├── 01_prepare_data.R              # Data loading, validation, symbol standardization
│   ├── 02_association_analysis.R      # LM regression, FDR correction, filtering
│   ├── 03_build_networks.R            # STRING/signaling graph construction, lookup indices
│   ├── 04_compute_network_features.R  # Shortest path, path confidence, hub penalty
│   ├── 05_compute_empirical_evidence.R # Pathway, signaling, complex, SL prior, drug scores
│   ├── 06_score_candidates.R          # Weighted integration, classification (I–IV)
│   ├── 07_null_model_permutation.R    # Degree-matched empirical p-values
│   ├── 08_visualization.R             # Dot plots, heatmaps, network path figures
│   ├── 09_generate_report.R           # Mechanism cards, HTML/Markdown reports
│   ├── simulate_data.R                # Realistic synthetic test data generator
│   └── run_pipeline.R                 # CLI orchestrator with full step control
│
├── data/
│   ├── raw/                           # Input data (populated by user or simulate_data.R)
│   │   ├── mutation_matrix.csv
│   │   ├── dependency_matrix.csv
│   │   ├── expression_matrix.csv
│   │   ├── sample_metadata.csv
│   │   ├── STRING_edges.csv
│   │   ├── Reactome_gene_pathway.csv
│   │   ├── OmniPath_edges.csv
│   │   ├── CORUM_complexes.csv
│   │   ├── SynLethDB_pairs.csv
│   │   └── druggable_targets.csv
│   │
│   └── processed/                     # Pipeline outputs (auto-generated)
│       ├── association_results.csv
│       ├── candidate_pairs.csv
│       ├── network_features.csv
│       ├── empirical_features.csv
│       ├── final_ranked_candidates.csv
│       └── permutation_results.csv
│
├── results/
│   ├── tables/                        # Final ranked tables
│   ├── figures/                       # Publication-quality PDF/PNG figures
│   │   ├── ranking_dotplot.pdf
│   │   ├── ranking_dotplot.png
│   │   ├── evidence_heatmap.pdf
│   │   ├── score_distributions.pdf
│   │   └── network_path_*.pdf
│   └── reports/                       # Auto-generated reports
│       ├── synthetic_dependency_report.html
│       └── synthetic_dependency_report.md
│
├── logs/
│   └── run_log.txt
│
├── README.md
└── project.md                         # Full specification document
```

---

## Quick Start

```bash
# Clone and enter project
cd network_synthetic_dependency_project

# Generate simulated test data (500 samples, 200 mut genes, 300 dep genes, 50 true SL pairs)
Rscript R/run_pipeline.R --simulate --stop-at 1

# Run the full pipeline (association → networks → features → scoring → report)
Rscript R/run_pipeline.R --verbose

# Or do everything in one command:
Rscript R/run_pipeline.R --simulate --verbose

# View top candidates
head data/processed/final_ranked_candidates.csv

# Open the auto-generated report
open results/reports/synthetic_dependency_report.html
```

---

## Installation & Dependencies

### System Requirements

- **R** >= 4.0
- **OS**: Linux, macOS, or Windows
- **Memory**: >= 8 GB recommended for genome-scale analysis
- **Disk**: ~500 MB for intermediate files (depends on input size)

### R Packages

**Required:**
```r
install.packages(c("dplyr", "ggplot2", "tidyr", "purrr", "readr", "igraph", "yaml", "optparse"))
```

**Optional but recommended:**
```r
install.packages(c("pheatmap"))   # Enhanced heatmaps with row annotations
install.packages(c("furrr"))       # Parallel processing for association step
```

- `igraph` — network analysis (graph construction, shortest paths, centrality)
- `yaml` — configuration file parsing
- `optparse` — command-line argument parsing
- `dplyr`, `tidyr`, `purrr`, `readr` — data manipulation
- `ggplot2` — visualization

---

## Input Data

### Required Files

#### 1. Mutation Matrix (`mutation_matrix.csv`)

Binary matrix (samples × genes). 1 = nonsynonymous mutation present, 0 = absent.

| sample_id | TP53 | KRAS | ATM | ARFGEF3 | ... |
|-----------|------|------|-----|---------|-----|
| S1        | 1    | 0    | 0   | 0       |     |
| S2        | 0    | 1    | 1   | 1       |     |

**Recommendations:**
- Include nonsynonymous, frameshift, nonsense, and splice-site mutations
- Avoid silent mutations unless specifically needed
- Loss-of-function-only matrices can be used for stricter analysis

#### 2. Dependency Matrix (`dependency_matrix.csv`)

Continuous matrix (samples × genes). More negative values = stronger dependency (standard CERES/Chronos gene effect scores).

| sample_id | EGFR  | ATM   | PARP1 | CDK6  | ... |
|-----------|-------|-------|-------|-------|-----|
| S1        | -0.45 | -0.12 | -0.78 | -0.20 |     |
| S2        | -0.10 | -0.35 | -0.15 | -0.60 |     |

#### 3. Sample Metadata (`sample_metadata.csv`)

| sample_id | lineage | TMB  | MSI_status | CNV_burden | batch |
|-----------|---------|------|------------|------------|-------|
| S1        | LUAD    | 5.2  | 0          | 0.12       | 1     |
| S2        | BRCA    | 3.1  | 1          | 0.08       | 2     |

**Required columns:** `sample_id`, `lineage`, `TMB`, `MSI_status`, `CNV_burden`

#### 4. STRING Edges (`STRING_edges.csv`)

Functional association network edges with confidence scores.

| geneA | geneB | combined_score |
|-------|-------|---------------|
| TP53  | MDM2  | 999           |
| BRCA1 | BARD1 | 998           |

- `combined_score`: 0–1000 (STRING combined score)
- **Important:** Gene symbols must match mutation/dependency matrices (Bug #1)
- **Important:** If using STRING protein IDs (e.g., `9606.ENSP...`), map to gene symbols first (Bug #2)

#### 5. Reactome Pathways (`Reactome_gene_pathway.csv`)

| gene   | pathway                                |
|--------|----------------------------------------|
| BRCA1  | DNA Double-Strand Break Repair         |
| BARD1  | DNA Double-Strand Break Repair         |

#### 6. SynLethDB Pairs (`SynLethDB_pairs.csv`)

Known synthetic lethal gene pairs.

| geneA | geneB  | evidence_type |
|-------|--------|---------------|
| BRCA1 | PARP1  | literature    |
| KRAS  | GATA2  | database      |

#### 7. Druggable Targets (`druggable_targets.csv`)

| gene  | drug_category  | drug_name      |
|-------|----------------|-----------------|
| EGFR  | FDA_approved   | Osimertinib     |
| PARP1 | FDA_approved   | Olaparib        |
| CDK6  | clinical       | Palbociclib     |

**Categories:** `FDA_approved`, `clinical`, `preclinical`, `druggable_class`, `unknown`

### Optional Files

#### 8. Expression Matrix (`expression_matrix.csv`)

Samples × genes. Used as covariate in regression to control for target gene expression confounding (Bug #10).

#### 9. OmniPath Signaling Edges (`OmniPath_edges.csv`)

| source_gene | target_gene | effect      | confidence |
|-------------|-------------|-------------|------------|
| EGFR        | KRAS        | activation  | 0.95       |

#### 10. CORUM Complexes (`CORUM_complexes.csv`)

| gene  | complex_id |
|-------|------------|
| BRCA1 | COMPLEX_1  |
| BARD1 | COMPLEX_1  |

---

## Pipeline Workflow

The pipeline consists of 9 sequential steps. Each step saves intermediate results to `data/processed/`, enabling checkpoint-based resume.

| Step | Module | Description | Output |
|------|--------|-------------|--------|
| **0** | `simulate_data.R` | Generate synthetic test data | 10 CSV files in `data/raw/` |
| **1** | `01_prepare_data.R` | Load, validate, standardize gene symbols | `data_list` (in-memory) |
| **2** | `02_association_analysis.R` | LM regression on all mutation×dependency pairs, BH-FDR correction, candidate filtering | `association_results.csv`, `candidate_pairs.csv` |
| **3** | `03_build_networks.R` | Build STRING igraph, signaling igraph, pathway/complex/SL/drug lookup indices | `network_objects` (in-memory) |
| **4** | `04_compute_network_features.R` | STRING shortest path, path confidence, hub penalty for each candidate | `network_features.csv` |
| **5** | `05_compute_empirical_evidence.R` | Pathway concordance, directional signaling, complex, SL prior, drug scores | `empirical_features.csv` |
| **6** | `06_score_candidates.R` | Weighted score integration, evidence classification (I–IV), ranking | `final_ranked_candidates.csv` |
| **7** | `07_null_model_permutation.R` | Degree-matched permutation test for network significance | `permutation_results.csv` |
| **8** | `08_visualization.R` | Dot plot, heatmap, network path plots, score distributions | PDF/PNG figures |
| **9** | `09_generate_report.R` | Mechanism cards for top candidates, HTML + Markdown report | Report files |

### Step Dependency Graph

```
Step 0 (simulate) ──> Step 1 (prepare) ──> Step 2 (association) ──> Step 3 (networks)
                         │                                              │
                         │                    ┌─────────────────────────┘
                         │                    │
                         └──> Step 4 (network features) ──> Step 5 (empirical) ──> Step 6 (scoring)
                                                                                     │
                                                                          ┌──────────┼──────────┐
                                                                          │          │          │
                                                                     Step 7    Step 8      Step 9
                                                                   (permutation) (viz)    (report)
```

---

## Output & Results

### Main Output: `final_ranked_candidates.csv`

Contains all candidate pairs ranked by `final_SL_score` with complete evidence annotations.

#### Full Column Reference

| Column | Type | Description |
|--------|------|-------------|
| `rank` | int | Final ranking (1 = top candidate) |
| `mutation_gene` | chr | Mutated gene |
| `target_gene` | chr | Dependency target gene |
| `beta` | num | Regression coefficient for mutation |
| `FDR` | num | Benjamini-Hochberg adjusted p-value |
| `delta_dependency` | num | Mean dependency difference (mutant − WT) |
| `n_mut` | int | Number of mutant samples |
| `stat_score` | num | Scaled statistical score (0–1) |
| `stat_score_raw` | num | Raw: −log₁₀(FDR) × |beta| |
| `shortest_distance` | int | STRING shortest path edge count |
| `path_genes` | chr | Path genes (semicolon-separated) |
| `path_confidence` | num | Product of edge confidence scores |
| `hub_penalty` | num | 1 / log₂(mean intermediate degree + 2) |
| `distance_score` | num | Distance → score mapping (1→1.0, 2→0.7, ...) |
| `network_score` | num | distance_score × path_confidence × hub_penalty |
| `network_score_scaled` | num | Min-max scaled to 0–1 |
| `pathway_score` | num | 0, 0.5 (parent), or 1 (exact match) |
| `shared_pathways` | chr | Semicolon-separated shared pathway names |
| `signaling_score` | num | 0, 0.5 (weak), 0.7 (reverse), or 1 (direct) |
| `signaling_direction` | chr | Direction of regulatory path |
| `signaling_path` | chr | Directed path genes |
| `complex_score` | num | 1 if same complex, else 0 |
| `complex_evidence` | chr | Shared complex ID(s) |
| `sl_prior_score` | num | 0, 0.3 (pathway), 0.5 (neighbor), or 1 (known pair) |
| `sl_prior_type` | chr | Evidence tier description |
| `codep_r` | num | Pearson correlation of dependency profiles |
| `codependency_score` | num | abs(codep_r), 0–1 |
| `drug_score` | num | 0, 0.5 (druggable class), 0.7 (preclinical), 0.8 (clinical), 1 (FDA) |
| `drug_category` | chr | Drug development stage |
| `final_SL_score` | num | Balanced weighted score |
| `final_mechanism_score` | num | Mechanism-focused weighted score |
| `final_therapeutic_score` | num | Therapeutic-focused weighted score |
| `evidence_class` | chr | Class I–IV or Low Priority |
| `empirical_network_p` | num | Degree-matched permutation p-value (Step 7) |

### Intermediate Files

| File | Description |
|------|-------------|
| `association_results.csv` | All tested pairs (before filtering) with beta, p-value, FDR |
| `candidate_pairs.csv` | Filtered candidates passing statistical thresholds |
| `network_features.csv` | Candidates enriched with STRING network features |
| `empirical_features.csv` | Candidates enriched with all empirical evidence |
| `permutation_results.csv` | Top candidates with empirical network p-values |

---

## Configuration Reference

All parameters are in `config/config.yaml`. Key sections:

### Association Analysis
```yaml
association:
  min_mut_count: 5          # Minimum mutant samples per gene
  fdr_method: "BH"          # Multiple testing correction method
  covariates:               # Regression covariates
    - "lineage"
    - "TMB"
    - "MSI_status"
    - "CNV_burden"
```

### Candidate Filtering
```yaml
candidate_filtering:
  fdr_threshold: 0.05       # FDR cutoff
  beta_direction: "negative" # Expect mutation to increase dependency
  min_n_mut: 5              # Minimum mutant count
  min_abs_delta_dependency: 0.2  # Effect size threshold
```

### Network Construction
```yaml
networks:
  string:
    score_cutoff: 700       # STRING combined score threshold (0–1000)
                            # 700 = high confidence, 900 = very high
```

### Scoring Weights
```yaml
scoring:
  balanced:                 # Default scheme
    stat_weight: 0.35
    network_weight: 0.25
    pathway_weight: 0.15
    signaling_weight: 0.10
    complex_weight: 0.05
    sl_prior_weight: 0.05
    drug_weight: 0.05
  mechanism:                # Emphasizes biological evidence
    stat_weight: 0.30
    network_weight: 0.30
    pathway_weight: 0.20
    signaling_weight: 0.10
    sl_prior_weight: 0.10
  therapeutic:              # Emphasizes translational potential
    stat_weight: 0.30
    drug_weight: 0.15
    ...
```

### Evidence Classification Thresholds
```yaml
classification:
  class_I:
    fdr_max: 0.05
    network_distance_max: 2
    pathway_score_min: 0.5
    final_score_min: 0.75
```

### Permutation Testing
```yaml
permutation:
  enabled: true
  n_perm: 1000              # Number of permutations per pair
  top_n_pairs: 200          # Only test top pairs (expensive)
  match_on: ["degree"]      # Degree-matched null model
```

---

## Scoring System

### Component Scoring Details

**Statistical Score:**
```
stat_score_raw = −log₁₀(FDR + 10⁻³⁰⁰) × |beta|
stat_score = scale₀_₁(stat_score_raw)
```

**Network Score:**
```
distance_score ∈ {1.0, 0.7, 0.4, 0.1} for distances {1, 2, 3, 4–5}
path_confidence = ∏(edge_confidences) for paths ≤ 5 (capped)
hub_penalty = 1 / log₂(mean_degree(intermediate_genes) + 2)
network_score = distance_score × path_confidence × hub_penalty
```

The hub penalty prevents paths through highly-connected genes (TP53, MYC, UBC, etc.) from dominating the ranking (Bug #3).

**Pathway Score (Reactome):**
```
1.0 = same exact pathway
0.5 = same parent pathway
0.0 = no shared pathway
```

**Signaling Score (OmniPath):**
```
1.0 = directed path A→B, length ≤ 3
0.7 = directed path B→A, length ≤ 3
0.5 = any directed path, length ≤ 5
0.0 = no directed path
```

**SL Prior Score (SynLethDB):**
```
1.0 = known SL pair in database
0.5 = neighbor-supported (one gene's STRING neighbor is SL with the other)
0.3 = pathway-level support (genes share pathway with known SL pair)
0.0 = no evidence
```

**Drug Score:**
```
1.0 = FDA-approved drug target
0.8 = clinical trial target
0.7 = preclinical inhibitor target
0.5 = druggable class (kinase, enzyme, receptor)
0.0 = unknown
```

---

## Evidence Classification

| Class | Label | Criteria | Interpretation |
|-------|-------|----------|---------------|
| **Class I** | High-confidence synthetic dependency | FDR < 0.05, beta < 0, distance ≤ 2, pathway ≥ 0.5, score ≥ 0.75 | Best candidates for mechanism writing and experimental validation |
| **Class II** | Network-supported candidate | FDR < 0.05, beta < 0, distance ≤ 3, score ≥ 0.60 | Good statistical + network support, useful for hypothesis generation |
| **Class III** | SL-prior-supported candidate | FDR < 0.05, beta < 0, sl_prior_score > 0 | Supported by known SL evidence, useful for validation/rediscovery |
| **Class IV** | Statistical-only candidate | FDR < 0.05, beta < 0, network_score = 0 | Statistically significant but no network evidence — may represent novel biology or false positives |
| **Low Priority** | — | Fails all class criteria | Requires weaker thresholds or additional evidence |

**Important:** Class IV candidates are NOT discarded. They may represent genuinely novel biology not yet captured in existing databases (Bug #11).

---

## Visualization

The pipeline generates four types of publication-quality figures:

### 1. Ranking Dot Plot
Shows top 30 candidates with final SL score on x-axis, colored by evidence class, sized by −log₁₀(FDR).

### 2. Evidence Heatmap
Rows = top candidates, columns = component scores. Gradient from white (0) to dark blue (1). Row annotations indicate evidence class. Clustered by score similarity.

### 3. Network Path Plots
For each Class I/II candidate with a STRING path, a subgraph visualization shows:
- Mutation gene (red node)
- Target gene (blue node)
- Intermediate genes (gray nodes)
- Edge width proportional to STRING confidence

### 4. Score Distributions
Density plots for each component score and final score, faceted for distribution inspection.

---

## Simulated Data

`simulate_data.R` generates realistic synthetic test data with known ground truth:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_samples` | 500 | Number of cell lines / samples |
| `n_mutation_genes` | 200 | Number of mutated genes |
| `n_dependency_genes` | 300 | Number of dependency target genes |
| `n_true_sl_pairs` | 50 | Injected true positive SL signals |
| `true_sl_effect_size` | −0.4 | Dependency shift in mutant samples |
| `string_coverage` | 0.85 | Fraction of genes present in STRING |
| `n_pathway_terms` | 80 | Number of Reactome pathways |
| `n_complexes` | 40 | Number of CORUM protein complexes |

**Key features of simulated data:**
- True positive SL pairs have their dependency values shifted downward in mutant samples
- True pairs are connected by 1-hop or 2-hop paths in the STRING network
- True pairs share Reactome pathways and some have directed signaling paths
- A subset of true pairs appear in the SynLethDB pair list
- Realistic noise: random mutation rates, pan-essential genes, lineage effects

The ground truth pairs are saved to `data/raw/true_sl_pairs.csv` for benchmarking pipeline recovery.

---

## Null Model & Permutation Testing

Network proximity scores can be inflated for high-degree genes (hubs are close to everything). To assess whether observed network proximity is stronger than expected by chance, the pipeline implements degree-matched permutation testing (Bug #3):

**Procedure:**
1. For each candidate pair (A, B), compute observed network score
2. Keep mutation gene A fixed
3. Sample 1000 random target genes matched by STRING degree to gene B
4. Compute network_score(A, B_random) for each permutation
5. Empirical p = (count(random ≥ observed) + 1) / (n_perm + 1)

**Degree matching** ensures the null distribution accounts for the fact that high-degree genes tend to have shorter paths to everything.

---

## CLI Usage

### Basic Commands

```bash
# Full pipeline with simulated data
Rscript R/run_pipeline.R --simulate --verbose

# Generate simulated data only (no analysis)
Rscript R/run_pipeline.R --simulate --stop-at 1

# Core pipeline without expensive steps
Rscript R/run_pipeline.R --skip-permutation --skip-visualization --skip-report

# Resume from a specific step
Rscript R/run_pipeline.R --start-at 6                # Resume from scoring
Rscript R/run_pipeline.R --start-at 8 --stop-at 9    # Just viz + report

# Use custom configuration
Rscript R/run_pipeline.R --config my_custom_config.yaml
```

### CLI Options

| Flag | Description |
|------|-------------|
| `--config PATH` | Path to YAML config (default: `config/config.yaml`) |
| `--simulate` | Generate simulated data before running |
| `--skip-association` | Skip Step 2 (use existing results) |
| `--skip-permutation` | Skip Step 7 (permutation is slow) |
| `--skip-visualization` | Skip Step 8 |
| `--skip-report` | Skip Step 9 |
| `--start-at N` | Start from step N (1–9) |
| `--stop-at N` | Stop after step N (1–9) |
| `--verbose` | Enable detailed logging |

---

## Biological Interpretation

For each top-ranked candidate, the pipeline generates a structured mechanism card:

```
Candidate: ARFGEF3 mutation → EGFR dependency

1. Statistical association:
   ARFGEF3-mutant samples showed stronger dependency on EGFR
   (beta = −0.38, FDR = 1e−3, delta_dependency = −0.34)

2. Network support:
   ARFGEF3 → intermediate_gene → EGFR
   distance = 2, path_confidence = 0.64

3. Pathway support:
   Shared pathway: membrane trafficking / receptor signaling

4. Synthetic lethality prior:
   Known SL evidence is absent / neighbor-supported

5. Therapeutic relevance:
   EGFR is classified as FDA_approved_target

6. Working hypothesis:
   Mutation of ARFGEF3 may rewire membrane trafficking, increasing
   cellular reliance on EGFR signaling. This supports EGFR as a
   candidate synthetic dependency in ARFGEF3-altered contexts.

7. Suggested validation:
   - Compare EGFR dependency in ARFGEF3-mutant vs wild-type models
   - Knockout or inhibit EGFR in mutant vs wild-type cells
   - Evaluate EGFR inhibitor sensitivity in ARFGEF3-mutant cells
```

### Language Guidelines

- Say **"candidate synthetic dependency"**, not "confirmed synthetic lethality"
- Say **"embedded in a known signaling cascade"**, not "directly regulates" (Bug #13)
- Say **"network proximity suggests functional coupling"**, not "physical interaction"

---

## Experimental Validation Guide

### Computational Validation

1. Check association in independent dependency datasets (e.g., Project DRIVE vs DepMap)
2. Stratify by cancer lineage (leave-one-lineage-out analysis)
3. Adjust for target gene CNV and expression
4. Validate in both CRISPR and RNAi dependency data separately
5. Compare with drug sensitivity data (CTRP, GDSC, PRISM)
6. Verify target dependency is selective, not pan-essential (Bug #18)
7. Compare with orthogonal SL prediction methods

### Wet-Lab Validation Design

```
mutation-positive cell line    |    mutation-negative matched control
         ↓                     |              ↓
   target_gene knockout        |    target_gene knockout
   or inhibitor treatment      |    or inhibitor treatment
         ↓                     |              ↓
   cell viability / colony     |    cell viability / colony
   formation / apoptosis       |    formation / apoptosis
         ↓                     |              ↓
   Compare: IC50 / AUC / viability difference between groups
         ↓
   Rescue experiment or pathway biomarker assay
```

**For druggable targets:** Compare inhibitor dose-response curves (IC50) between mutation-positive and mutation-negative cell lines.

**For non-druggable targets:** CRISPR knockout or siRNA knockdown followed by viability or proliferation assay.

---

## Troubleshooting

### Common Issues & Solutions

| Symptom | Likely Cause | Solution |
|---------|-------------|----------|
| Most genes not found in network | Gene symbol mismatch (Bug #1) | Standardize symbols to HGNC-approved names; check case, aliases |
| STRING IDs look like `9606.ENSP...` | Protein IDs, not gene symbols (Bug #2) | Map via STRING alias file; use `map_protein_to_gene()` in utils |
| Top paths all go through TP53, MYC, UBC | Hub gene bias (Bug #3) | Hub penalty is enabled by default; check hub_penalty values |
| Shortest path seems biologically irrelevant | Shortest path limitation (Bug #4) | Use higher STRING score cutoff (900); inspect path genes manually |
| Path confidence = 0 for long paths | Product of edge scores vanishes (Bug #5) | Capped at distance=5; use geometric_mean method for longer paths |
| −log₁₀(FDR) = Inf | FDR numerically zero (Bug #6) | Handled by safe_log10() with pseudocount 10⁻³⁰⁰ |
| scale_0_1 returns NaN | All values identical (Bug #7) | Returns 0.5 for constant vectors |
| lm coefficients are NA | Singular model matrix (Bug #8) | Increase min_mut_count; check covariate variance |
| Top pairs are lineage markers | Lineage confounding (Bug #9) | Lineage is included as covariate; perform within-lineage analysis |
| Dependency driven by expression | Expression confounding (Bug #10) | Include target expression as covariate; enabled by default |
| Network score always 0 for novel pairs | Network over-penalizes novel biology (Bug #11) | Class IV preserves these; inspect manually |
| Using same SL DB for scoring and validation | Circular validation (Bug #12) | SL prior is one component among many; external validation kept separate |
| NA values in score columns | Missing data after merging (Bug #15) | All NA scores replaced with 0 automatically |
| Too many pairwise tests | Combinatorial explosion (Bug #20) | Pre-filter genes by mutation frequency and dependency variance |

---

## Known Limitations

1. **Network evidence is biological plausibility, not proof.** Network proximity supports but does not confirm synthetic dependency. Class IV candidates may represent novel biology not yet in databases.

2. **STRING includes functional associations beyond physical PPIs.** STRING edges encompass co-expression, text mining, and genomic context, not only direct physical interactions.

3. **Shortest path may not be biologically most relevant.** The path-minimizing algorithm may choose hub-mediated routes. Hub penalty mitigates but does not eliminate this.

4. **Pan-cancer analysis may confound lineage-specific effects.** A pair may appear significant because the mutation occurs predominantly in one cancer type where the target is also essential. Always perform lineage-specific sensitivity analysis.

5. **Known SL databases are incomplete and biased toward well-studied genes.** Novel SL pairs may lack database support even when real.

6. **Co-dependency does not equal synthetic lethality.** Strongly correlated dependency profiles may reflect shared essentiality, complex co-membership, or technical artifacts (Bug #18).

7. **Directed signaling paths do not prove causality.** A directed edge in OmniPath indicates literature evidence of regulation, not necessarily direct biochemical causation (Bug #13).

8. **Druggability scores reflect current knowledge.** A target classified as "unknown" today may become druggable in the future. Use the therapeutic scoring scheme to up-weight druggability when prioritizing for drug development.

---

## Novelty & Positioning

This framework differs from existing approaches in several key ways:

### 1. Mutation-Conditioned Pair Discovery

The output is not just "EGFR is a dependency gene" but **"ARFGEF3-mutant tumors may be more vulnerable to EGFR inhibition."** This is directly aligned with precision oncology biomarker strategies.

### 2. Multi-Layer Empirical Evidence Integration

Rather than using a single network or database, the framework integrates 7 independent evidence layers, each measuring a distinct aspect of biological plausibility. This creates a transparent, auditable ranking.

### 3. Interpretable Biological Paths

Each candidate is accompanied by an explicit mechanism path (mutation_gene → intermediate → target_gene) that can directly support mechanism figures and experimental hypotheses.

### 4. Confounder-Aware Statistical Framework

The regression model adjusts for lineage, TMB, MSI, CNV burden, and target expression — the major known confounders in CRISPR dependency analysis.

### 5. Degree-Matched Null Model

The permutation test controls for the fact that hub genes are close to everything, providing empirical significance estimates rather than relying on arbitrary network distance cutoffs.

### 6. Therapeutic Prioritization

Druggability scoring enables stratification by translational potential, with separate mechanism-focused and therapeutic-focused scoring schemes.

---

## References

- **STRING:** Szklarczyk D, et al. The STRING database in 2023: protein-protein association networks and functional enrichment analyses. *Nucleic Acids Research*, 2023.
- **Reactome:** Gillespie M, et al. The reactome pathway knowledgebase 2022. *Nucleic Acids Research*, 2022.
- **OmniPath:** Turei D, et al. Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. *Molecular Systems Biology*, 2021.
- **CORUM:** Giurgiu M, et al. CORUM: the comprehensive resource of mammalian protein complexes — 2019. *Nucleic Acids Research*, 2019.
- **SynLethDB:** Guo J, et al. SynLethDB 2.0: a web-based knowledge graph database on synthetic lethality. *Nucleic Acids Research*, 2023.
- **DepMap:** Tsherniak A, et al. Defining a cancer dependency map. *Cell*, 2017.
- **DGIdb:** Freshour SL, et al. Integration of the Drug-Gene Interaction Database. *Nucleic Acids Research*, 2021.
- **Open Targets:** Ochoa D, et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. *Nucleic Acids Research*, 2021.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

*For the complete specification, methodology details, and design rationale, see [project.md](project.md).*
