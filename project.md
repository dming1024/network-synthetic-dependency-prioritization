# Project: Empirical Network-Constrained Synthetic Dependency Prioritization

## 1. Project Overview

### 1.1 Project name

**Empirical Network-Constrained Synthetic Dependency Prioritization**

Alternative names:

- Network-weighted Mutation-Conditioned Dependency Ranking
- Network-Constrained Synthetic Lethality Prioritization
- Empirical Evidence-Weighted Synthetic Dependency Discovery

### 1.2 Core objective

This project aims to prioritize potential **synthetic lethal / synthetic dependency gene pairs** from large-scale mutation and CRISPR dependency data by integrating:

1. Statistical mutation–dependency association;
2. Empirical biological network proximity;
3. Pathway concordance;
4. Directional signaling evidence;
5. Protein complex/module evidence;
6. Known synthetic lethality prior evidence;
7. Target druggability;
8. Confounder correction and network null models.

The central biological question is:

> If mutation of gene A is associated with increased dependency on gene B, can empirical biological networks help determine whether this association is mechanistically plausible and therapeutically actionable?

---

## 2. Biological Rationale

### 2.1 Starting point

The input statistical association usually has the following structure:

```text
mutation_gene    target_gene    beta    p_value    FDR    delta_dependency    n_mut    n_wt
ARFGEF3          EGFR           -0.38   1e-6       ...    -0.34              10       1102
FIRRM            EGFR           -0.43   2e-6       ...    -0.43               8       1104
```

Typical interpretation:

- `mutation_gene`: the mutated gene, such as `ARFGEF3`;
- `target_gene`: the CRISPR dependency gene, such as `EGFR`;
- `beta < 0`: mutation-positive samples show stronger dependency on the target gene, assuming more negative CRISPR gene effect means stronger dependency;
- `FDR`: multiple-testing adjusted significance;
- `delta_dependency`: dependency difference between mutant and wild-type groups;
- `n_mut`: number of mutant samples.

The synthetic dependency hypothesis is:

```text
mutation_gene alteration
        ↓
target_gene dependency becomes stronger
        ↓
target_gene may represent a conditional vulnerability
```

### 2.2 Why network constraint is needed

Pure statistical associations can be misleading because of:

1. Cancer lineage effects;
2. Tumor mutation burden / genomic instability;
3. Microsatellite instability;
4. Copy number burden;
5. Target gene expression differences;
6. Low mutation frequency;
7. Hub gene artifacts;
8. Batch or cohort effects;
9. Passenger mutation effects.

Therefore, empirical biological networks are used as **biological plausibility constraints**, not as final proof.

The logic is:

```text
Statistical association
        +
Network proximity
        +
Pathway concordance
        +
Known synthetic lethality evidence
        +
Targetability
        ↓
Higher-priority synthetic dependency candidate
```

---

## 3. Main Deliverable

The project should produce a ranked candidate table:

```text
mutation_gene
target_gene
beta
FDR
delta_dependency
n_mut
stat_score
STRING_distance
STRING_path
STRING_path_confidence
hub_penalty
network_score
shared_Reactome_pathway
pathway_score
OmniPath_directional_path
signaling_score
complex_score
known_SL_evidence
SL_prior_score
depmap_codependency_score
drug_score
final_SL_score
evidence_class
mechanism_summary
```

Example final output:

```text
ARFGEF3 mutation → EGFR dependency

Statistical evidence:
beta = -0.38, FDR = 0.001, delta_dependency = -0.34

Network evidence:
ARFGEF3 → intermediate_gene → EGFR
distance = 2
path_confidence = 0.64

Pathway evidence:
Shared pathway = membrane trafficking / receptor signaling

Targetability:
EGFR is druggable

Final classification:
Class I or Class II candidate
```

---

## 4. Recommended Project Structure

```text
network_synthetic_dependency_project/
│
├── data/
│   ├── raw/
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
│   ├── processed/
│   │   ├── association_results.csv
│   │   ├── candidate_pairs.csv
│   │   ├── network_features.csv
│   │   ├── pathway_features.csv
│   │   ├── sl_prior_features.csv
│   │   └── final_ranked_candidates.csv
│
├── R/
│   ├── 00_utils.R
│   ├── 01_prepare_data.R
│   ├── 02_association_analysis.R
│   ├── 03_build_networks.R
│   ├── 04_compute_network_features.R
│   ├── 05_compute_empirical_evidence.R
│   ├── 06_score_candidates.R
│   ├── 07_null_model_permutation.R
│   ├── 08_visualization.R
│   └── 09_generate_report.R
│
├── results/
│   ├── tables/
│   ├── figures/
│   └── reports/
│
├── logs/
│   └── run_log.txt
│
├── config/
│   └── config.yaml
│
├── README.md
└── project.md
```

---

## 5. Input Data

### 5.1 Mutation matrix

Rows are samples, columns are genes.

```text
sample_id    TP53    KRAS    ATM    BRCA1    ARFGEF3
S1           1       0       0      1        0
S2           0       1       1      0        1
```

Recommended mutation definition:

- nonsynonymous mutation;
- frameshift insertion/deletion;
- nonsense mutation;
- splice-site mutation;
- damaging missense mutation;
- optionally loss-of-function-only matrix.

Avoid mixing silent mutations unless specifically needed.

### 5.2 Dependency matrix

Rows are samples/cell lines, columns are genes.

```text
sample_id    EGFR    ATM    PARP1    CDK6
S1           -0.45   -0.12  -0.78    -0.20
S2           -0.10   -0.35  -0.15    -0.60
```

Interpretation:

- More negative gene effect indicates stronger dependency.
- `beta < 0` in regression generally suggests mutation-positive samples are more dependent on the target gene.

### 5.3 Sample metadata

Recommended columns:

```text
sample_id
lineage
cancer_type
TMB
MSI_status
CNV_burden
purity
batch
source_dataset
```

For cell-line dependency data, at least include:

```text
sample_id
lineage
subtype
TMB
MSI_status
CNV_burden
```

### 5.4 Expression matrix

Optional but highly recommended.

Purpose:

- Correct for target gene expression;
- Avoid false dependency associations caused by target gene not being expressed;
- Improve biological interpretability.

Regression covariate:

```text
dependency_target ~ mutation_gene + expression_target + lineage + TMB + MSI + CNV_burden
```

---

## 6. Recommended Empirical Networks

### 6.1 STRING functional association network

Purpose:

- Broad functional relationship detection;
- Shortest path calculation;
- Path confidence scoring;
- Common-neighbor calculation;
- Network propagation.

Recommended filters:

```text
combined_score >= 700  # high confidence
combined_score >= 900  # very high confidence, stricter but sparser
```

Important note:

STRING is not a pure physical PPI network. It contains physical interactions and broader functional associations.

### 6.2 Reactome pathway network

Purpose:

- Curated pathway concordance;
- Mechanism-level interpretation;
- Shared biological process identification.

Scoring:

```text
same exact pathway = 1
same parent pathway = 0.5
no shared pathway = 0
```

### 6.3 OmniPath / SIGNOR directional signaling network

Purpose:

- Directional signal flow;
- Activation/inhibition relationship;
- Mechanistic chain inference.

Example:

```text
mutation_gene → signaling intermediate → target_gene
```

Scoring:

```text
directed path length <= 3 = 1
directed path length 4–5 = 0.5
no directed path = 0
```

### 6.4 CORUM / protein complex network

Purpose:

- Detect same-complex relationships;
- Detect complex-level dependency;
- Especially useful for chromatin remodeling, DNA repair, proteasome, spliceosome, ribosome, cohesin, and mediator complex.

Scoring:

```text
same protein complex = 1
same complex family/module = 0.5
no complex relationship = 0
```

### 6.5 Synthetic lethality prior databases

Potential sources:

- SynLethDB;
- SLKG;
- literature-curated SL pairs;
- internal experimentally validated SL pairs;
- DepMap-derived co-dependency or conditional dependency evidence.

Scoring:

```text
known SL pair = 1
neighbor-supported SL = 0.5
pathway-level SL support = 0.3
no evidence = 0
```

### 6.6 Druggability resources

Potential sources:

- Open Targets;
- DGIdb;
- DrugBank;
- Tclin/Tchem target lists;
- kinase/receptor/enzyme annotations;
- FDA-approved drug target lists.

Scoring:

```text
FDA-approved drug target = 1
clinical/preclinical inhibitor = 0.7
drug target class, such as kinase/enzyme/receptor = 0.5
unknown = 0
```

---

## 7. Module 1: Statistical Association Analysis

### 7.1 Goal

Identify mutation–dependency associations:

```text
mutation_gene A is associated with dependency change of target_gene B
```

### 7.2 Recommended regression model

For each pair `A mutation → B dependency`:

```r
dependency_B ~ mutation_A + expression_B + lineage + TMB + MSI + CNV_burden
```

If expression is unavailable:

```r
dependency_B ~ mutation_A + lineage + TMB + MSI + CNV_burden
```

### 7.3 Pseudocode

```r
for each mutation_gene A:
  for each dependency_gene B:

    df = merge mutation_A, dependency_B, expression_B, metadata

    if number of A-mut samples < min_mut_count:
        skip

    fit = lm(
      dependency_B ~ mutation_A + expression_B + lineage + TMB + MSI + CNV_burden,
      data = df
    )

    beta = coefficient of mutation_A
    p_value = p value of mutation_A
    delta_dependency = mean(dependency_B in A-mut) - mean(dependency_B in A-wt)

    save A, B, beta, p_value, delta_dependency, n_mut, n_wt

adjust all p values by BH-FDR
```

### 7.4 Initial candidate filtering

```r
candidate_df <- assoc_df %>%
  filter(
    FDR < 0.05,
    beta < 0,
    n_mut >= 5,
    abs(delta_dependency) >= 0.2
  )
```

Recommended thresholds can be adjusted:

```text
n_mut >= 5   # exploratory
n_mut >= 10  # more reliable
abs(delta_dependency) >= 0.2 or 0.3
FDR < 0.05 or FDR < 0.1 for discovery
```

---

## 8. Module 2: Build Biological Networks

### 8.1 Build STRING graph

```r
build_string_graph <- function(string_edges, score_cutoff = 700) {
  string_use <- string_edges %>%
    filter(combined_score >= score_cutoff) %>%
    mutate(
      edge_confidence = combined_score / 1000,
      distance_weight = 1 / edge_confidence
    ) %>%
    distinct(geneA, geneB, .keep_all = TRUE)

  graph <- igraph::graph_from_data_frame(
    string_use %>% select(geneA, geneB, edge_confidence, distance_weight),
    directed = FALSE
  )

  return(graph)
}
```

### 8.2 Build directional signaling graph

```r
build_signaling_graph <- function(signaling_edges) {
  signaling_use <- signaling_edges %>%
    mutate(
      edge_confidence = ifelse(is.na(confidence), 1, confidence),
      distance_weight = 1 / edge_confidence
    )

  graph <- igraph::graph_from_data_frame(
    signaling_use %>% select(source_gene, target_gene, effect, edge_confidence, distance_weight),
    directed = TRUE
  )

  return(graph)
}
```

### 8.3 Important preprocessing checks

Before building networks:

```r
# Harmonize gene symbols
# Remove duplicated edges
# Remove self-loops
# Ensure all gene names are uppercase if using human symbols
# Check if mutation genes and dependency genes are covered by the network
```

---

## 9. Module 3: Network Feature Calculation

### 9.1 Features to calculate

For each pair:

```text
connected
shortest_distance
path_genes
path_confidence
common_neighbors
direct_interaction
hub_penalty
network_score
empirical_network_p
```

### 9.2 Pseudocode

```r
compute_network_features <- function(graph, geneA, geneB, max_path_len = 5) {

  if geneA not in graph or geneB not in graph:
      return disconnected features

  path = shortest path from geneA to geneB using distance_weight

  if no path:
      return disconnected features

  path_genes = genes along path
  path_edges = edges along path
  distance = number of edges
  edge_scores = confidence scores of path edges

  if distance > max_path_len:
      path_confidence = 0
  else:
      path_confidence = product(edge_scores)

  intermediate_genes = path_genes excluding geneA and geneB

  if no intermediate genes:
      hub_penalty = 1
  else:
      mean_degree = mean degree of intermediate genes
      hub_penalty = 1 / log2(mean_degree + 2)

  distance_score = convert distance to score
  network_score = distance_score * path_confidence * hub_penalty

  common_neighbors = number of shared direct neighbors
  direct_interaction = distance == 1

  return all features
}
```

### 9.3 Distance score

```r
convert_distance_to_score <- function(distance) {
  dplyr::case_when(
    distance == 1 ~ 1.0,
    distance == 2 ~ 0.7,
    distance == 3 ~ 0.4,
    distance <= 5 ~ 0.1,
    TRUE ~ 0
  )
}
```

### 9.4 Hub penalty

Formula:

```text
hub_penalty = 1 / log2(mean_degree_of_intermediate_nodes + 2)
```

Reason:

Paths through hub genes such as `TP53`, `AKT1`, `MYC`, `EGFR`, `SRC`, or `UBC` may be biologically less specific.

---

## 10. Module 4: Pathway Concordance

### 10.1 Goal

Determine whether mutation gene and dependency gene are involved in the same pathway or parent biological process.

### 10.2 Pseudocode

```r
compute_pathway_score <- function(geneA, geneB, reactome_df) {

  pathways_A = reactome_df %>% filter(gene == geneA) %>% pull(pathway)
  pathways_B = reactome_df %>% filter(gene == geneB) %>% pull(pathway)

  shared = intersect(pathways_A, pathways_B)

  if length(shared) > 0:
      pathway_score = 1
      shared_pathways = paste(shared, collapse = '; ')
  else if share_parent_pathway(geneA, geneB):
      pathway_score = 0.5
      shared_pathways = parent_pathway_name
  else:
      pathway_score = 0
      shared_pathways = NA

  return pathway_score and shared_pathways
}
```

### 10.3 Potential bug

Reactome pathway names may be too specific. Two genes may share a parent pathway but not the exact child pathway. Therefore, if possible, keep both:

```text
pathway_id
pathway_name
parent_pathway_id
parent_pathway_name
```

---

## 11. Module 5: Directional Signaling Evidence

### 11.1 Goal

Use OmniPath/SIGNOR-like directed networks to determine whether a directed relationship exists between `mutation_gene` and `target_gene`.

### 11.2 Pseudocode

```r
compute_signaling_score <- function(signaling_graph, geneA, geneB, max_len = 5) {

  path_A_to_B = directed shortest path from geneA to geneB
  path_B_to_A = directed shortest path from geneB to geneA

  if path_A_to_B exists and length <= 3:
      signaling_score = 1
      signaling_direction = 'mutation_gene_to_target_gene'
      signaling_path = path_A_to_B

  else if path_B_to_A exists and length <= 3:
      signaling_score = 0.7
      signaling_direction = 'target_gene_to_mutation_gene'
      signaling_path = path_B_to_A

  else if any directed path exists and length <= 5:
      signaling_score = 0.5
      signaling_direction = 'weak_directional_support'

  else:
      signaling_score = 0
      signaling_direction = 'none'
      signaling_path = NA

  return signaling_score, direction, path
}
```

### 11.3 Interpretation warning

A directed signaling path does not prove causality. It only indicates that both genes are embedded in a known regulatory chain.

---

## 12. Module 6: Protein Complex Evidence

### 12.1 Goal

Detect whether mutation and target genes belong to the same complex or related complex module.

### 12.2 Pseudocode

```r
compute_complex_score <- function(geneA, geneB, corum_df) {

  complex_A = corum_df %>% filter(gene == geneA) %>% pull(complex_id)
  complex_B = corum_df %>% filter(gene == geneB) %>% pull(complex_id)

  shared_complex = intersect(complex_A, complex_B)

  if length(shared_complex) > 0:
      complex_score = 1
      complex_evidence = paste(shared_complex, collapse = '; ')
  else:
      complex_score = 0
      complex_evidence = NA

  return complex_score, complex_evidence
}
```

### 12.3 Important biological note

Same-complex relationships are especially informative for:

- SWI/SNF complex;
- DNA repair complexes;
- spliceosome;
- proteasome;
- ribosome biogenesis;
- cohesin;
- mediator complex;
- chromatin modifier complexes.

---

## 13. Module 7: Synthetic Lethality Prior Evidence

### 13.1 Goal

Check whether the candidate pair is already known or indirectly supported by synthetic lethality databases.

### 13.2 Scoring strategy

```text
known SL pair = 1
neighbor-supported SL = 0.5
pathway-level SL support = 0.3
no prior evidence = 0
```

### 13.3 Pseudocode

```r
compute_sl_prior_score <- function(geneA, geneB, sl_db, string_graph, reactome_df) {

  if unordered_pair(geneA, geneB) in sl_db:
      sl_prior_score = 1
      sl_prior_type = 'known_SL_pair'

  else:
      neighbors_A = neighbors of geneA in STRING
      neighbors_B = neighbors of geneB in STRING

      if any unordered_pair(neighbor_A, geneB) in sl_db or
         any unordered_pair(geneA, neighbor_B) in sl_db:
          sl_prior_score = 0.5
          sl_prior_type = 'neighbor_supported_SL'

      else if pathway-level SL support exists:
          sl_prior_score = 0.3
          sl_prior_type = 'pathway_level_SL_support'

      else:
          sl_prior_score = 0
          sl_prior_type = 'none'

  return sl_prior_score, sl_prior_type
}
```

---

## 14. Module 8: DepMap Co-dependency Evidence

### 14.1 Goal

Determine whether mutation gene and target gene have related dependency profiles or whether target gene belongs to a dependency module relevant to mutation gene.

### 14.2 Possible features

```text
correlation between dependency profiles of geneA and geneB
co-essentiality module membership
co-dependency network distance
lineage-specific co-dependency
```

### 14.3 Pseudocode

```r
compute_codependency_score <- function(geneA, geneB, dependency_matrix) {

  if geneA and geneB are both in dependency_matrix columns:
      r = correlation(dependency_matrix[, geneA], dependency_matrix[, geneB], use = 'pairwise.complete.obs')
      codependency_score = rescale(abs(r), 0, 1)
  else:
      r = NA
      codependency_score = 0

  return r, codependency_score
}
```

### 14.4 Interpretation warning

Co-dependency does not necessarily indicate synthetic lethality. It may indicate same complex membership or shared essentiality.

---

## 15. Module 9: Druggability Evidence

### 15.1 Goal

Prioritize target genes with potential therapeutic value.

### 15.2 Pseudocode

```r
compute_druggability_score <- function(target_gene, drug_df) {

  if target_gene in FDA_approved_targets:
      drug_score = 1
      drug_category = 'FDA_approved_target'

  else if target_gene in clinical_trial_targets:
      drug_score = 0.8
      drug_category = 'clinical_target'

  else if target_gene in preclinical_targets:
      drug_score = 0.7
      drug_category = 'preclinical_target'

  else if target_gene is kinase or enzyme or receptor:
      drug_score = 0.5
      drug_category = 'druggable_class'

  else:
      drug_score = 0
      drug_category = 'unknown'

  return drug_score, drug_category
}
```

---

## 16. Module 10: Final Scoring

### 16.1 Statistical score

```r
stat_score_raw = -log10(FDR + 1e-300) * abs(beta)
stat_score = scale_0_1(stat_score_raw)
```

### 16.2 Network score

```r
network_score = distance_score * path_confidence * hub_penalty
network_score_scaled = scale_0_1(network_score)
```

### 16.3 Final synthetic dependency score

Balanced version:

```text
Final_SL_score =
  0.35 × statistical_score
+ 0.25 × network_score_adjusted
+ 0.15 × pathway_score
+ 0.10 × signaling_score
+ 0.05 × complex_score
+ 0.05 × SL_prior_score
+ 0.05 × druggability_score
```

Mechanism-focused version:

```text
Final_mechanism_score =
  0.30 × statistical_score
+ 0.30 × network_score_adjusted
+ 0.20 × pathway_score
+ 0.10 × signaling_score
+ 0.10 × SL_prior_score
```

Therapeutic-focused version:

```text
Final_therapeutic_score =
  0.30 × statistical_score
+ 0.20 × network_score_adjusted
+ 0.15 × pathway_score
+ 0.10 × signaling_score
+ 0.10 × SL_prior_score
+ 0.15 × druggability_score
```

### 16.4 R-style pseudocode

```r
feature_df <- feature_df %>%
  mutate(
    stat_score = scale_0_1(-log10(FDR + 1e-300) * abs(beta)),
    network_score_scaled = scale_0_1(network_score),

    final_SL_score =
      0.35 * stat_score +
      0.25 * network_score_scaled +
      0.15 * pathway_score +
      0.10 * signaling_score +
      0.05 * complex_score +
      0.05 * sl_prior_score +
      0.05 * druggability_score
  ) %>%
  arrange(desc(final_SL_score))
```

---

## 17. Evidence Classification

### 17.1 Recommended classes

```r
feature_df <- feature_df %>%
  mutate(
    evidence_class = case_when(

      FDR < 0.05 &
      beta < 0 &
      network_distance <= 2 &
      pathway_score >= 0.5 &
      final_SL_score >= 0.75 ~
      'Class I: high-confidence synthetic dependency candidate',

      FDR < 0.05 &
      beta < 0 &
      network_distance <= 3 &
      final_SL_score >= 0.60 ~
      'Class II: network-supported candidate',

      FDR < 0.05 &
      beta < 0 &
      sl_prior_score > 0 ~
      'Class III: SL-prior-supported candidate',

      FDR < 0.05 &
      beta < 0 &
      network_score_scaled == 0 ~
      'Class IV: statistical-only candidate',

      TRUE ~ 'Low priority'
    )
  )
```

### 17.2 Interpretation

```text
Class I:
Strong statistical association plus strong biological network/pathway support.
Best candidates for mechanism writing and experimental validation.

Class II:
Good statistical and network support, but less complete evidence.
Useful for hypothesis generation.

Class III:
Supported by known SL prior, but may not be close in current network.
Useful for validation or rediscovery.

Class IV:
Statistically significant but no network evidence.
Treat carefully; may represent novel biology or false positives.
```

---

## 18. Null Model and Permutation Test

### 18.1 Why this is important

Network scores can be biased because:

- Some genes are hubs;
- Some targets are heavily studied;
- Some networks are dense in cancer-related genes;
- High-degree genes are close to almost everything.

Therefore, calculate empirical network significance.

### 18.2 Null model idea

For each real pair `A–B`:

```text
Keep A fixed.
Randomly sample B_random from genes matched by degree or dependency frequency.
Calculate network_score(A, B_random).
Repeat 1000 times.
Compare observed network_score(A, B) with random scores.
```

### 18.3 Pseudocode

```r
compute_empirical_network_p <- function(graph, geneA, geneB, candidate_targets, n_perm = 1000) {

  observed_score = compute_network_score(graph, geneA, geneB)

  random_scores = c()

  for i in 1:n_perm:
      B_random = sample(candidate_targets, size = 1)
      random_score = compute_network_score(graph, geneA, B_random)
      random_scores = c(random_scores, random_score)

  empirical_p = (sum(random_scores >= observed_score) + 1) / (n_perm + 1)

  return observed_score, empirical_p
}
```

### 18.4 Better matched null model

Instead of random target sampling, match by:

```text
target gene degree
target dependency variance
target mean dependency
target expression level
target essentiality class
```

Recommended later version:

```r
B_random <- sample_targets_matched_by_degree_and_dependency(geneB)
```

---

## 19. Visualization Plan

### 19.1 Ranking dot plot

X-axis:

```text
final_SL_score
```

Y-axis:

```text
mutation_gene → target_gene
```

Color:

```text
evidence_class
```

Size:

```text
abs(beta) or -log10(FDR)
```

### 19.2 Network path plot

For top candidates:

```text
mutation_gene → intermediate_gene_1 → intermediate_gene_2 → target_gene
```

Node color:

```text
mutation gene / target gene / intermediate gene
```

Edge width:

```text
STRING confidence score
```

### 19.3 Evidence heatmap

Rows:

```text
candidate pairs
```

Columns:

```text
stat_score
network_score
pathway_score
signaling_score
complex_score
SL_prior_score
drug_score
```

### 19.4 Candidate mechanism card

For each top candidate:

```text
Candidate: ARFGEF3 mutation → EGFR dependency
Statistical evidence: beta, FDR, delta
Network path: ARFGEF3 → X → EGFR
Pathway evidence: shared pathway
SL prior: yes/no
Druggability: yes/no
Interpretation: one-paragraph mechanism hypothesis
Suggested validation: CRISPR, inhibitor, rescue assay
```

---

## 20. Recommended MVP Development Plan

### Version 0.1: Minimal viable prototype

Use only:

```text
association_results.csv
STRING network
Reactome pathway
SynLethDB pairs
drug target list
```

Output:

```text
final_ranked_candidates.csv
candidate mechanism summaries
basic plots
```

### Version 0.2: Improved biological network support

Add:

```text
OmniPath directional signaling
CORUM protein complex
DepMap co-dependency
```

### Version 0.3: Robustness and null model

Add:

```text
permutation empirical network p-value
degree-matched randomization
lineage-specific sensitivity analysis
```

### Version 0.4: Report generation

Add:

```text
automatic markdown/html report
candidate mechanism cards
network plots
publication-style figure legends
```

### Version 1.0: Reusable package/pipeline

Add:

```text
config file
command-line interface
standardized input checks
unit tests
example dataset
full documentation
```

---

## 21. Possible Bugs and Troubleshooting

### Bug 1: Gene symbols do not match

Symptoms:

```text
Many genes are not found in the network.
Most pairs are disconnected.
```

Causes:

- ENSG IDs mixed with gene symbols;
- ENSP IDs from STRING not mapped to symbols;
- old gene symbols;
- lowercase/uppercase inconsistency;
- aliases not resolved.

Solutions:

```r
# Standardize symbols before analysis
# Use org.Hs.eg.db, biomaRt, HGNChelper, or STRING mapping file
# Convert all symbols to approved HGNC gene symbols
```

Checklist:

```text
Percentage of mutation genes in network > 70%
Percentage of target genes in network > 70%
```

---

### Bug 2: STRING IDs are protein IDs, not gene symbols

Symptoms:

```text
9606.ENSP00000354587 appears instead of TP53.
```

Solution:

Use STRING protein alias mapping:

```text
protein_id → preferred_name / gene symbol
```

Then collapse duplicated gene-pair edges by keeping the maximum combined score.

---

### Bug 3: Too many pairs are connected through hub genes

Symptoms:

```text
Most significant pairs show paths through TP53, AKT1, MYC, UBC, SRC.
```

Solutions:

1. Apply hub penalty;
2. Remove extreme hub genes in sensitivity analysis;
3. Use degree-matched permutation;
4. Report hub-mediated paths separately;
5. Avoid overinterpreting hub-only links.

---

### Bug 4: Shortest path returns biologically meaningless route

Cause:

Shortest path optimizes distance, not biological specificity.

Solutions:

- Use high-confidence STRING edges only;
- Use pathway-restricted paths;
- Penalize hub nodes;
- Compare multiple shortest paths;
- Use random walk with restart as an alternative.

---

### Bug 5: Path confidence becomes extremely small for long paths

Cause:

Multiplying many edge scores quickly approaches zero.

Solutions:

```r
# Option 1: Use geometric mean
path_confidence = prod(edge_scores)^(1 / length(edge_scores))

# Option 2: Use log score
path_confidence_log = sum(log(edge_scores + 1e-10))

# Option 3: Truncate paths longer than 5
if distance > 5: path_confidence = 0
```

Recommended first version:

```text
Use product for paths <= 5 and set longer paths to zero.
```

---

### Bug 6: FDR is zero due to numerical underflow

Symptoms:

```r
-log10(FDR) returns Inf
```

Solution:

```r
-log10(FDR + 1e-300)
```

---

### Bug 7: Scaling fails when all values are identical

Symptoms:

```text
scale_0_1 returns NaN.
```

Solution:

```r
scale_0_1 <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}
```

---

### Bug 8: Regression fails because mutation group is too small

Symptoms:

```text
lm coefficients are NA.
Model matrix is singular.
```

Solutions:

- Require `n_mut >= 5` or `n_mut >= 10`;
- Remove genes with extremely low mutation frequency;
- Use lineage-stratified analysis only when enough samples exist;
- Consider Firth logistic/robust regression for special cases.

---

### Bug 9: Lineage confounding dominates associations

Symptoms:

```text
Top pairs are actually lineage markers.
```

Solutions:

1. Include lineage as covariate;
2. Perform within-lineage analysis;
3. Use meta-analysis across lineages;
4. Require effect consistency across multiple lineages;
5. Compare pan-cancer versus lineage-specific results.

---

### Bug 10: Target dependency is driven by target expression

Symptoms:

```text
Mutation-positive samples have higher target expression, and dependency also appears stronger.
```

Solution:

Include target expression as covariate:

```r
dependency_B ~ mutation_A + expression_B + lineage + TMB + MSI + CNV_burden
```

Also report:

```text
target_expression_beta
target_expression_p
```

---

### Bug 11: Network evidence over-penalizes novel biology

Problem:

Novel synthetic lethal pairs may not be well represented in existing databases.

Solution:

Do not discard network-unsupported candidates. Classify them as:

```text
Class IV: statistical-only candidate
```

Then manually inspect high-effect candidates.

---

### Bug 12: Known SL database creates circular validation

Problem:

If known SL pairs are used in scoring, they cannot be used again as independent validation.

Solution:

Separate:

```text
SL_prior_score: used for ranking
External validation: independent dataset or held-out SL database
```

---

### Bug 13: Directional signaling path is interpreted too strongly

Problem:

A directed path does not prove mutation A causally regulates dependency B.

Solution:

Use cautious language:

```text
The candidate pair is embedded in a known directional signaling cascade.
```

Avoid saying:

```text
A directly regulates B.
```

---

### Bug 14: Druggability dominates final ranking

Problem:

Famous drug targets always rank high.

Solution:

Use small druggability weight in the main score, such as 0.05–0.10.

Optionally create separate scores:

```text
mechanism_score
therapeutic_score
```

---

### Bug 15: Many missing values after joining feature tables

Symptoms:

```text
NA pathway_score, NA network_score, NA drug_score.
```

Solution:

After joining:

```r
feature_df <- feature_df %>%
  mutate(
    pathway_score = replace_na(pathway_score, 0),
    network_score = replace_na(network_score, 0),
    signaling_score = replace_na(signaling_score, 0),
    complex_score = replace_na(complex_score, 0),
    sl_prior_score = replace_na(sl_prior_score, 0),
    druggability_score = replace_na(druggability_score, 0)
  )
```

---

### Bug 16: Multiple networks use different edge meanings

Problem:

STRING, Reactome, OmniPath, and CORUM edges do not mean the same thing.

Solution:

Do not merge all networks naively as one graph in the first version. Keep separate layer-specific scores:

```text
STRING_score
Reactome_score
OmniPath_score
CORUM_score
SL_prior_score
```

Then integrate at the scoring level.

---

### Bug 17: Pan-cancer analysis produces misleading synthetic lethality

Problem:

A pair may be significant only because one mutation occurs in one cancer type and the dependency target is essential in that same lineage.

Solution:

Use:

```r
dependency ~ mutation + lineage + covariates
```

And perform:

```text
within-lineage sensitivity analysis
leave-one-lineage-out analysis
```

---

### Bug 18: Candidate target is universally essential

Problem:

Some targets are essential in almost all cell lines, not conditionally dependent.

Solution:

Add selectivity filter:

```text
target dependency variance
difference between mutation-positive and mutation-negative groups
selective dependency score
```

Avoid prioritizing pan-essential genes unless the mutation-specific effect is strong.

---

### Bug 19: Copy number confounding in CRISPR dependency

Problem:

CRISPR gene effect can be biased by copy number amplification.

Solution:

Include:

```text
target_gene_CNV
CNV_burden
```

Or remove target genes located in highly amplified regions for sensitivity analysis.

---

### Bug 20: Too many pairwise tests

Problem:

All mutation genes × all dependency genes can be huge.

Solutions:

1. Pre-filter mutation genes by frequency;
2. Pre-filter dependency genes by variance/selectivity;
3. Run parallel computation;
4. Save intermediate results;
5. Use chunk-based processing.

---

## 22. Suggested Utility Functions

### 22.1 Scaling function

```r
scale_0_1 <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}
```

### 22.2 Pair key function

```r
make_pair_key <- function(geneA, geneB) {
  apply(cbind(geneA, geneB), 1, function(x) {
    paste(sort(x), collapse = '__')
  })
}
```

### 22.3 Safe shortest path function

```r
safe_shortest_path <- function(graph, from, to, weights = NULL) {
  if (!(from %in% V(graph)$name) || !(to %in% V(graph)$name)) {
    return(NULL)
  }

  sp <- tryCatch({
    igraph::shortest_paths(graph, from = from, to = to, weights = weights, output = 'both')
  }, error = function(e) NULL)

  return(sp)
}
```

---

## 23. Example Full Pipeline Pseudocode

```r
# 0. Load packages
library(tidyverse)
library(igraph)

# 1. Load data
mutation_mat <- read.csv('data/raw/mutation_matrix.csv')
dependency_mat <- read.csv('data/raw/dependency_matrix.csv')
expression_mat <- read.csv('data/raw/expression_matrix.csv')
metadata <- read.csv('data/raw/sample_metadata.csv')
string_edges <- read.csv('data/raw/STRING_edges.csv')
reactome_df <- read.csv('data/raw/Reactome_gene_pathway.csv')
signaling_edges <- read.csv('data/raw/OmniPath_edges.csv')
corum_df <- read.csv('data/raw/CORUM_complexes.csv')
sl_db <- read.csv('data/raw/SynLethDB_pairs.csv')
drug_df <- read.csv('data/raw/druggable_targets.csv')

# 2. Run association analysis
assoc_df <- run_mutation_dependency_association(
  mutation_mat = mutation_mat,
  dependency_mat = dependency_mat,
  expression_mat = expression_mat,
  metadata = metadata,
  min_mut_count = 5
)

# 3. Filter candidates
candidate_df <- assoc_df %>%
  filter(FDR < 0.05, beta < 0, n_mut >= 5, abs(delta_dependency) >= 0.2)

# 4. Build networks
string_graph <- build_string_graph(string_edges, score_cutoff = 700)
signaling_graph <- build_signaling_graph(signaling_edges)

# 5. Compute STRING network features
network_features <- candidate_df %>%
  select(mutation_gene, target_gene) %>%
  distinct() %>%
  pmap_dfr(function(mutation_gene, target_gene) {
    compute_network_features(string_graph, mutation_gene, target_gene)
  })

# 6. Compute pathway features
pathway_features <- candidate_df %>%
  select(mutation_gene, target_gene) %>%
  distinct() %>%
  pmap_dfr(function(mutation_gene, target_gene) {
    compute_pathway_score(mutation_gene, target_gene, reactome_df)
  })

# 7. Compute signaling features
signaling_features <- candidate_df %>%
  select(mutation_gene, target_gene) %>%
  distinct() %>%
  pmap_dfr(function(mutation_gene, target_gene) {
    compute_signaling_score(signaling_graph, mutation_gene, target_gene)
  })

# 8. Compute complex features
complex_features <- candidate_df %>%
  select(mutation_gene, target_gene) %>%
  distinct() %>%
  pmap_dfr(function(mutation_gene, target_gene) {
    compute_complex_score(mutation_gene, target_gene, corum_df)
  })

# 9. Compute SL prior features
sl_features <- candidate_df %>%
  select(mutation_gene, target_gene) %>%
  distinct() %>%
  pmap_dfr(function(mutation_gene, target_gene) {
    compute_sl_prior_score(mutation_gene, target_gene, sl_db, string_graph, reactome_df)
  })

# 10. Compute druggability features
drug_features <- candidate_df %>%
  select(target_gene) %>%
  distinct() %>%
  rowwise() %>%
  do(compute_druggability_score(.$target_gene, drug_df))

# 11. Merge features
feature_df <- candidate_df %>%
  left_join(network_features, by = c('mutation_gene', 'target_gene')) %>%
  left_join(pathway_features, by = c('mutation_gene', 'target_gene')) %>%
  left_join(signaling_features, by = c('mutation_gene', 'target_gene')) %>%
  left_join(complex_features, by = c('mutation_gene', 'target_gene')) %>%
  left_join(sl_features, by = c('mutation_gene', 'target_gene')) %>%
  left_join(drug_features, by = 'target_gene') %>%
  replace_na(list(
    network_score = 0,
    pathway_score = 0,
    signaling_score = 0,
    complex_score = 0,
    sl_prior_score = 0,
    druggability_score = 0
  ))

# 12. Calculate final score
ranked_df <- feature_df %>%
  mutate(
    stat_score = scale_0_1(-log10(FDR + 1e-300) * abs(beta)),
    network_score_scaled = scale_0_1(network_score),
    final_SL_score =
      0.35 * stat_score +
      0.25 * network_score_scaled +
      0.15 * pathway_score +
      0.10 * signaling_score +
      0.05 * complex_score +
      0.05 * sl_prior_score +
      0.05 * druggability_score
  ) %>%
  arrange(desc(final_SL_score))

# 13. Classify candidates
ranked_df <- classify_candidates(ranked_df)

# 14. Save results
write.csv(ranked_df, 'data/processed/final_ranked_candidates.csv', row.names = FALSE)
```

---

## 24. Biological Interpretation Template

For each top candidate, generate a structured interpretation:

```text
Candidate: [mutation_gene] mutation → [target_gene] dependency

1. Statistical association:
[mutation_gene]-mutant samples showed stronger dependency on [target_gene]
(beta = ..., FDR = ..., delta_dependency = ...).

2. Network support:
[mutation_gene] and [target_gene] are connected in the STRING network through:
[path_genes]
with path confidence = ... and distance = ... .

3. Pathway support:
The two genes share or converge on [pathway_name], suggesting potential functional coupling.

4. Synthetic lethality prior:
Known SL evidence is [present/absent/neighbor-supported].

5. Therapeutic relevance:
[target_gene] is classified as [drug_category], suggesting potential druggability.

6. Working hypothesis:
Mutation of [mutation_gene] may rewire [pathway/process], increasing cellular reliance on [target_gene].
This supports [target_gene] as a candidate synthetic dependency in [mutation_gene]-altered contexts.

7. Suggested validation:
- Compare dependency in mutant versus wild-type models;
- Knockout or inhibit target_gene in mutation-positive and mutation-negative cells;
- Rescue mutation_gene or pathway intermediate;
- Test pathway biomarkers;
- Evaluate drug sensitivity if inhibitor exists.
```

---

## 25. Suggested Experimental Validation Strategy

### 25.1 Computational validation

1. Check association in independent dependency datasets;
2. Stratify by lineage;
3. Adjust for target expression and CNV;
4. Validate in DepMap CRISPR and RNAi separately;
5. Compare with drug sensitivity data;
6. Check whether target dependency is selective, not pan-essential;
7. Compare with known SL databases.

### 25.2 Wet-lab validation

Recommended design:

```text
mutation-positive cell line
mutation-negative matched control cell line
        ↓
target_gene knockout or inhibitor treatment
        ↓
cell viability / colony formation / apoptosis
        ↓
rescue experiment or pathway biomarker assay
```

For druggable target:

```text
mutation-positive vs mutation-negative cells
        ↓
target inhibitor dose-response
        ↓
IC50 / AUC comparison
```

For non-druggable target:

```text
CRISPR knockout / siRNA knockdown
        ↓
viability or proliferation assay
```

---

## 26. How to Position the Novelty

The novelty should not be described as simply “using PPI network.” Many studies already use PPI or pathway networks.

The stronger novelty is:

### 26.1 Mutation-conditioned dependency pair discovery

The project focuses on:

```text
mutation_gene alteration → target_gene dependency change
```

rather than simply identifying universally essential targets.

### 26.2 Pair-level vulnerability prioritization

The output is not only:

```text
EGFR is a dependency gene
```

but:

```text
ARFGEF3-mutant tumors may be more vulnerable to EGFR dependency/inhibition
```

This is more directly aligned with precision oncology.

### 26.3 Network-constrained interpretability

Each candidate pair is accompanied by an interpretable biological path:

```text
mutation_gene → intermediate genes/pathway → target_gene
```

This can directly support mechanism figures and experimental hypotheses.

### 26.4 Multi-layer empirical evidence integration

The method integrates multiple evidence layers:

```text
statistical association
STRING proximity
Reactome pathway
OmniPath directionality
CORUM complex
known SL prior
DepMap co-dependency
drug target evidence
```

This creates a more transparent ranking framework than pure statistical screening.

### 26.5 Therapeutic prioritization

By adding druggability and target selectivity, the method prioritizes candidates with translational potential.

---

## 27. Recommended Manuscript/Report Method Description

```text
We developed a network-constrained prioritization framework to rank mutation-conditioned CRISPR dependency candidates. Candidate mutation–dependency pairs were first identified using regression-based association analysis after adjusting for lineage, target expression, tumor mutation burden, microsatellite instability, and copy-number burden. Each candidate pair was then projected onto multiple empirical biological networks, including functional association, curated pathway, directional signaling, protein complex, and synthetic lethality prior networks. For each pair, we calculated network distance, path confidence, hub-adjusted network proximity, pathway concordance, directional signaling support, complex membership, known synthetic lethality evidence, and target druggability. These features were integrated with statistical effect size and FDR to generate a final synthetic dependency score, enabling prioritization of mutation–target pairs with both statistical support and mechanistic plausibility.
```

---

## 28. Key Design Principles

1. **Do not use network evidence as proof.** Use it as biological plausibility weighting.
2. **Do not discard statistical-only candidates.** They may represent novel biology.
3. **Always correct for lineage and genomic confounders.**
4. **Penalize hub-mediated network paths.**
5. **Keep network layers separate in early versions.** Do not merge all edges blindly.
6. **Use empirical permutation tests to evaluate whether network proximity is stronger than random expectation.**
7. **Distinguish mechanism score from therapeutic score.**
8. **Prioritize interpretable outputs.** Each top pair should have a clear path and mechanism hypothesis.
9. **Validate top candidates in independent data and experiments.**
10. **Be cautious with language.** Say “candidate synthetic dependency” rather than “confirmed synthetic lethality.”

---

## 29. Minimal Checklist Before Running

```text
[ ] Mutation matrix is binary and sample IDs are clean.
[ ] Dependency matrix uses the same sample IDs.
[ ] Metadata contains lineage or cancer type.
[ ] Gene symbols are standardized.
[ ] Mutation genes with too few mutant samples are removed.
[ ] Dependency genes with very low variance are removed.
[ ] STRING IDs are converted to gene symbols.
[ ] Network coverage is checked.
[ ] Regression covariates are available.
[ ] FDR correction is applied globally.
[ ] Missing feature scores are replaced with zero.
[ ] Hub penalty is enabled.
[ ] Final score weights are documented.
[ ] Top candidates are manually inspected.
```

---

## 30. Final Suggested First Development Task

Build the MVP in this order:

```text
1. Clean association result table.
2. Load STRING network and compute shortest path features.
3. Add Reactome pathway score.
4. Add SynLethDB known SL score.
5. Add druggability score.
6. Calculate final score.
7. Export ranked table.
8. Generate candidate mechanism cards for top 20 pairs.
```

After MVP works, add:

```text
1. OmniPath directionality.
2. CORUM complex evidence.
3. DepMap co-dependency.
4. Empirical permutation p-value.
5. Lineage-specific robustness analysis.
```

---

## 31. Final Conceptual Summary

This project transforms large-scale mutation–dependency associations into biologically interpretable and therapeutically meaningful synthetic dependency candidates.

The central framework is:

```text
Mutation matrix
      +
CRISPR dependency matrix
      +
cell line metadata
      ↓
Regression-based mutation–dependency association
      ↓
Candidate mutation-conditioned dependencies
      ↓
Network-constrained evidence integration
      ├── STRING functional proximity
      ├── Reactome pathway concordance
      ├── OmniPath directional signaling
      ├── CORUM protein complex support
      ├── SynLethDB known SL prior
      ├── DepMap co-dependency
      └── Druggability evidence
      ↓
Weighted synthetic dependency score
      ↓
Prioritized mutation–target pairs
      ↓
Mechanistic interpretation + experimental validation
```

Best project positioning:

> A network-constrained, empirical evidence-weighted framework for prioritizing mutation-conditioned synthetic dependencies from CRISPR gene dependency screens.

