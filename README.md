# cRowflow

## Overview
`cRowflow` is an R package designed for assessing clustering stability through repeated stochastic clustering. It is compatible with any clustering algorithm that outputs labels or returns results containing cluster assignments. By running clustering multiple times with different seeds, `cRowflow` quantifies clustering consistency using Element-Centric Similarity (ECS) and Element-Centric Consistency (ECC), offering insights into the robustness and reproducibility of cluster assignments. The package enables users to optimize feature subsets, fine-tune clustering parameters, and evaluate clustering robustness against perturbations.

`cRowflow` generalizes the [`ClustAssess`](https://github.com/Core-Bioinformatics/ClustAssess) package, which focuses on parameter selection for community-detection clustering in single-cell analysis. It extends this approach to any clustering task, enabling a data-driven identification of robust and reproducible clustering solutions across diverse applications.

## Features

### 1. Stochastic Clustering Runner
Runs a stochastic clustering algorithm multiple times with different random seeds and evaluates the stability of results using ECC. It identifies, with element-wise precision, the stability of clustering results and provides majority voting labels.

**Important Note:** In R, the kmeans algorithm stores cluster labels in the "cluster" field, so we must set labels_name = "cluster" to extract them correctly. If the clustering algorithm directly returns the labels, labels_name can be left as NULL. Additionally, we must specify the number of clusters (e.g., centers = 5), as it does not have a default k value.

#### Function: `stochastic_clustering_runner()`
```r
set.seed(42)
my_data <- matrix(rnorm(200 * 10), ncol = 10)

# Run stochastic clustering
result <- stochastic_clustering_runner(
  data = my_data,
  clustering_algo = kmeans,
  labels_name = "cluster",
  n_runs = 30,
  centers = 5
)
```
#### Returns:
- `partitions`: List of clustering assignments for each run
- `majority_voting_labels`: Consensus labels from majority voting
- `ecc`: ECC scores indicating clustering stability

---

### 2. Genetic Algorithm Feature Selector
Uses a genetic algorithm to iteratively optimize feature selection for clustering stability. It repeatedly applies stochastic clustering with different feature subsets and evaluates stability using ECC. The algorithm evolves through selection, crossover, and mutation, converging on the feature set that maximizes clustering robustness.

#### Function: `genetic_algorithm_feature_selector()`
```r
set.seed(42)
my_data <- matrix(rnorm(200 * 10), ncol = 10)

# Run genetic algorithm feature selection
result <- genetic_algorithm_feature_selector(
  data = my_data,
  clustering_algo = kmeans,
  labels_name = "cluster",
  n_runs = 30,
  population_size = 20,
  generations = 50,
  centers = 5
)
```
#### Returns:
- `best_features`: Optimal feature subset
- `best_ecc`: Highest median ECC achieved
- `history`: Evolution of fitness across generations

---

### 3. Parameter Optimizer
Systematically tunes each hyperparameter separately by performing repeated clustering and evaluating stability using ECC.

#### Function: `parameter_optimizer()`
```r
set.seed(42)
my_data <- matrix(rnorm(200 * 10), ncol = 10)

# Run parameter optimization
# find number of clusters (between 2 to 5) resultin in most stable partitions.
result <- parameter_optimizer(
  data = my_data,
  clustering_algo = kmeans,
  labels_name = "cluster",
  parameters_optimise_list = list(centers = 2:5),
  n_runs = 30
)
```
#### Returns:
- `results_df`: ECC scores for each hyperparameter setting
- `stochastic_clustering_results`: Full clustering results for each setting

---

### 4. Parameter Searcher
Evaluates all possible combinations (exhaustive grid search) of specified parameters, running repeated clustering and computing ECC for each combination. The purpose is to find the configuration (set of hyperparameter values) that provides the most stable clustering results.

#### Function: `parameter_searcher()`
```r
set.seed(42)
my_data <- matrix(rnorm(200 * 10), ncol = 10)

# Run exhaustive parameter search
# Find combination of centers and nstart that results in the most stable partitions.
result <- parameter_searcher(
  data = my_data,
  clustering_algo = kmeans,
  labels_name = "cluster",
  param_grid = list(centers = 2:4, nstart = c(1, 5)),
  n_runs = 30
)
```
#### Returns:
- `results_df`: A dataframe of ECC values for different parameter combinations
- `stochastic_clustering_results`: Detailed clustering results

---

### 5. K-Fold Clustering Validator
Evaluates how stable clustering assignments remain across different data partitions by comparing clustering results on k-fold subsets with those from the full dataset. ECS is used to quantify similarity/stability between fold-level clustering and the baseline (full dataset).

#### Function: `kfold_clustering_validator()`
```r
# Simulate dataset
set.seed(42)
my_data <- matrix(rnorm(200 * 10), ncol = 10)

# Run k-fold clustering validation
result <- kfold_clustering_validator(
  data = my_data,
  clustering_algo = kmeans,
  labels_name = "cluster",
  k_folds = 5,
  n_runs = 30,
  centers = 5
)
```
#### Returns:
- `baseline_results`: Clustering results on full dataset
- `kfolds_robustness_results`: Clustering robustness metrics for each fold

---

### 6. Perturbation Robustness Tester
Tests how stable clustering results are when features are altered/perturbed. The user must provide a perturbation function, which modifies the dataset before clustering is re-run. Stability is assessed using ECS between the baseline clustering and perturbation-induced clusterings.

#### Function: `perturbation_robustness_tester()`
```r
# Simulate dataset
set.seed(42)
my_data <- matrix(rnorm(200 * 10), ncol = 10)

# Define a feature-shuffling perturbation function
shuffle_features <- function(data) {
  data[sample(nrow(data)), ]
}

# Run perturbation robustness test
result <- perturbation_robustness_tester(
  data = my_data,
  clustering_algo = kmeans,
  labels_name = "cluster",
  perturbation_func = shuffle_features,
  n_perturbations = 10,
  n_runs = 30,
  centers = 5
)
```
#### Returns:
- `perturbation_el_sim_scores`: Similarity scores after perturbation
- `mean_score`: Mean robustness score

---

## Installation
To install the latest version from GitHub, use:

```r
# install.packages("devtools")
devtools::install_github("Core-Bioinformatics/cRowflow")
```

Load the package and necessary dependencies:
```r
library(cRowflow)
```

### Dependencies
- caret
- ggplot2
- dplyr
- tidyr
- viridis
- purrr
- rlang
- ClustAssess

## License
This package is released under the MIT License.

Developed by Rafael Kollyfas (rk720@cam.ac.uk), Core Bioinformatics (Mohorianu Lab) group, University of Cambridge. February 2025.
