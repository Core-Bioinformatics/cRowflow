#' K-Fold Clustering Validator
#'
#' Performs k-fold validation to assess clustering stability.
#'
#' This function first clusters the full dataset to obtain baseline clustering labels.
#' It then performs k-fold cross-validation by splitting the data into k folds, clustering
#' on each training set (k-1 folds), and comparing the resulting majority-voting labels with
#' those obtained from the full dataset. Element-Centric Similarity (ECS) is computed to quantify
#' the consistency between fold-level clustering and the baseline.
#'
#' @param data A data frame or matrix of shape \code{(n_samples, n_features)} to cluster.
#' @param clustering_algo A clustering function or callable (e.g., \code{kmeans}) that is stochastic and
#'   compatible with \code{\link{stochastic_clustering_runner}}.
#' @param labels_name A string with the name of the label containing the clusters returned by the algorithm.
#' @param k_folds Integer; number of folds for cross-validation (default is 5).
#' @param n_runs Integer; number of repeated runs per clustering call (default is 30).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default is \code{FALSE}).
#' @param ... Additional parameters passed to \code{clustering_algo}.
#'
#' @return A named list with two components:
#' \describe{
#'   \item{\code{baseline_results}}{
#'     A list containing the clustering results on the full dataset (e.g., \code{majority_voting_labels},
#'     \code{ecc}, \code{partitions}, etc.).
#'   }
#'   \item{\code{kfolds_robustness_results}}{
#'     A list with one entry per fold. Each entry is a list containing:
#'     \describe{
#'       \item{\code{stochastic_clustering_results}}{
#'         The clustering results on the training subset of the current fold.
#'       }
#'       \item{\code{used_indices}}{
#'         The indices of the data used for clustering in the current fold.
#'       }
#'       \item{\code{leave_out_indices}}{
#'         The indices of the data left out in the current fold.
#'       }
#'       \item{\code{el_score_vector}}{
#'         A numeric vector of element-wise similarity scores comparing the fold's majority-voting labels
#'         with the corresponding baseline labels, as computed by \code{ClustAssess::element_sim_elscore}.
#'       }
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#'   results <- kfold_clustering_validator(
#'     data = data,
#'     clustering_algo = kmeans,
#'     labels_name = "cluster",
#'     k_folds = 5,
#'     n_runs = 30,
#'     verbose = TRUE,
#'     centers = 3
#'   )
#'   baseline <- results$baseline_results
#'   kfold_results <- results$kfolds_robustness_results
#' }
#'
#' @export
kfold_clustering_validator <- function(
    data,
    clustering_algo,
    labels_name,
    k_folds = 5,
    n_runs = 30,
    verbose = FALSE,
    ...
) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is needed for k-fold splitting. Please install it.")
  }
  if (!requireNamespace("ClustAssess", quietly = TRUE)) {
    stop("Package 'ClustAssess' is needed for element_sim_elscore(). Please install it.")
  }

  baseline_results <- stochastic_clustering_runner(
    data = data,
    clustering_algo = clustering_algo,
    labels_name = labels_name,
    n_runs = n_runs,
    verbose = verbose,
    ...
  )
  baseline_majority_labels <- baseline_results$majority_voting_labels

  # Create k-fold splits via caret
  folds <- caret::createFolds(1:nrow(data), k = k_folds, returnTrain = TRUE)
  kfolds_robustness_results <- vector("list", length(folds))
  names(kfolds_robustness_results) <- paste0("fold_", seq_along(folds))

  # Loop over each fold
  for (i in seq_along(folds)) {
    if (verbose) {
      message(sprintf("Fold %d:", i))
    }
    used_index <- folds[[i]]
    leave_out_index <- setdiff(seq_len(nrow(data)), used_index)

    # Subset the data for current fold
    used_data <- data[used_index, , drop = FALSE]

    # Cluster only on the "used_data" subset
    fold_result <- stochastic_clustering_runner(
      data = used_data,
      clustering_algo = clustering_algo,
      labels_name = labels_name,
      n_runs = n_runs,
      verbose = verbose,
      ...
    )

    # Compare with baseline majority labels on the same indices
    labels_intersection <- baseline_majority_labels[used_index]
    fold_labels <- fold_result$majority_voting_labels

    # Compute element-wise label similarity
    el_score_vector <- ClustAssess::element_sim_elscore(
      labels_intersection,
      fold_labels
    )

    kfolds_robustness_results[[i]] <- list(
      stochastic_clustering_results = fold_result,
      used_indices = used_index,
      leave_out_indices = leave_out_index,
      el_score_vector = el_score_vector
    )
  }

  list(
    baseline_results = baseline_results,
    kfolds_robustness_results = kfolds_robustness_results
  )
}


#' Perturbation Robustness Tester
#'
#' Evaluates the robustness of clustering solutions under feature perturbations.
#'
#' This function first runs a baseline clustering on the unperturbed data, then repeatedly perturbs the data
#' using a user-supplied perturbation function and performs repeated stochastic clustering on each perturbed dataset.
#' Element-Centric Similarity (ECS) is computed to compare the majority-voting labels of the perturbed clusterings
#' against the baseline, quantifying overall robustness.
#'
#' @param data A data frame or matrix of shape \code{(n_samples, n_features)} representing the dataset to cluster.
#' @param clustering_algo A clustering function or callable (e.g., \code{kmeans}) that is stochastic and compatible with
#'   \code{\link{stochastic_clustering_runner}}.
#' @param labels_name A string with the name of the label containing the clusters returned by the algorithm.
#' @param perturbation_func A function that takes \code{data} and returns a perturbed version of the dataset
#'   (of the same shape/type).
#' @param n_perturbations Integer; number of perturbation trials (default is 10).
#' @param n_runs Integer; number of repeated clustering runs per perturbation (default is 30).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default is \code{FALSE}).
#' @param ... Additional parameters passed to \code{clustering_algo} via \code{stochastic_clustering_runner}.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{baseline_majority_labels}}{The majority-voting labels from the baseline clustering.}
#'   \item{\code{baseline_ecc}}{The Element-Centric Consistency (ECC) score of the baseline clustering.}
#'   \item{\code{perturbation_el_sim_scores}}{A numeric vector of element-wise similarity scores comparing the baseline
#'         labels to the labels obtained after each perturbation.}
#'   \item{\code{mean_score}}{The mean similarity score across all perturbations, indicating overall robustness.}
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   data <- as.data.frame(matrix(rnorm(2000), nrow = 200, ncol = 10))
#'   colnames(data) <- paste0("feature_", 1:10)
#'   shuffle_features <- function(X) {
#'     X_shuffled <- X
#'     for (col in colnames(X_shuffled)) {
#'       X_shuffled[[col]] <- sample(X_shuffled[[col]])
#'     }
#'     X_shuffled
#'   }
#'   perturb_results <- perturbation_robustness_tester(
#'     data = data,
#'     clustering_algo = kmeans,
#'     labels_name = "cluster",
#'     perturbation_func = shuffle_features,
#'     n_perturbations = 5,
#'     n_runs = 30,
#'     verbose = TRUE,
#'     centers = 3
#'   )
#' }
#'
#' @export
perturbation_robustness_tester <- function(
    data,
    clustering_algo,
    labels_name,
    perturbation_func,
    n_perturbations = 10,
    n_runs = 30,
    verbose = FALSE,
    ...
) {
  if (!requireNamespace("ClustAssess", quietly = TRUE)) {
    stop("Package 'clustassess' is needed for `element_sim()`. Please install it.")
  }

  runner <- function(data) {
    stochastic_clustering_runner(
      data = data,
      clustering_algo = clustering_algo,
      labels_name = labels_name,
      n_runs = n_runs,
      verbose = verbose,
      ...
    )
  }

  # Baseline clustering
  baseline_results <- runner(data)
  baseline_majority_labels <- baseline_results$majority_voting_labels
  baseline_ecc <- baseline_results$ecc

  # For each perturbation, cluster the perturbed data and compare to baseline
  perturbation_el_sim_scores <- numeric(n_perturbations)

  for (i in seq_len(n_perturbations)) {
    if (verbose) {
      message(sprintf("Perturbation %d / %d", i, n_perturbations))
    }

    set.seed(i)
    perturbed_data <- perturbation_func(data)
    perturbed_results <- runner(perturbed_data)
    perturbed_majority_labels <- perturbed_results$majority_voting_labels

    # Calculate mean element centric similarity score.
    similarity_score <- ClustAssess::element_sim(
      baseline_majority_labels,
      perturbed_majority_labels
    )
    perturbation_el_sim_scores[i] <- similarity_score
  }

  # Mean similarity across all perturbations
  mean_score <- mean(perturbation_el_sim_scores, na.rm = TRUE)

  list(
    baseline_majority_labels = baseline_majority_labels,
    baseline_ecc = baseline_ecc,
    perturbation_el_sim_scores = perturbation_el_sim_scores,
    mean_score = mean_score
  )
}

