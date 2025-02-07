#' Stochastic Clustering Runner
#'
#' Repeats stochastic clustering multiple times to assess stability.
#'
#' This function runs a clustering algorithm multiple times with different random seeds
#' and evaluates the stability of results using Element-Centric Consistency (ECC).
#' It identifies in an element-wise precision the stability of clustering results and
#' provides a consensus labeling through majority voting.
#'
#' @param data A data matrix or data frame on which clustering is performed.
#' @param clustering_algo A stochastic clustering function that returns a vector of cluster labels.
#'   The function should produce different results when the seed is changed (via \code{set.seed}).
#' @param labels_name (Optional) A string indicating the name of the element in the output of
#'   \code{clustering_algo} that contains the cluster labels. If \code{NULL}, the function is assumed
#'   to return the labels directly.
#' @param n_runs Number of times to repeat the clustering process. Defaults to 30.
#' @param verbose Logical; if \code{TRUE}, prints progress messages showing the seed used in each run.
#' @param ... Additional arguments to be passed to \code{clustering_algo}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{partitions}}{A list of integer vectors representing clustering assignments for each run.}
#'   \item{\code{partition_frequencies}}{A list mapping each unique partition (as a concatenated string)
#'         to its frequency of occurrence.}
#'   \item{\code{majority_voting_labels}}{An integer vector of consensus cluster labels derived from majority voting.}
#'   \item{\code{ecc}}{A numeric vector of element-centric consistency (ECC) scores indicating clustering stability.}
#'   \item{\code{seeds}}{An integer vector of the random seeds used in each run.}
#' }
#'
#' @details
#' The function generates seeds starting from 100 and increasing by 100 for each run. For each seed, it sets
#' the random number generator using \code{set.seed} to ensure reproducibility of the stochastic clustering.
#' The clustering results from each run are reconciled through majority voting and evaluated using ECC.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#'   results <- stochastic_clustering_runner(
#'     data = data,
#'     clustering_algo = kmeans,
#'     labels_name = "cluster",
#'     n_runs = 30,
#'     verbose = TRUE,
#'     centers = 3
#'   )
#'   print("Majority Voting Labels:")
#'   print(results$majority_voting_labels)
#'   print("ECC:")
#'   print(results$ecc)
#' }
#'
#' @export
stochastic_clustering_runner <- function(data,
                                         clustering_algo,
                                         labels_name = NULL,
                                         n_runs = 30,
                                         verbose = FALSE,
                                         ...) {
  # Validate clustering_algo is callable
  if (!is.function(clustering_algo)) {
    stop("clustering_algo must be an R function that returns cluster labels.")
  }

  seeds <- seq(100, by = 100, length.out = n_runs)

  partitions_list <- vector("list", n_runs)
  for (i in seq_along(seeds)) {
    seed_value <- seeds[i]
    if (verbose) {
      message("  Seed = ", seed_value)
    }

    algo_params <- list(...)

    # Set seed before running clustering algorithm.
    set.seed(seed_value)

    # Get cluster labels from the helper.
    labels <- get_clustering_labels(data, clustering_algo, algo_params, labels_name)
    labels <- as.integer(labels)

    partitions_list[[i]] <- labels
  }

  # Convert each partition to a string and check frequency using table.
  partition_strings <- vapply(
    partitions_list,
    FUN = function(x)
      paste(x, collapse = "_"),
    FUN.VALUE = character(1)
  )
  partition_table <- table(partition_strings)
  partition_frequencies <- as.list(partition_table)


  # Reconcile partitions and majority voting
  majority_res <- reconcile_partitions_and_majority_voting(partitions_list)
  majority_voting_labels <- majority_res$majority_voting_labels

  ecc <- ClustAssess::element_consistency(partitions_list)

  list(
    partitions = partitions_list,
    partition_frequencies = partition_frequencies,
    majority_voting_labels = majority_voting_labels,
    ecc = ecc,
    seeds = seeds
  )
}
