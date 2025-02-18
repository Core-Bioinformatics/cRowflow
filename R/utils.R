#' Get clustering labels (helper function)
#'
#' @param data The dataset to cluster.
#' @param clustering_algo A function that returns cluster labels (vector).
#' @param algo_params A list of parameters passed to `clustering_algo`.
#' @param labels_name A string with the name of the parameter indicating where the cluster labels are stored in
#'   the output of \code{clustering_algo}.
#'   
#' @return Integer or numeric vector of cluster labels.
get_clustering_labels <- function(data, clustering_algo, algo_params, labels_name) {

  if (!is.null(labels_name)) {
    labels <- do.call(clustering_algo, c(list(data), algo_params))[[labels_name]]
  } else {
    labels <- do.call(clustering_algo, c(list(data), algo_params))

    if (!is.vector(labels) || length(labels) != nrow(data)) {
      stop("Clustering function must return a vector of cluster labels, one per row of data.")
    }
  }

  labels
}

#' Reconcile partitions and derive majority-voting labels (helper function)
#'
#' @param partitions_list A list of integer vectors (cluster labels).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{majority_voting_labels}: integer vector
#'   \item \code{reconciled_partitions}: list of integer vectors after label alignment
#' }
#'
reconcile_partitions_and_majority_voting <- function(partitions_list) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("The 'clue' package is required for reconcile_partitions_and_majority_voting.")
  }
  if (length(partitions_list) == 0) {
    stop("No partitions provided.")
  }

  # Find the most frequent partition
  partition_strings <- vapply(partitions_list, function(x)
    paste(x, collapse = "_"), FUN.VALUE = character(1))
  freq_table <- sort(table(partition_strings), decreasing = TRUE)
  sorted_partitions <- names(freq_table)

  reference_str <- sorted_partitions[1]
  reference_partition <- as.integer(strsplit(reference_str, "_")[[1]])
  unique_ref <- unique(reference_partition)

  # Store the final mapping for each partition_string -> "aligned" vector
  mapping_dict <- list()
  mapping_dict[[reference_str]] <- reference_partition

  reconciled_partitions <- list(reference_partition)

  # For each of the other partitions, align them to the reference
  for (current_str in sorted_partitions[-1]) {
    current_partition <- as.integer(strsplit(current_str, "_")[[1]])
    unique_current <- unique(current_partition)

    # Build cost matrix using 1 - Jaccard
    cost_matrix <- matrix(0,
                          nrow = length(unique_current),
                          ncol = length(unique_ref))

    for (i in seq_along(unique_current)) {
      c1 <- unique_current[i]
      mask1 <- (current_partition == c1)
      sum1 <- sum(mask1)
      for (j in seq_along(unique_ref)) {
        c2 <- unique_ref[j]
        mask2 <- (reference_partition == c2)
        intersection <- sum(mask1 & mask2)
        union_ <- sum1 + sum(mask2) - intersection
        jaccard_index <- if (union_ > 0)
          intersection / union_
        else
          0
        cost_matrix[i, j] <- 1 - jaccard_index
      }
    }

    # Solve assignment using clue::solve_LSAP
    assignment <- clue::solve_LSAP(cost_matrix)
    # assignment[i] = column chosen for row i
    mapping_from <- unique_current
    mapping_to <- unique_ref[assignment]

    # Remap current_partition
    adjusted_partition <- current_partition
    for (k in seq_along(mapping_from)) {
      adjusted_partition[current_partition == mapping_from[k]] <- mapping_to[k]
    }

    mapping_dict[[current_str]] <- adjusted_partition
    reconciled_partitions[[length(reconciled_partitions) + 1]] <- adjusted_partition
  }

  # Majority voting across the aligned partitions
  n_samples <- length(reference_partition)
  n_runs <- length(partitions_list)
  stacked_partitions <- matrix(0, nrow = n_samples, ncol = n_runs)

  for (i in seq_along(partitions_list)) {
    p_str <- partition_strings[i]
    aligned <- mapping_dict[[p_str]]
    stacked_partitions[, i] <- aligned
  }

  # For each sample, find the most common label
  majority_voting_labels <- apply(stacked_partitions, 1, function(sample_labels) {
    tab <- table(sample_labels)
    as.integer(names(which.max(tab)))
  })

  list(majority_voting_labels = majority_voting_labels,
       reconciled_partitions = reconciled_partitions)
}
