library(testthat)
library(cRowflow)  # Replace with your package's name

# Load required packages for the clustering algorithms used in tests
if (!requireNamespace("e1071", quietly = TRUE)) {
  skip("e1071 not installed, skipping tests")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  skip("igraph not installed, skipping tests")
}

# Define a helper for the Louvain clustering algorithm.
louvain_clustering <- function(data, ...) {
  # Compute Euclidean distance matrix
  dmat <- as.matrix(dist(data))
  # Create a similarity matrix: higher similarity for closer points.
  sim_mat <- 1 / (1 + dmat)
  diag(sim_mat) <- 0
  # Build an undirected weighted graph
  g <- igraph::graph_from_adjacency_matrix(sim_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  # Apply the Louvain algorithm. It returns a membership vector.
  cl <- igraph::cluster_louvain(g)
  return(cl$membership)
}


set.seed(123)
test_data <- matrix(rnorm(250), nrow = 50, ncol = 5)

algorithms <- list(
  kmeans = list(
    func = kmeans,
    labels_name = "cluster",         # kmeans returns a list with the element 'cluster'
    extra_args = list(centers = 3, nstart = 1)  # nstart=1 preserves randomness
  ),
  cmeans = list(
    func = e1071::cmeans,
    labels_name = "cluster",         # cmeans returns an object with a 'cluster' component
    extra_args = list(centers = 3, iter.max = 100, m = 2)
  ),
  louvain = list(
    func = louvain_clustering,
    labels_name = NULL,              # louvain_clustering (my helper function) returns membership vector directly
    extra_args = list()
  )
)

test_that("stochastic_clustering_runner works with stochastic algorithms", {
  for (algo_name in names(algorithms)) {
    algo_info <- algorithms[[algo_name]]
    clustering_algo <- algo_info$func
    lbl_name <- algo_info$labels_name
    extra_args <- algo_info$extra_args
    
    message("Testing algorithm: ", algo_name)
    
    # Build argument list for the runner, including extra arguments.
    args_list <- c(
      list(
        data = test_data,
        clustering_algo = clustering_algo,
        labels_name = lbl_name,
        n_runs = 10,
        verbose = FALSE
      ),
      extra_args
    )
    
    # Call the runner with the given algorithm using do.call.
    result <- do.call(stochastic_clustering_runner, args_list)
    
    print(summary(result$ecc))
    # Check that the result is a list with the expected components.
    expected_components <- c("partitions", "partition_frequencies", 
                             "majority_voting_labels", "ecc", "seeds")
    expect_true(is.list(result))
    expect_true(all(expected_components %in% names(result)))
    
    # Verify that the number of partitions equals n_runs (5).
    expect_equal(length(result$partitions), 10)
    
    # Verify that the seeds are generated as expected.
    expect_equal(result$seeds, seq(100, by = 100, length.out = 10))
    
    # Ensure each partition is an integer vector with length equal to the number of observations.
    lapply(result$partitions, function(partition) {
      expect_true(is.integer(partition))
      expect_equal(length(partition), nrow(test_data))
    })
    
    # Check that the majority voting labels and ECC scores match the number of rows.
    expect_equal(length(result$majority_voting_labels), nrow(test_data))
    expect_equal(length(result$ecc), nrow(test_data))
  }
})
