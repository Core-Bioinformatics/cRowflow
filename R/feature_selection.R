#' Genetic Algorithm Feature Selector
#'
#' Selects the most stable subset of features using a genetic algorithm.
#'
#' This function uses a genetic algorithm to iteratively optimize feature selection for clustering stability.
#' It repeatedly applies a stochastic clustering algorithm on different feature subsets and evaluates their
#' stability using Element-Centric Consistency (ECC). The algorithm evolves through selection, crossover, and
#' mutation, converging on the feature set that maximizes clustering robustness.
#'
#' @param data A data frame or matrix of shape \code{(n_samples, n_features)}.
#' @param clustering_algo A clustering function or callable (e.g., \code{kmeans}) that is stochastic and
#'   compatible with \code{\link{stochastic_clustering_runner}}.
#' @param labels_name A string with the name of the parameter indicating where the cluster labels are stored in
#'   the output of \code{clustering_algo}.
#' @param n_runs Integer, number of repeated clustering runs per feature subset. Defaults to 30.
#' @param population_size Integer, number of candidate feature subsets in each generation. Defaults to 20.
#' @param generations Integer, number of generations to evolve. Defaults to 50.
#' @param tournament_selection_k Integer, tournament size for parent selection. Defaults to 3.
#' @param mutation_rate Numeric, probability of flipping a feature selection bit in an offspring. Defaults to 0.01.
#' @param crossover_rate Numeric, probability of performing crossover between two parents. Defaults to 0.8.
#' @param elite_size Integer, number of top-performing subsets preserved per generation. Defaults to 2.
#' @param n_generations_no_change Integer, early stopping criterion: generations without improvement. Defaults to 10.
#' @param verbose Logical; if \code{TRUE}, prints progress messages for each generation. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed on to \code{clustering_algo}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{best_features}}{A character vector (if \code{data} is a data frame) or integer indices (if \code{data}
#'         is a matrix) representing the optimal feature subset found.}
#'   \item{\code{best_ecc}}{Numeric, the highest median ECC achieved during evolution.}
#'   \item{\code{history}}{A data frame tracking the best fitness (median ECC) by generation.}
#'   \item{\code{best_fitness_scr_result}}{Clustering results for the optimal feature subset.}
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#'   ga_results <- genetic_algorithm_feature_selector(
#'     data = data,
#'     clustering_algo = kmeans,
#'     labels_name = "cluster",
#'     n_runs = 30,
#'     population_size = 20,
#'     generations = 50,
#'     tournament_selection_k = 3,
#'     mutation_rate = 0.01,
#'     crossover_rate = 0.8,
#'     elite_size = 2,
#'     n_generations_no_change = 10,
#'     verbose = TRUE,
#'     centers = 3
#'   )
#'   print("Best Features:")
#'   print(ga_results$best_features)
#' }
#'
#' @export
genetic_algorithm_feature_selector <- function(
    data,
    clustering_algo,
    labels_name,
    n_runs = 30,
    population_size = 20,
    generations = 50,
    tournament_selection_k = 3,
    mutation_rate = 0.01,
    crossover_rate = 0.8,
    elite_size = 2,
    n_generations_no_change = 10,
    verbose = FALSE,
    ...
) {

  if (!is.function(clustering_algo)) {
    stop("The provided 'clustering_algo' must be a function.")
  }

  feature_names <- NULL
  best_fitness_scr_result <- NULL

  if (is.data.frame(data)) {
    feature_names <- colnames(data)
    data <- as.matrix(data)
  }
  n_samples <- nrow(data)
  n_features <- ncol(data)
  if (n_features < 1) {
    stop("No features found in 'data'.")
  }

  evaluate_fitness <- function(candidate) {
    chosen <- which(candidate == 1L)
    if (length(chosen) == 0L) {
      return(list(scr_result = NULL, median_ecc = -Inf))
    }
    subset_data <- data[, chosen, drop = FALSE]

    results <- stochastic_clustering_runner(
      data = subset_data,
      clustering_algo = clustering_algo,
      labels_name = labels_name,
      n_runs = n_runs,
      verbose = FALSE,
      ...
    )
    list(scr_result = results, median_ecc = median(results$ecc))
  }


  population <- matrix(
    sample(c(0L, 1L), population_size * n_features, replace = TRUE),
    nrow = population_size,
    ncol = n_features
  )

  scr_results_list <- vector("list", population_size)
  fitness_scores <- numeric(population_size)

  for (i in seq_len(population_size)) {
    fitness_score_results <- evaluate_fitness(population[i, ])
    fitness_scores[i] <- fitness_score_results$median_ecc
    scr_results_list[[i]] <- fitness_score_results$scr_result
  }

  best_idx <- which.max(fitness_scores)
  best_fitness <- fitness_scores[best_idx]
  best_candidate <- population[best_idx, , drop = FALSE]
  best_fitness_scr_result <- scr_results_list[[best_idx]]

  history_list <- list(
    list(generation = 0L, best_fitness = best_fitness)
  )
  if (verbose) {
    message(sprintf("Gen 0 - Best ECC: %.4f", best_fitness))
  }

  iters_no_improvement <- 0L

  for (gen in seq_len(generations)) {
    if (iters_no_improvement >= n_generations_no_change) {
      break
    }

    tournament_selection <- function(k) {
      selected <- sample.int(population_size, size = k, replace = FALSE)
      best_sel <- selected[which.max(fitness_scores[selected])]
      population[best_sel, ]
    }

    new_population <- list()
    elite_indices <- order(fitness_scores, decreasing = TRUE)[seq_len(elite_size)]
    for (idx in elite_indices) {
      new_population[[length(new_population) + 1L]] <- population[idx, ]
    }

    while (length(new_population) < population_size) {
      parent1 <- tournament_selection(tournament_selection_k)
      parent2 <- tournament_selection(tournament_selection_k)

      if (runif(1) < crossover_rate) {
        point <- sample.int(n_features - 1, size = 1)
        offspring1 <- c(parent1[1:point], parent2[(point + 1):n_features])
        offspring2 <- c(parent2[1:point], parent1[(point + 1):n_features])
      } else {
        offspring1 <- parent1
        offspring2 <- parent2
      }

      mutate_one <- function(offspring) {
        for (i in seq_len(n_features)) {
          if (runif(1) < mutation_rate) {
            offspring[i] <- 1L - offspring[i]
          }
        }
        offspring
      }
      offspring1 <- mutate_one(offspring1)
      offspring2 <- mutate_one(offspring2)

      new_population[[length(new_population) + 1L]] <- offspring1
      if (length(new_population) < population_size) {
        new_population[[length(new_population) + 1L]] <- offspring2
      }
    }

    population <- do.call(rbind, new_population)

    generation_scr_results_list <- vector("list", population_size)
    fitness_scores <- numeric(population_size)

    for (i in seq_len(population_size)) {
      fitness_score_results <- evaluate_fitness(population[i, ])
      fitness_scores[i] <- fitness_score_results$median_ecc
      generation_scr_results_list[[i]] <- fitness_score_results$scr_result
    }

    generation_best_idx <- which.max(fitness_scores)
    generation_best_fitness <- fitness_scores[generation_best_idx]

    if (generation_best_fitness > best_fitness) {
      best_fitness <- generation_best_fitness
      best_candidate <- population[generation_best_idx, , drop = FALSE]
      iters_no_improvement <- 0L
      best_fitness_scr_result <- generation_scr_results_list[[generation_best_idx]]
    } else {
      iters_no_improvement <- iters_no_improvement + 1L
    }

    history_list[[length(history_list) + 1L]] <- list(
      generation = gen,
      best_fitness = best_fitness
    )
    if (verbose) {
      message(sprintf("Gen %d - Best ECC: %.4f", gen, best_fitness))
    }
  }

  history_df <- do.call(
    rbind,
    lapply(history_list, function(x) as.data.frame(x, stringsAsFactors = FALSE))
  )

  best_candidate_vec <- as.integer(best_candidate[1, ])
  chosen_features <- which(best_candidate_vec == 1L)

  if (!is.null(feature_names)) {
    chosen_features <- feature_names[chosen_features]
  }

  return(list(
    best_features = chosen_features,
    best_ecc = best_fitness,
    history = history_df,
    best_fitness_scr_result = best_fitness_scr_result
  ))
}







