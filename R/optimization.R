#' Parameter Optimizer for Stochastic Clustering
#'
#' Optimizes individual hyperparameters of a stochastic clustering algorithm.
#'
#' This function systematically tunes each hyperparameter separately by performing repeated clustering
#' and evaluating stability using Element-Centric Consistency (ECC). It's purpose is to help identify the
#' value of each parameter that maximizes clustering robustness.
#'
#' @param data  A data frame or matrix of shape \code{(n_samples, n_features)}.
#' @param clustering_algo A clustering function or callable (e.g., \code{kmeans}) that is stochastic and
#'   compatible with \code{\link{stochastic_clustering_runner}}.
#' @param labels_name  A string with the name of the parameter indicating where the cluster labels are stored in
#'   the output of \code{clustering_algo}.
#' @param parameters_optimise_list A named list specifying the parameters and their candidate values, e.g.,
#'   \code{list(n_clusters = c(2, 3, 4))}.
#' @param n_runs Integer; number of clustering runs per parameter setting (default is 30).
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' @param ... Additional parameters to be passed to \code{clustering_algo} (beyond the one being tuned).
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{results_df}}{A data frame where each row represents a parameter-value combination.
#'      Each row includes:
#'      \itemize{
#'        \item \code{param}: a string of the form \code{"parameterName_value"},
#'        \item \code{ecc}: the numeric vector of ECC values from each run,
#'        \item \code{median_ecc}: the median ECC across runs.
#'      }}
#'   \item{\code{stochastic_clustering_results}}{A named list containing the full output from
#'      \code{stochastic_clustering_runner} for each parameter setting tested, keyed by \code{"parameterName_value"}.}
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#'   results <- parameter_optimizer(
#'     data = data,
#'     clustering_algo = kmeans,
#'     labels_name = "cluster",
#'     parameters_optimise_list = list(centers = 2:4),
#'     n_runs = 30,
#'     verbose = TRUE
#'   )
#'   results_df <- results$results_df
#'   clustering_details <- results$stochastic_clustering_results
#'   print(results_df)
#' }
#'
#' @export
parameter_optimizer <- function(data,
                                clustering_algo,
                                labels_name,
                                parameters_optimise_list,
                                n_runs = 30,
                                verbose = FALSE,
                                ...) {
  if (!is.function(clustering_algo)) {
    stop("The provided clustering_algo must be an R function.")
  }


  results_list <- list()
  row_list <- list()

  # Loop over each parameter, then each candidate value
  for (param_name in names(parameters_optimise_list)) {
    values_to_try <- parameters_optimise_list[[param_name]]
    for (val in values_to_try) {
      if (verbose) {
        message(sprintf("Running with %s = %s", param_name, val))
      }

      arg_list <- list(
        data = data,
        clustering_algo = clustering_algo,
        labels_name = labels_name,
        n_runs = n_runs,
        verbose = FALSE
      )

      # Add user arguments
      dot_args <- list(...)
      if (length(dot_args) > 0) {
        arg_list <- c(arg_list, dot_args)
      }

      arg_list[[param_name]] <- val
      runner_result <- do.call(stochastic_clustering_runner, arg_list)

      ecc_vec <- runner_result$ecc
      median_ecc <- median(ecc_vec)

      if (verbose) {
        message(sprintf("   Median ECC: %f", median_ecc))
        message("--------------------------------------------------------------")
      }

      current_row <- list(
        param = paste0(param_name, "_", val),
        ecc = ecc_vec,
        median_ecc = median_ecc
      )

      current_row[[param_name]] = val
      row_list[[length(row_list) + 1]] <- current_row

      results_list[[paste0(param_name, "_", val)]] <- runner_result
    }
  }

  # Create results dataframe
  results_df <- do.call(rbind, lapply(row_list, function(x) {
    row_df <- as.data.frame(x[names(x) != "ecc"])
    row_df
  }))

  results_df$ecc <- I(lapply(row_list, `[[`, "ecc"))

  list(
    results_df = results_df,
    stochastic_clustering_results = results_list
  )

}

#' Full Grid Search for Stochastic Clustering Parameters
#'
#' Performs a full grid search over multiple clustering hyperparameters.
#'
#' This function evaluates all possible combinations of specified parameters by running repeated stochastic clustering
#' for each combination and computing Element-Centric Consistency (ECC) for each setting. Its purpose is to help identify
#' the optimal parameter set that maximizes clustering stability.
#'
#' @param data A data frame or matrix of shape \code{(n_samples, n_features)} to cluster.
#' @param clustering_algo A clustering function or callable (e.g., \code{kmeans}) that is stochastic and compatible with
#'   \code{\link{stochastic_clustering_runner}}.
#' @param labels_name A string with the name of the parameter indicating where the cluster labels are stored in
#'   the output of \code{clustering_algo}.
#' @param param_grid A named list of the form \code{list(paramA = c(values), paramB = c(values), ...)} specifying all parameter
#'   values to be combined.
#' @param n_runs Integer; number of repeated runs per parameter combination (default is 30).
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' @param ... Additional parameters passed to \code{clustering_algo}.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{results_df}}{A data frame where each row corresponds to one parameter combination.
#'      Each row includes:
#'      \itemize{
#'        \item \code{params}: a (list-)column with a named list of parameter values,
#'        \item \code{ecc}: the numeric vector of ECC values for that combination,
#'        \item \code{median_ecc}: the median ECC.
#'      }}
#'   \item{\code{stochastic_clustering_results}}{A named list of the full outputs from
#'      \code{stochastic_clustering_runner} for each parameter combination, keyed by a string such as
#'      \code{"paramA_valA_paramB_valB"}.}
#' }
#'
#' @examples
#' \dontrun{
#'    set.seed(42)
#'    data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#'    param_grid <- list(
#'     centers = 2:4,
#'     algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")
#'    )
#'    results <- parameter_searcher(
#'     data = data,
#'     clustering_algo = kmeans,
#'     labels_name = "cluster",
#'     param_grid = param_grid,
#'     n_runs = 30,
#'     verbose = TRUE
#'     )
#'     results_df <- results$results_df
#'     clustering_details <- results$stochastic_clustering_results
#' }
#'
#' @export
parameter_searcher <- function(data,
                               clustering_algo,
                               labels_name,
                               param_grid,
                               n_runs = 30,
                               verbose = FALSE,
                               ...) {

  if (!is.function(clustering_algo)) {
    stop("The provided clustering_algo must be an R function.")
  }

  # Create all combinations of parameters via expand.grid
  combinations_df <- do.call(expand.grid, c(param_grid, list(stringsAsFactors = FALSE)))

  full_results_list <- list()
  row_list <- list()

  # Create all combinations of parameters via expand.grid
  for (i in seq_len(nrow(combinations_df))) {
    param_set <- combinations_df[i, , drop = FALSE]

    # Convert it into a proper named list
    param_dict <- as.list(param_set)

    if (verbose) {
      message("Testing parameters: ", paste0(
        paste(names(param_dict), param_dict, sep="="), collapse=", "
      ))
    }

    # Call the stochastic clustering runner with the current combination
    arg_list <- list(
      data = data,
      clustering_algo = clustering_algo,
      labels_name = labels_name,
      n_runs = n_runs,
      verbose = FALSE
    )

    dot_args <- list(...)
    if (length(dot_args) > 0) {
      arg_list <- c(arg_list, dot_args)
    }

    arg_list <- c(arg_list, param_dict)
    runner_result <- do.call(stochastic_clustering_runner, arg_list)

    ecc_vec <- runner_result$ecc
    median_ecc <- median(ecc_vec)

    if (verbose) {
      message(sprintf("   Median ECC: %f", median_ecc))
      message("--------------------------------------------------------------")
    }

    # Build a string key for storing results
    param_str_parts <- mapply(
      function(k, v) paste0(k, "_", v),
      names(param_dict), param_dict
    )
    param_str <- paste(param_str_parts, collapse = "_")

    full_results_list[[param_str]] <- runner_result

    # Store a row for the final data frame
    current_row <- list(
      params = param_dict,
      ecc = ecc_vec,
      median_ecc = median_ecc
    )

    for (param_name in names(param_dict)) {
      current_row[[param_name]] <- param_dict[[param_name]]
    }

    row_list[[length(row_list) + 1]] <- current_row
  }

  results_df <- do.call(rbind, lapply(row_list, function(x) {
    row_df <- as.data.frame(x[names(x) != "ecc"])
    row_df
  }))

  results_df$params <- I(lapply(row_list, `[[`, "params"))
  results_df$ecc <- I(lapply(row_list, `[[`, "ecc"))

  results_cols <- c("params", "ecc", names(results_df)[!(names(results_df) %in% c("params", "ecc"))])
  results_df <- results_df[, results_cols]

  list(
    results_df = results_df,
    stochastic_clustering_results = full_results_list
  )
}

