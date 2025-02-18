#' Scatter Plot - ECC
#'
#' Generate a scatter plot with color-coded points. Intended to be used for visualizing the
#' stability of individual points using element-centric consistency (i.e., where the color
#' column represents ECC).
#'
#' @param df A data frame containing the input data.
#' @param x_col A character string specifying the name of the column to plot on the x-axis.
#' @param y_col A character string specifying the name of the column to plot on the y-axis.
#' @param color_col A character string specifying the name of the column to use for coloring the points.
#' @param title A character string specifying the title of the plot (default: "Scatter Plot - ECC").
#' @param alpha A numeric value specifying the transparency level for the points (default: 0.7).
#' @param size A numeric value specifying the size of the points (default: 2).
#' @param cmap A character string specifying the color map used for the points (default: "viridis").
#'
#' @return A ggplot object representing the generated scatter plot.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#'   results <- stochastic_clustering_runner(
#'     data = data,
#'     clustering_algo = kmeans,
#'     labels_name = 'cluster',
#'     n_runs = 30,
#'     verbose = TRUE,
#'     centers = 3
#'   )
#'   library(umap)
#'   umap_res <- umap::umap(data)
#'   umap_df <- as.data.frame(umap_res$layout)
#'   colnames(umap_df) <- c("UMAP_1", "UMAP_2")
#'   umap_df$ecc <- results$ecc
#'   plot_scatter_ecc(umap_df, "UMAP_1", "UMAP_2", "ecc")
#' }
#'
#' @export
plot_scatter_ecc <- function(df,
                             x_col,
                             y_col,
                             color_col,
                             title = "Scatter Plot - ECC",
                             alpha = 0.7,
                             size = 2,
                             cmap = "viridis") {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_col]])) +
    ggplot2::geom_point(size = size, alpha = alpha) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = x_col, y = y_col, color = color_col) +
    ggplot2::scale_color_viridis_c(option = cmap) +
    ggplot2::ggtitle(title)

  return(p)
}

#' Plot Optimization of Parameter (Boxplot)
#'
#' Generate a boxplot visualization of element-centric consistency (ECC) values for different parameter settings.
#'
#' Expected input is a data frame with at least two columns:
#' \describe{
#'   \item{\code{param}}{A string representing the parameter and its value (e.g., "n_clusters_2") tested during optimization.}
#'   \item{\code{ecc}}{A list-column where each entry is a numeric vector of ECC values corresponding to that parameter setting.}
#' }
#'
#' @param parameter_optimizer_results_df A data frame with columns "param" and "ecc".
#' @param name_parameter A character string for the x-axis label (default: "param").
#' @param title A character string for the plot title (default: "Parameter Optimization - Boxplot Distribution").
#'
#' @return A ggplot object representing the boxplot visualization.
#
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   # Generate synthetic data: 200 samples with 10 features
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
#'   plot_optimization_of_parameter_boxplot(results_df)
#' }
#' @export
plot_optimization_of_parameter_boxplot <- function(parameter_optimizer_results_df,
                                                   name_parameter = "param",
                                                   title = "Parameter Optimization - Boxplot Distribution") {
  # Check that required columns exist
  if (!("param" %in% names(parameter_optimizer_results_df)))
    stop("The input data frame must contain a 'param' column.")
  if (!("ecc" %in% names(parameter_optimizer_results_df)))
    stop("The input data frame must contain an 'ecc' column.")

  # Convert the list-column 'ecc' into long format using tidyr::unnest
  df_long <- tidyr::unnest(parameter_optimizer_results_df, cols = c(ecc))

  # Ensure the ECC column is numeric
  df_long$ecc <- as.numeric(df_long$ecc)

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = param, y = ecc)) +
    ggplot2::geom_boxplot(fill = "lightblue", color = "black") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = name_parameter, y = "ECC") +
    ggplot2::ggtitle(title)

  return(p)
}

#' Heatmap of Median ECC
#'
#' Plot a heatmap of median element-centric consistency (ECC) scores for combinations of two parameters.
#' Expected input is a data frame with a 'params' column (a list of parameter settings) and a 'median_ecc' column.
#'
#' @param parameter_search_results_df A data frame containing a 'params' column (list of parameter settings)
#'        and a 'median_ecc' column.
#' @param key1 (Optional) A character string specifying the first parameter for the heatmap.
#'        If NULL, the first parameter in the data is used.
#' @param key2 (Optional) A character string specifying the second parameter for the heatmap.
#'        If NULL, the second parameter in the data is used.
#' @param agg_func (Optional) An aggregation function (e.g., mean, median). Default is median.
#' @param fixed_params (Optional) A named list of parameters to fix at specific values.
#'
#' @return A ggplot object representing the generated heatmap.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#'   param_grid <- list(centers = 2:4, nstart = c(1, 5))
#'   results <- parameter_searcher(
#'     data = data,
#'     clustering_algo = kmeans,
#'     labels_name = "cluster",
#'     param_grid = param_grid,
#'     n_runs = 30,
#'     verbose = TRUE
#'   )
#'   results_df <- results$results_df
#'   plot_heatmap(results_df, key1 = "centers", key2 = "nstart", agg_func = median)
#' }
#'
#' @export
plot_heatmap <- function(parameter_search_results_df,
                         key1 = NULL,
                         key2 = NULL,
                         agg_func = median,
                         fixed_params = NULL) {

  # Check if 'params' and 'median_ecc' columns exist
  if (!all(c("params", "median_ecc") %in% names(parameter_search_results_df))) {
    stop("The data frame must contain 'params' and 'median_ecc' columns.")
  }

  # Extract parameter names from the first element
  first_params <- parameter_search_results_df$params[[1]]
  param_keys <- names(first_params)

  # Select keys for the heatmap
  if (length(param_keys) >= 2) {
    if (is.null(key1)) key1 <- param_keys[1]
    if (is.null(key2)) key2 <- param_keys[2]
    message(sprintf("Selected keys for visualization: %s, %s", key1, key2))
  } else {
    stop("You must have at least two parameters to create a heatmap.")
  }

  # Verify that key1 and key2 are in param_keys
  if (!(key1 %in% param_keys) || !(key2 %in% param_keys)) {
    stop("Specify valid 'key1' and 'key2' parameters for visualization.")
  }

  message("Creating DataFrame from parameter search results...")

  # Determine the data type of each key from the first row
  key1_type <- class(first_params[[key1]])
  key2_type <- class(first_params[[key2]])

  get_map_func <- function(type) {
    if (type %in% c("numeric", "integer", "double")) {
      return(purrr::map_dbl)
    } else if (type %in% c("character", "factor", "logical")) {
      return(purrr::map_chr)
    } else {
      stop(sprintf("Unsupported data type: %s", type))
    }
  }

  map_func_key1 <- get_map_func(key1_type)
  map_func_key2 <- get_map_func(key2_type)

  # Instead of piping, we'll do step-by-step calls:
  params_df <- dplyr::mutate(
    parameter_search_results_df,
    !!rlang::sym(key1) := map_func_key1(.data$params, key1),
    !!rlang::sym(key2) := map_func_key2(.data$params, key2),
    median_ecc = as.numeric(.data$median_ecc)
  )
  params_df <- dplyr::select(params_df, -params)

  if (key1_type %in% c("character", "factor", "logical")) {
    params_df[[key1]] <- as.factor(params_df[[key1]])
  }
  if (key2_type %in% c("character", "factor", "logical")) {
    params_df[[key2]] <- as.factor(params_df[[key2]])
  }

  # Apply fixed parameters if provided
  if (!is.null(fixed_params)) {
    message(sprintf(
      "Applying filters for fixed parameters: %s",
      paste(names(fixed_params), collapse = ", ")
    ))
    for (param in names(fixed_params)) {
      if (param %in% colnames(params_df)) {
        before_filter <- nrow(params_df)
        params_df <- dplyr::filter(params_df, .data[[param]] == fixed_params[[param]])
        after_filter <- nrow(params_df)
        message(sprintf(
          "Filtered %s to %s. Rows: %d -> %d",
          param, fixed_params[[param]], before_filter, after_filter
        ))
      } else {
        warning(sprintf("Parameter '%s' not in dataset. Skipping.", param))
      }
    }
  }

  # Check for duplicates
  unique_combinations <- dplyr::distinct(params_df, dplyr::across(dplyr::all_of(c(key1, key2))))
  if (nrow(params_df) != nrow(unique_combinations)) {
    message("Duplicate combinations found. Aggregating values...")
    heatmap_data <- dplyr::group_by(params_df, dplyr::across(dplyr::all_of(c(key1, key2)))) %>%
      dplyr::summarize(median_ecc = agg_func(.data$median_ecc, na.rm = TRUE), .groups = 'drop')
    message(sprintf("Aggregation with '%s' done.", deparse(substitute(agg_func))))
  } else {
    message("No duplicates found. Proceeding without aggregation.")
    heatmap_data <- params_df
  }

  heatmap_plot <- ggplot2::ggplot(heatmap_data, ggplot2::aes_string(x = key2, y = key1, fill = "median_ecc")) +
    ggplot2::geom_tile(color = "white") +
    viridis::scale_fill_viridis(name = "ECC") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Heatmap of Median ECC",
      x = key2,
      y = key1
    ) +
    ggplot2::theme(ggplot2::element_text(angle = 45, hjust = 1))

  return(heatmap_plot)
}

#' Plot GA Fitness Evolution
#'
#' Plot the best fitness evolution over generations from a genetic algorithm.
#'
#' This function is intended to be used after running the genetic_algorithm_feature_selector,
#' where the history of best fitness per generation is returned.
#'
#' @param ga_fs_history A data frame containing the history of the feature selection GA,
#'        with columns "generation" and "best_fitness".
#'
#' @return A ggplot object representing the best fitness evolution plot.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#'
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
#'   plot_ga_fitness_evolution(ga_results$history)
#' }
#'
#' @export
plot_ga_fitness_evolution <- function(ga_fs_history) {

  p <- ggplot2::ggplot(ga_fs_history, ggplot2::aes(x = generation, y = best_fitness)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Generation", y = "Best Fitness") +
    ggplot2::ggtitle("Feature Selection GA - Best Fitness Evolution")

  return(p)
}


#' Plot K-Fold Similarity Scores
#'
#' Plot aggregated element-centric similarity (ECS) scores across k-folds.
#'
#' This function aggregates ECS scores (using a specified aggregation function) from multiple k-folds
#' and plots them along with an overall aggregated score.
#'
#' @param kfolds_results A named list containing k-fold results, where each element includes an "el_score_vector"
#'        of similarity scores.
#' @param agg_func A function to aggregate scores (default: median).
#'
#' @return A ggplot object representing the k-fold similarity scores.
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
#'   plot_kfold_scores(kfold_results, agg_func = median)
#' }
#'
#' @export
plot_kfold_scores <- function(kfolds_results, agg_func = median) {

  data_list <- list()

  for (fold in names(kfolds_results)) {
    score <- agg_func(kfolds_results[[fold]]$el_score_vector)
    data_list <- append(data_list, list(data.frame(Fold = fold, Score = score)))
  }
  df <- do.call(rbind, data_list)

  # Sort folds by the numeric portion after an underscore
  fold_sort_key <- function(fold_name) {
    as.numeric(strsplit(fold_name, "_")[[1]][2])
  }
  sorted_folds <- names(sort(sapply(unique(df$Fold), fold_sort_key)))
  df$Fold <- factor(df$Fold, levels = sorted_folds, ordered = TRUE)

  # Overall score
  overall_score <- agg_func(df$Score)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Fold, y = Score, group = 1)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(
      yintercept = overall_score,
      linetype = 'dashed',
      color = 'red',
      size = 1
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = 'Fold',
      y = 'Element Similarity Score',
      title = sprintf('%s Element Similarity Scores Across K-Folds',
                            tools::toTitleCase(deparse(substitute(agg_func))))
    ) +
    ggplot2::annotate(
      'text',
      x = length(levels(df$Fold)),
      y = min(overall_score - 0.1, 1.0),
      label = sprintf('Overall %s: %.3f',
                            deparse(substitute(agg_func)), overall_score),
      color = 'red',
      hjust = 1
    ) +
    ggplot2::ylim(0, 1.0)

  return(p)
}
