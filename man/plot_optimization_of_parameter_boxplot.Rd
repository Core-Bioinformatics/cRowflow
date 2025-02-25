% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_optimization_of_parameter_boxplot}
\alias{plot_optimization_of_parameter_boxplot}
\title{Plot Optimization of Parameter (Boxplot)}
\usage{
plot_optimization_of_parameter_boxplot(
  parameter_optimizer_results_df,
  name_parameter = "param",
  title = "Parameter Optimization - Boxplot Distribution"
)
}
\arguments{
\item{parameter_optimizer_results_df}{A data frame with columns "param" and "ecc".}

\item{name_parameter}{A character string for the x-axis label (default: "param").}

\item{title}{A character string for the plot title (default: "Parameter Optimization - Boxplot Distribution").}
}
\value{
A ggplot object representing the boxplot visualization.
}
\description{
Generate a boxplot visualization of element-centric consistency (ECC) values for different parameter settings.
}
\details{
Expected input is a data frame with at least two columns:
\describe{
  \item{\code{param}}{A string representing the parameter and its value (e.g., "n_clusters_2") tested during optimization.}
  \item{\code{ecc}}{A list-column where each entry is a numeric vector of ECC values corresponding to that parameter setting.}
}
}
\examples{
\dontrun{
  set.seed(42)
  # Generate synthetic data: 200 samples with 10 features
  data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
  results <- parameter_optimizer(
    data = data,
    clustering_algo = kmeans,
    labels_name = "cluster",
    parameters_optimise_list = list(centers = 2:4),
    n_runs = 30,
    verbose = TRUE
  )
  results_df <- results$results_df
  plot_optimization_of_parameter_boxplot(results_df)
}
}
