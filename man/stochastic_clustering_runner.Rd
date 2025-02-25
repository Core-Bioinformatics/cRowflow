% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stochastic_clustering_runner.R
\name{stochastic_clustering_runner}
\alias{stochastic_clustering_runner}
\title{Stochastic Clustering Runner}
\usage{
stochastic_clustering_runner(
  data,
  clustering_algo,
  labels_name = NULL,
  n_runs = 30,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{data}{A data matrix or data frame on which clustering is performed.}

\item{clustering_algo}{A stochastic clustering function that returns a vector of cluster labels.
The function should produce different results when the seed is changed (via \code{set.seed}).}

\item{labels_name}{(Optional) A string indicating the name of the element in the output of
\code{clustering_algo} that contains the cluster labels. If \code{NULL}, the function is assumed
to return the labels directly.}

\item{n_runs}{Number of times to repeat the clustering process. Defaults to 30.}

\item{verbose}{Logical; if \code{TRUE}, prints progress messages showing the seed used in each run.}

\item{...}{Additional arguments to be passed to \code{clustering_algo}.}
}
\value{
A list with components:
\describe{
  \item{\code{partitions}}{A list of integer vectors representing clustering assignments for each run.}
  \item{\code{partition_frequencies}}{A list mapping each unique partition (as a concatenated string)
        to its frequency of occurrence.}
  \item{\code{majority_voting_labels}}{An integer vector of consensus cluster labels derived from majority voting.}
  \item{\code{ecc}}{A numeric vector of element-centric consistency (ECC) scores indicating clustering stability.}
  \item{\code{seeds}}{An integer vector of the random seeds used in each run.}
}
}
\description{
Repeats stochastic clustering multiple times to assess stability.
}
\details{
This function runs a clustering algorithm multiple times with different random seeds
and evaluates the stability of results using Element-Centric Consistency (ECC).
It identifies in an element-wise precision the stability of clustering results and
provides a consensus labeling through majority voting.


The function generates seeds starting from 100 and increasing by 100 for each run. For each seed, it sets
the random number generator using \code{set.seed} to ensure reproducibility of the stochastic clustering.
The clustering results from each run are reconciled through majority voting and evaluated using ECC.
}
\examples{
\dontrun{
  set.seed(42)
  data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
  results <- stochastic_clustering_runner(
    data = data,
    clustering_algo = kmeans,
    labels_name = "cluster",
    n_runs = 30,
    verbose = TRUE,
    centers = 3
  )
  print("Majority Voting Labels:")
  print(results$majority_voting_labels)
  print("ECC:")
  print(results$ecc)
}

}
