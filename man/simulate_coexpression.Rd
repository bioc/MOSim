% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_coexpression.R
\name{simulate_coexpression}
\alias{simulate_coexpression}
\title{simulate coexpression}
\usage{
simulate_coexpression(
  sim_matrix,
  feature_no,
  cellTypes,
  patterns,
  cluster_size = NULL
)
}
\arguments{
\item{sim_matrix}{Matrix with rows as features and columns as cells}

\item{feature_no}{Total number of features to be distributed between the 
coexpression clusters}

\item{cellTypes}{list where the i-th element of the list contains the column 
indices for i-th experimental conditions. List must be a named list.}

\item{patterns}{Tibble with TRUE FALSE depicting the cluster patterns to
simulate. Generated by the user or by \code{make_cluster_patterns}.}

\item{cluster_size}{OPTIONAL. It may be inputted by the user. By default, 
its the number of features divided by the number of patterns to generate.}
}
\value{
the simulated coexpression
}
\description{
Adapted from ACORDE (https://github.com/ConesaLab/acorde) to adapt to our
data input type. Simulates coexpression of genes along celltypes
}
