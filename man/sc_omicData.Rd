% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_MOSim.R
\name{sc_omicData}
\alias{sc_omicData}
\title{sc_omicData}
\usage{
sc_omicData(omics_types, data = NULL)
}
\arguments{
\item{omics_types}{A list of strings which can be either "scRNA-seq" or "scATAC-seq"}

\item{data}{A user input matrix with genes (peaks in case of scATAC-seq) as 
rows and cells as columns. By default, it loads the example data.
If a user input matrix is included, cell columns must be sorted by cell t
ype.}
}
\value{
a named list with omics type as name and the count matrix as value
}
\description{
Checks if the user defined data is in the correct format, or loads the default
multiomics pbmc dataset, a subset from SeuratData package
}
\examples{
# Simulate from PBMC
omicsList <- sc_omicData(list("scRNA-seq", "scATAC-seq"))

}
