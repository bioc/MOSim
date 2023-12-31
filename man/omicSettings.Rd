% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MOSim.R
\name{omicSettings}
\alias{omicSettings}
\title{Retrieves the settings used in a simulation}
\usage{
omicSettings(
  simulation,
  omics = NULL,
  association = FALSE,
  reverse = FALSE,
  only.linked = FALSE,
  prefix = FALSE,
  include.lagged = TRUE
)
}
\arguments{
\item{simulation}{A MOSimulation object.}

\item{omics}{List of omics to retrieve the settings.}

\item{association}{A boolean indicating if the association must also be
returned for the regulators.}

\item{reverse}{A boolean, swap the column order in the association list in
case we want to use the output directly and the program requires a
different ordering.}

\item{only.linked}{Return only the interactions that have an effect.}

\item{prefix}{Logical indicating if the name of the omic should prefix the
name of the regulator.}

\item{include.lagged}{Logical indicating if interactions with transitory
profile and different minimum/maximum time point between gene and regulator
should be included or not.}
}
\value{
A list containing a data frame with the settings used to simulate
  each of the indicated omics. If association is TRUE, it will be a list with
  3 keys: 'associations', 'settings' and 'regulators', with the first two
  keys being a list containing the information for the selected omics and the
  last one a global data frame giving the merged information.
}
\description{
Retrieves the settings used in a simulation
}
\examples{

omic_list <- c("RNA-seq", "miRNA-seq")
multi_simulation <- mosim(omics = omic_list)

# This will be a data frame with RNA-seq settings (DE flag, profiles)
rnaseq_settings <- omicSettings(multi_simulation, "RNA-seq")

# This will be a list containing all the simulated omics (RNA-seq
# and DNase-seq in this case)
all_settings <- omicSettings(multi_simulation)

}
