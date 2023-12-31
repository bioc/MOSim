% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClass.R
\docType{class}
\name{MOSimulation-class}
\alias{MOSimulation-class}
\title{This class manages the global simulation process, like associating genes with
gene classes, regulatory programs and other settings. Finally it will
initialize the simulators with their options that will use the previously
generated settings to simulate the data.}
\description{
This class manages the global simulation process, like associating genes with
gene classes, regulatory programs and other settings. Finally it will
initialize the simulators with their options that will use the previously
generated settings to simulate the data.
}
\section{Slots}{

\describe{
\item{\code{simulators}}{Vector containing either S4 initialized classes of
simulators or a list with the class name as keys, and its options as value,
see example.}

\item{\code{totalGenes}}{A number with the total number of genes including not
expressed. Overwritten if a genome reference is provided. Currently not used
as we force to provide real data.}

\item{\code{diffGenes}}{A number with the total number of differential genes (if
value > 1) or \% or total genes (if value < 1).}

\item{\code{numberReps}}{Number of replicates of the experiment.}

\item{\code{numberGroups}}{Number of samples considered on the experiment.}

\item{\code{times}}{Numeric vector containing the measured times. If numberGroups <
2, the number of times must be at least 2.}

\item{\code{geneNames}}{Read only. List containing the IDs of the genes. Overwrited
by the genome reference if provided. Currently not used as we force to
provide real data.}

\item{\code{simSettings}}{List of settings that overrides initializing the
configuration of the simulation by passing a previously generated list.
This could be used to tweak by hand the assigned profiles, genes,
regulatory programs, etc.}

\item{\code{noiseFunction}}{Noise function to apply when simulating counts. Must
accept the parameter 'n' and return a vector of the same length. Defaults
to `rnorm`}

\item{\code{profiles}}{Named list containing the patterns with their coefficients.}

\item{\code{profileProbs}}{Numeric vector with the probabilities to assign each of
the patterns. Defaults to 0.2 for each.}

\item{\code{noiseParams}}{Default noise parameters to be used with noise function.}

\item{\code{depth}}{Default depth to simulate.}

\item{\code{TFtoGene}}{Boolean (for default data) or 3 column data frame containing
Symbol-TFGene-LinkedGene}

\item{\code{minMaxQuantile}}{Numeric vector of length 2 indicating the quantiles to
use in order to retrieve the absolute minimum and maximum value that a
differentially expressed feature can have.}

\item{\code{minMaxFC}}{Numeric vector of length 2 indicating the minimum and maximum
fold-change that a differentially expressed feature can have.}
}}

\keyword{internal}
