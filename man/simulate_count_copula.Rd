% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scDesign2_simulate_revised.R
\name{simulate_count_copula}
\alias{simulate_count_copula}
\title{Simulate a count matrix for a single cell type based on a copula model}
\usage{
simulate_count_copula(copula_result, n = 100, marginal = c("nb", "Gamma"))
}
\arguments{
\item{copula_result}{A list that contains the parameters of a copula model.}

\item{n}{An integer value that indicates the number of cells to generate.}

\item{marginal}{A character string that indicates whether the generated values should
stay as discrete or switch to continuous. Default value is 'nb', which
should be used for generating a count marix. The alternative 'Gamma' is
only needed when this function is being called by other functions that
generate data with a user-specified sequencing depth. Normally, users
do not need to change this value.}
}
\value{
A matrix of shape p by n that contains the simulated count values. p is derived from
\code{copula_result}
}
\description{
Simulate a count matrix for a single cell type based on a copula model
}
