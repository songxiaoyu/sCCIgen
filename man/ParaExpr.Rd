% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParameterDigest.R
\name{ParaExpr}
\alias{ParaExpr}
\title{Simualte gene expression data based on parameters}
\usage{
ParaExpr(
  para,
  cell_loc_list,
  expr,
  region,
  seed_list,
  model_params = NULL,
  ncores = 1
)
}
\arguments{
\item{para}{Parameters loaded and cleaned from the parameter file using function
\code{ParaDigest}.}

\item{cell_loc_list}{Simulated cell location data}

\item{expr}{Input or pre-cleaned expression data}

\item{region}{Input region of cells}

\item{seed_list}{Seeds for all simulated data}

\item{model_params}{The fitted models of genes.}

\item{ncores}{No. of cores for simulation}
}
\value{
Simulated gene expression data for each cell.
}
\description{
Simualte gene expression data based on parameters
}
