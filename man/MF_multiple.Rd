% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mainfun.R
\name{MF_multiple}
\alias{MF_multiple}
\title{\code{MF_multiple}:\cr
Multi-functionality measures in multiple ecosystems. This is a function in the framework of Hill-Chao numbers}
\usage{
MF_multiple(func_data, species_data = NULL, q = c(0, 1, 2), by_group = NULL)
}
\arguments{
\item{func_data}{the ecosystem function data can be input as a vector of functions (for a single assemblage), matrix/data.frame (assemblages by functions).\cr
The data input must be normalized between 0 and 1 already.\cr
For \code{by_group = NULL} the func_data must contain only the ecosystem function columns. Otherwise, you must add the \code{by_group} column in data.
For \code{species_data} is not \code{NULL}, the rownames of func_data should be names of plotID.}

\item{species_data}{the species abundance data must include three columns: 'plotID', 'species' and 'abundance'. Default is \code{NULL}.}

\item{q}{a numerical vector specifying the diversity orders. Default is 0, 1 and 2.}

\item{by_group}{name of the column to be paired by group. Default is \code{NULL}.}
}
\value{
a data.frame with columns 'Order.q' and 'MF_gamma', 'MF_alpha' and 'MF_beta' for uncorrelated output, 'MF_gamma_cor', 'MF_alpha_cor' and 'MF_beta_cor' for correlated output.\cr
For \code{species_data} is not \code{NULL}, the data.frame will contain 'Species_gamma', 'Species_alpha' and 'Species_beta'.
}
\description{
\code{MF_multiple}:\cr
Multi-functionality measures in multiple ecosystems. This is a function in the framework of Hill-Chao numbers
}
