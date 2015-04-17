\name{getVariableGenes}
\alias{getVariableGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calls variable genes based on technical noise fit.
}
\description{
The fit estimates characterizing the baseline variability (typically obtained using the fitTechnicalNoise function), are used to call the most variable genes. 
}
\usage{
getVariableGenes(nCountsEndo, fit, method = "fit", threshold = 0.1, fit_type = NULL, sfEndo = NULL, sfERCC = NULL, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nCountsEndo}{
contains the normalized gene expression data
}
  \item{fit}{
the fit parameterizing the technical noise/baseline variability. Should be generated with fitTechnicalNoise.
}
  \item{method}{
method to call variable genes. Either 'fdr' for the statistical test following Brennecke et al. or 'fit' if only based on the fitted values. 
}
  \item{threshold}{
FDR cut-off. Only used if method is 'fdr'.
}
  \item{fit_type}{
how the technical noise/baseline variability was estimated, has to be one of 'counts', 'log' or 'logvar'. If not provided, it is guessed from the structure of the fit. 
}
  \item{sfEndo}{
size factor endogenous genes (only needed if method is 'fdr').
}
  \item{sfERCC}{
size factor ERCCs (only needed if method is 'fdr').
}
  \item{plot}{
whether the results should be plotted. 
}
}

\value{
A vector inidicaitng whether a gene is variable above technical noise/baseline.
}
\references{
Brennecke, P., Anders, S., Kim, J. K., Kołodziejczyk, A. A., Zhang, X., Proserpio, V., Baying, B., Benes, V., Teichmann, S. A., Marioni, J.C. and Heisler, M. G. (2013). Accounting for technical noise in single-cell RNA-seq experiments. \emph{Nature methods}, \bold{10}(11), 1093--1095.
}

