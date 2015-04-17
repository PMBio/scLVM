\name{fitTechnicalNoise}
\alias{fitTechnicalNoise}

\title{
Estimates the technical noise/baseline variability.
}
\description{
Based on ERCC spike-ins (or endogeneous read counts if not available), the technical noise is estimated for normalized read counts.
}
\usage{
fitTechnicalNoise(nCountsEndo, nCountsERCC = NULL, use_ERCC = TRUE, fit_type = "counts", plot = TRUE, fit_opts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nCountsEndo}{
the matrix of normalised read counts mapped to endogenous genes.
}
  \item{nCountsERCC}{
the matirx of nromalised read counts mapped to ERCC spike-ins (NULL is not available).
}
  \item{use_ERCC}{
whether ERCC spike-ins should be used to estimate the technical noise.
}
  \item{fit_type}{
how the technical noise/baseline variability should be estimated, has to be one of 'counts', 'log' or 'logvar'. 'counts' will use the Brenencke et al. parameterization of the mean-squared ocefficient of variation relation (cv2) in the count space. 'log' uses a log-linear fit of mean-cv2 in the log space and 'logvar' uses Local Polynomial Regression Fitting (degree 2) to fit the mean-variance relation in the log-space.
}
  \item{plot}{
whether to plot the fit.
}
  \item{fit_opts}{
named list with options for the fit parameters, 'mincv2', 'quan', 'span' and 'minmean'. 'mincv2' and  'quan' specify a minimum mean value to exclude genes with low mean and hence high CV for fit_type 'counts'. 'span' and 'minmean' specify the minimum mean value and span for the loess fit (fit_type 'logvar'). If not provided default values are used (0.3 and 0.8 for 'mincv2' and  'quan'; 0.8 and 0.5 for 'span' and 'minmean').
}
}
\details{
When ERCC spike-ins are available, we recommend following Brennecke et al. When no spik-ins are available, we recommend using the log-linear fit based on all genes. 
}
\value{
a list with
\item{fit}{parameters of the fit}
\item{techNoiseLog}{estimates for technical noise/baseline variability for all genes in log-space}
}
\references{
Brennecke, P., Anders, S., Kim, J. K., Kołodziejczyk, A. A., Zhang, X., Proserpio, V., Baying, B., Benes, V., Teichmann, S. A., Marioni, J.C. and Heisler, M. G. (2013). Accounting for technical noise in single-cell RNA-seq experiments. \emph{Nature methods}, \bold{10}(11), 1093--1095.
}

