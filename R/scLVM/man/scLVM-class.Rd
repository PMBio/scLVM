\name{scLVM-class}
\Rdversion{1.1}
\docType{class}
\alias{scLVM-class}
\alias{scLVM}


\title{Class \code{"scLVM"}}
\description{
The main class used by scLVM to hold processed single-cell RNA-Seq data along with technical noise estimates and results from the variance decomposition.
}
\section{Details}{
This class is initialized from a matirx of log-transformed read counts, estimates of baseline variability (technical noise) and gene identifiers.

}
\section{Slots}{
  \describe{
    \item{\code{Y}:}{Object of class \code{"matrix"}, containing the normalized gene expression data }
    \item{\code{geneID}:}{Object of class \code{"character"}, contianing the gene idenitifiers }
    \item{\code{tech_noise}:}{Object of class \code{"numeric"}, contianing estimates of baseline variability (technical noise)  }
    \item{\code{Ycorr}:}{Object of class \code{"matrix"},  contianing the (cell-cycle-) corrected gene expression data }
    \item{\code{var}:}{Object of class \code{"matrix"}, containing results from the variance decompisition (fraction of variance explained by all factors)}
    \item{\code{conv}:}{Object of class \code{"logical"}, containing information on whether the optimization routine as part of the variance decomposition has converged.}
  }
}

