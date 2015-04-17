\name{gpCLVM-class}
\Rdversion{1.1}
\docType{class}
\alias{gpCLVM-class}
\alias{gpCLVM}


\title{Class \code{"gpCLVM"}}
\description{
The class used to fit conditional GPLVM models.
}
\section{Details}{
This class is initialized from a matrix of log-transformed read counts and a matirx of factors on which should be conditioned on. 
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Y}:}{Object of class \code{"matrix"}, containing the normalized and (log-)transformed read counts on which a GPLVM is fitted. These should be only the genes annotated to the hidden factor which is modelled.}
    \item{\code{X}:}{Object of class \code{"matrix"}, containing the factors on which the model conditions on when fitting the GPLVM.}
    \item{\code{K}:}{Object of class \code{"matrix"}, containing the estimated covariance matrix form the GPLVM fit. }
    \item{\code{Ki}:}{Object of class \code{"matrix"}, containing the interaction covariance matirx between X and the hidden factor uderlying Y. }
    \item{\code{interaction}:}{Object of class \code{"logical"}, indicating whether an interaction term should be fitted. }
  }
}

\references{
Lawrence, N. D. (2004). Gaussian process latent variable models for visualisation of high dimensional data. \emph{Advances in neural information processing systems}, \bold{16}(3), 329--336.

}

