\name{getVarianceComponents}
\alias{getVarianceComponents}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get variance components
}
\description{
Gets variance components from sclvm object. 
}
\usage{
getVarianceComponents(sclvm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sclvm}{
scLVM object. varianceDecomposition has to be run before calling this function.
}

}
\value{
 A list containing
 \item{var}{a matrix containing for each gene the relative proportion of the variance explained by all factors}
 \item{conv}{a vector indicating whether the optimization for the variance decomposition has converged.}
}
