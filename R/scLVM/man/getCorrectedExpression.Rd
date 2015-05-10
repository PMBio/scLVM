\name{getCorrectedExpression}
\alias{getCorrectedExpression}

\alias{setCorrectedExpression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get corrected gene expression data
}
\description{
Gets corrected gene expression data where confounding hidden factors have been regressed out.
}
\usage{
getCorrectedExpression(sclvm, rand_eff_ids=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sclvm}{
sclvm object. varianceDecomposition has to be called before using this method.
}
  \item{rand_eff_ids}{
vector of integers indicating which hidden factors should be regressed out if multiple factors are present. If NULL (default), all hidden factors are regressed out.
}
}
\value{
Matrix of corrected gene expression data.
}
