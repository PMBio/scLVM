\name{varianceDecomposition}
\alias{varianceDecomposition}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance decomposition.
}
\description{
Variance is decomposed into technical noise, hidden biological (and technical) factors and residual biological noise. 
}
\usage{
varianceDecomposition(sclvm, K=NULL,i0=1,i1=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sclvm}{
sclvm object for which the variance decomposition is perfromed.
}
  \item{K}{
covariance matrix corresponding to hidden factor (usually estimated with fitGPLVM or gpCLVM). Can be a list if several hidden processes are present.
}
  \item{i0}{
Index of first gene for which the model is fit (correlations are calculated to all other genes).
}
  \item{i1}{
Index of last gene for which the model is fit.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An updated sclvm object.
}