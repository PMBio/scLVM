\name{fitLMM}
\alias{fitLMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit linear mixed model (lmm).
}
\description{
Fit lmm for calculating gene-gene correlation coefficients, controlling for hidden factors such as cell cycle which are modelled as random effcts. 
}
\usage{
fitLMM(sclvm, K=NULL,i0=i0,i1=i1,verbose=TRUE, geneID=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sclvm}{
sclvm object for which the lm is fitted. 
}
  \item{K}{
covariance matrix corresponding to hidden factor to be controlled for. If multiple hidden factors are present, the variance decompostion has to be performed before running fitLMM.
}
  \item{i0}{
Index of first gene for which the model is fit (correlations are calculated to all other genes).
}
  \item{i1}{
Index of last gene for which the model is fit.
}
  \item{verbose}{
Whether to display detailed status messages.
}
  \item{geneID}{
Gene identifiers. 
}
}

\value{
List with

\item{beta}{Matirx of correlation coefficients}
\item{pv}{Matrix of p-values}
}

