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
LMM(sclvm, expr = NULL, K=NULL, idx=NULL, geneSet=NULL, verbose=TRUE, recalc=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sclvm}{
sclvm object for which the lm is fitted. 
}
  \item{K}{
covariance matrix corresponding to hidden factor to be controlled for. If multiple hidden factors are present, the variance decompostion has to be performed before running fitLMM.
}
 \item{idx}{
indices of genes annotated to the process of interest (e.g. cell cycle). The GPLVM will be fit on these genes.
}
  \item{geneSet}{
character vector of gene IDs of genes annotated to the process of interest (e.g. cell cycle). The GPLVM will be fit on these genes.
}
\item{verbose}{
Whether to display detailed status messages.
}
\item{recalc}{
Whether to re-do the variance decomposition. recalc should only be set to FALSE by advanced users: scLVM then assumes that the random effects are the same as those for which the variance decompostion was performed earlier.
}

}

\value{
List with

\item{beta}{Matirx of correlation coefficients}
\item{pv}{Matrix of p-values}
}

