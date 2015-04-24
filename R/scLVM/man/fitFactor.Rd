\name{fitGPLVM}
\alias{fitGPLVM}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits a latent factor using a GPLVM to obtain a low-rank covariance structure. 
}
\description{
Fits a latent factor using a Gaussian Process Latent Variable Model, projecting the gene expression data annotated to a specific process into a low dimensional space. This results in a cell-cell covariance matrix reflecting similarities between cells in terms of the proces of interest.
}
\usage{
fitFactor(sclvm,idx=NULL,geneSet = NULL, XKnown = NULL, k=1,standardize=FALSE, use_ard=FALSE,interaction=TRUE, initMethod='fast')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sclvm}{
the sclvm object for which to fit the GPLVM
}
  \item{idx}{
indices of genes annotated to the process of interest (e.g. cell cycle). The GPLVM will be fit on these genes.
}
  \item{geneSet}{
character vector of gene IDs of genes annotated to the process of interest (e.g. cell cycle). The GPLVM will be fit on these genes.
}
  \item{XKnown}{
numeric vector or matrix containing known factor(s) on which the method should condition on when fitting the latent factor.
}
\item{k}{
the number of latent factors to be fit (has to be less than the number of cells)
}
  \item{standardize}{
 whether the gene expression data should be standardized before fitting the latent factor (ususally FALSE).
}
  \item{use_ard}{
whether automatic relevance determination (ard) should be used when fitign the GPLVM. ard automatically switches off factors which are not needed to explain the data and can be used to detrmine k. If TRUE, a large k should be chosen.
}

  \item{interaction}{
whether an interaction term should be fitted beween knwon facotr and latent factor (default TRUE when known factor is passed.)
}
  \item{initMethod}{
how to initialize the latent factor. Currently 'fast' and 'random' are supported. For method 'fast', the known factor is regressed out of gene set passed to the method and a PCA is perfromed on these data; the PCs are then used to iniitalize the latent factors. Method 'random' perfroms random initialization.

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A list containing
  \item{K }{The estimated covariance matirx}
  \item{X }{The estimated underlying latent factor}
  \item{X_ard}{The relevance of the indivisual factors (for use_ard=TRUE only)}
%% ...
}
\references{
Lawrence, N. D. (2004). Gaussian process latent variable models for visualisation of high dimensional data. \emph{Advances in neural information processing systems}, \bold{16}(3), 329--336.

Buettner, F. and Theis, F. J. (2012). A novel approach for resolving differences in single-cell gene expression patterns from zygote to blastocyst. \emph{Bioinformatics}, 28(18), i626-i632.
}