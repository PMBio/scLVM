\name{getSymbols}
\alias{getSymbols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get gene symbols from ENSEMBL IDs.
}
\description{
Get gene symbols from ENSEMBL IDs.
}
\usage{
getSymbols(ensIds, species = "mMus")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ensIds}{
Vector of ENSEMBL IDs.
}
  \item{species}{
Species, either 'mMus' or 'Hs'.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character vector containig the gene symblos corresponding to the supplied ENSEMBL IDs.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ensIds, species = "mMus") 
{
    if (!(species \%in\% c("mMus", "Hs"))) {
        stop("'species' needs to be either 'mMus' or 'Hs'")
    }
    if (species == "mMus") {
        require(org.Mm.eg.db)
        x <- org.Mm.egSYMBOL
        xxenseg <- AnnotationDbi::as.list(org.Mm.egENSEMBL2EG)
    }
    else {
        require(org.Hs.eg.db)
        x <- org.Hs.egSYMBOL
        xxenseg <- AnnotationDbi::as.list(org.Hs.egENSEMBL2EG)
    }
    gene_names = ensIds
    mapped_genes <- mappedkeys(x)
    xx <- as.list(x[mapped_genes])
    gene_syms = unlist(xx[unlist(xxenseg[gene_names])])
    gene_names_list <- (lapply(xxenseg[gene_names], function(x) {
        if (is.null(x)) {
            x = NA
        }
        else {
            x = x[1]
        }
    }))
    sym_names = unlist(lapply(xx[unlist(gene_names_list)], function(x) {
        if (is.null(x)) {
            x = NA
        }
        else {
            x = x[1]
        }
    }))
    sym_names[is.na(sym_names)] = gene_names[is.na(sym_names)]
    sym_names
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
