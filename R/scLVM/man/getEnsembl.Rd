\name{getEnsembl}
\alias{getEnsembl}
\title{
Get ENSEMBL IDs
}
\description{
Retrieve ENSEMBL Ids anootated to a GO term.
}
\usage{
getEnsembl(term, species = "mMus")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{term}{
GO term, e.g. 'GO:0007049' for cell cycle.
}
  \item{species}{
Species, either 'mMus' or 'Hs'.
}
}

\value{
Vector with all ENSEMBL Ids annotated to the specified GO term.
}

