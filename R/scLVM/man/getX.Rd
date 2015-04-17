\name{getX}
\alias{getX}
\alias{getK}
\alias{getKi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get K, Ki and X from gpCLVM object.
}
\description{
Retrieves latent factor, covariance matrix, and interaction covariance matrix from object of gpCLVM class.
}
\usage{
getX(gpclvm)
getK(gpclvm)
getKi(gpclvm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpclvm}{
Object of gpCLVM class.
}
}
\value{
The estimated latent factor X/covariance matrix K/interaction covariance matrix Ki.
}
