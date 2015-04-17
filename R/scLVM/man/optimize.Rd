\name{optimize}
\alias{optimize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimize conditional GPLVM
}
\description{
A GPLVM is fit to the data in the gpCLVM object, conditioning on some known factors.
}
\usage{
optimize(gpclvm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpclvm}{
Object of gpCLVM class.
}

}
\value{
Updated object of gpCLVM class with fitted covariance matrices and hidden factor.
}