\name{fit_cm}
\alias{fit_cm}
\title{Fit a conditional-moment function}
\usage{fit_cm(response, diagonal = NULL, target, instrument,
  method = c("sieve_md", "landweber", "pmmr"), control = list())}
\arguments{
\item{response}{Numeric response vector.}
\item{diagonal}{Numeric multiplicative/diagonal vector. Defaults to one.}
\item{target}{Variables indexing the unknown function.}
\item{instrument}{Variables defining the conditional moment.}
\item{method}{Estimation method.}
\item{control}{Method-specific control list.}
}
\description{Fits h in E[Y-D h(X) given Z]=0 using a common interface.}
\value{An object of class \code{cmbridge_fit}.}
