\name{truncnorm}
\alias{dtruncnorm}
\alias{ptruncnorm}
\alias{qtruncnorm}
\alias{rtruncnorm}
\alias{etruncnorm}
\alias{vtruncnorm}
\title{The Truncated Normal Distribution}
\description{
  Density, distribution function, quantile function, random generation
  and expected value function for the truncated normal distribution
  with mean equal to 'mean' and standard deviation equal to 'sd'.
}
\usage{
dtruncnorm(x, a, b, mean = 0, sd = 1)
ptruncnorm(q, a, b, mean = 0, sd = 1)
qtruncnorm(p, a, b, mean = 0, sd = 1)
rtruncnorm(n, a, b, mean = 0, sd = 1)
etruncnorm(a, b, mean, sd)
vtruncnorm(a, b, mean, sd)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilites.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.}
  \item{a}{vector of lower bounds.}
  \item{b}{vector of upper bounds.}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
}
\details{
  If \code{mean} or \code{sd} are not specified they assume the default values of
  \code{0} and \code{1}, respectively. The values of \code{a}, \code{b},
  \code{mean} and \code{sd} are recycled as needed.
}
\value{
  'dtruncnorm' gives the density, 'ptruncnorm' gives the distribution
  function, 'qtruncnorm' gives the quantile function, 'rtruncnorm'
  generates random deviates, 'etruncnorm' gives the expected value and
  'vtruncnorm' the variance of the distirbution.
}
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}.
  'qtruncnorm' was originally published in package 'msm' written by
  C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}
}
\keyword{distribution}
