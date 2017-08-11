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
dtruncnorm(x, a=-Inf, b=Inf, mean = 0, sd = 1)
ptruncnorm(q, a=-Inf, b=Inf, mean = 0, sd = 1)
qtruncnorm(p, a=-Inf, b=Inf, mean = 0, sd = 1)
rtruncnorm(n, a=-Inf, b=Inf, mean = 0, sd = 1)
etruncnorm(a=-Inf, b=Inf, mean=0, sd=1)
vtruncnorm(a=-Inf, b=Inf, mean=0, sd=1)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilites.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.}
  \item{a}{vector of lower bounds. These may be \code{-Inf}}
  \item{b}{vector of upper bounds. These may be \code{Inf}}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
}
\details{
  If \code{mean} or \code{sd} are not specified they assume the default values of
  \code{0} and \code{1}, respectively. The values of \code{a}, \code{b},
  \code{mean} and \code{sd} are recycled as needed.

  The numerical arguments other than \code{n} are recycled to the length
  of the result.
}
\value{
  'dtruncnorm' gives the density, 'ptruncnorm' gives the distribution
  function, 'qtruncnorm' gives the quantile function, 'rtruncnorm'
  generates random deviates, 'etruncnorm' gives the expected value and
  'vtruncnorm' the variance of the distribution.
}
\references{
  The accept-reject sampler follows the description given in

  Geweke, J. (1991). \emph{Efficient simulation from the multivariate normal
  and student-t distributions subject to linear constraints}. In
  Computing Science and Statistics: Proceedings of the 23rd Symposium on
  the Interface, Ed. E. Keramidas and S. Kaufman, pp. 571-8. Fairfax
  Station, VA: Interface Foundation of North America.
}
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de},
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de} and
  Björn Bornkamp \email{bornkamp@statistik.tu-dortmund.de} who donated a
  much improved \code{rtruncnorm} implementation using an accept-reject sampler.  
}
\keyword{distribution}
