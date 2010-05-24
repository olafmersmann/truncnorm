##
## truncnorm.R - Interface to truncnorm.c
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

dtruncnorm <- function(x, a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_dtruncnorm", x, a, b, mean, sd)

ptruncnorm <- function(q, a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_ptruncnorm", q, a, b, mean, sd)

qtruncnorm <- function(p, a=-Inf, b=Inf, mean=0, sd=1) {  
  ## Taken out of package 'msm'
  ret <- numeric(length(p))
  ret[p == 1] <- Inf
  ret[p == 0] <- -Inf
  ret[p < 0 | p > 1] <- NaN
  ind <- (p > 0 & p < 1)
  if (any(ind)) {
    hind <- seq(along = p)[ind]
    h <- function(y) {
      (ptruncnorm(y, a=-Inf, b=Inf, mean, sd) - p)[hind[i]]
    }
    ptmp <- numeric(length(p[ind]))
    for (i in 1:length(p[ind])) {
      interval <- c(-1, 1)
      while (h(interval[1]) * h(interval[2]) >= 0) interval <- interval + 
        c(-1, 1) * 0.5 * (interval[2] - interval[1])
      ptmp[i] <- uniroot(h, interval, tol = .Machine$double.eps)$root
    }
    ret[ind] <- ptmp
  }
  return (ret)
}

rtruncnorm <- function(n, a=-Inf, b=Inf, mean=0, sd=1) {
  if (length(n) > 1)
    n <- length(n)
  if (length(n) > 1)
    n <- length(n)
  else if (!is.numeric(n))
    stop("non-numeric argument n.")
  .Call("do_rtruncnorm", as.integer(n), a, b, mean, sd)
}

etruncnorm <- function(a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_etruncnorm", a, b, mean, sd)

vtruncnorm <- function(a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_vtruncnorm", a, b, mean, sd)
