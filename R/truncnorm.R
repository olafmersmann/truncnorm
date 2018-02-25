##
## truncnorm.R - Interface to truncnorm.c
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

dtruncnorm <- function(x, a=-Inf, b=Inf, mean=0, sd=1)
  .Call(C_do_dtruncnorm, as.numeric(x),
        as.numeric(a), as.numeric(b), as.numeric(mean), as.numeric(sd))

ptruncnorm <- function(q, a=-Inf, b=Inf, mean=0, sd=1)
  .Call(C_do_ptruncnorm, as.numeric(q),
        as.numeric(a), as.numeric(b), as.numeric(mean), as.numeric(sd))

qtruncnorm <- function(p, a=-Inf, b=Inf, mean=0, sd=1)
  .Call(C_do_qtruncnorm, as.numeric(p),
        as.numeric(a), as.numeric(b), as.numeric(mean), as.numeric(sd))

rtruncnorm <- function(n, a=-Inf, b=Inf, mean=0, sd=1) {
  stopifnot(length(a) > 0,
            length(b) > 0,
            length(mean) > 0,
            length(sd) > 0)
  if (length(n) > 1)
    n <- length(n)
  else if (!is.numeric(n))
    stop("non-numeric argument n.")
  else if (n == 0)
    return(NULL)
  .Call(C_do_rtruncnorm, as.integer(n),
        as.numeric(a), as.numeric(b), as.numeric(mean), as.numeric(sd))
}

etruncnorm <- function(a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_etruncnorm",
        as.numeric(a), as.numeric(b), as.numeric(mean), as.numeric(sd))

vtruncnorm <- function(a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_vtruncnorm", 
        as.numeric(a), as.numeric(b), as.numeric(mean), as.numeric(sd))
