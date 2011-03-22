##
## runit-bugs.r - Unit tests for reported bugs to avoid regressions
##

## Cast n argument from numeric to integer.
test.001 <- function() {
  data <- rtruncnorm(10, 0, 1, 0, 1)
  checkEquals(length(data), 10L)
  checkException(rtruncnorm("a", 0, 1, 0, 1))
}

## If n is a vector, return length(n) elements
test.002 <- function() {
  data <- rtruncnorm(1:10, 0, 1, 0, 1)
  checkEquals(length(data), length(1:10))
}

check_rtruncnorm <- function(a, b, mean, sd, n=2e7) {
  x <- rtruncnorm(n=n, a, b, mean, sd)
  mx <- mean(x)
  vx <- var(x)
  ev <- etruncnorm(a, b, mean, sd)
  vv <- vtruncnorm(a, b, mean, sd)
  checkTrue(all(x >= a))
  checkTrue(all(x <= b))
  checkEqualsNumeric(ev, mx, tolerance=1e-3)
  checkEqualsNumeric(vv, vx, tolerance=1e-3)
}

## These two regression tests are not checked by default because the
## are stochastic in nature and might fail even though everything is
## OK. If you change anything in rtruncnorm.c, be sure to run these
## and investigate if they fail.
if (FALSE) {  
  ## Wrong sampling strategy:
  test.003 <- function() {
    set.seed(42)
    check_rtruncnorm(-0.75, 0.75, 0, 1)
  }
  
  ## Wrong sampling strategy:
  test.004 <- function() {
    set.seed(42)
    check_rtruncnorm(-1, 4, 0, 1) ## Case 2 (a)
  }

  test.005 <- function() {
    set.seed(42)
    check_rtruncnorm(0.5, 1.5, 6.5, 0.8) ## Case 3s (b)
  }
  
  ## Magic constants for reference:
  ##
  ##   t1 = 0.15, t2 = 2.18, t3 = 0.725, t4 = 0.45
  ## 
  test.006 <- function() {
    set.seed(42)
    check_rtruncnorm(-3, -2, 0, 1) # 3s 
  }

  test.007 <- function() {
    set.seed(42)
    check_rtruncnorm(-2, -1, 0, 1)
  }

  test.008 <- function() {
    set.seed(42)
    check_rtruncnorm(1, 2, 0, 1)
  }
}
