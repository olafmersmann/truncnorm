/*
 * rtruncnorm.c - Random truncated normal number generator.
 *
 * Authors:
 *  Björn Bornkamp   <bornkamp@statistik.tu-dortmund.de>
 *  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
 */

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#include "sexp_macros.h"

#define ALLOC_REAL_VECTOR(S, D, N)		       \
    SEXP S;					       \
    PROTECT(S = allocVector(REALSXP, N));	       \
    double *D = REAL(S);

#ifndef MAX
#define MAX(A, B) ((A>B)?(A):(B))
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>

static const double t1 = 0.15;
static const double t2 = 2.18;
static const double t3 = 0.725;
static const double t4 = 0.45;

/* Exponential rejection sampling (a,inf) */
static R_INLINE double ers_a_inf(double a) {
    const double ainv = 1.0 / a;
    double x,rho;
    do{
	x = rexp(ainv) + a; /* rexp works with 1/lambda */
	rho = exp(-0.5 * pow((x - a), 2));
    } while (runif(0, 1) > rho);
    return x;
}

/* Exponential rejection sampling (a,b) */
static R_INLINE double ers_a_b(double a, double b) {
    const double ainv = 1.0 / a;
    double x, rho;
    do{
	x = rexp(ainv) + a; /* rexp works with 1/lambda */
	rho = exp(-0.5 * pow((x-a), 2));
    } while (runif(0, 1) > rho || x > b);
    return x;
}

/* Normal rejection sampling (a,b) */
static R_INLINE double nrs_a_b(double a, double b){
  double x = a - 1.0;
  while(x < a || x > b){
    x = rnorm(0, 1);
  }
  return x;
}

/* Normal rejection sampling (a,inf) */
static R_INLINE double nrs_a_inf(double a){
  double x = a - 1.0;
  while(x < a){
    x = rnorm(0, 1);
  }
  return x;
}

/* Half-normal rejection sampling */
double hnrs_a_b(double a, double b){
  double x = a - 1.0;
  while(x < a || x > b) {
    x = rnorm(0, 1);
    x = fabs(x);
  }
  return x;
}

/* Uniform rejection sampling */
static R_INLINE double urs_a_b(double a, double b){
    const double phi_a = dnorm(a, 0.0, 1.0, FALSE);
    double x = 0.0, u = 0.0;
    if (a < 0 && b > 0) {
	do{
	    x = runif(a, b);
	} while (runif(0, 1) * M_1_SQRT_2PI > dnorm(x,0,1,0));
    } else {
	do{
	    x = runif(a, b);
	    u = runif(0, 1);
	} while (runif(0, 1) * phi_a > dnorm(x, 0.0, 1.0, FALSE));
    }
    return x;
}

/* Previously this was refered to as type 1 sampling: */
static inline double r_lefttruncnorm(double a, double mean, double sd) {
    const double alpha = (a - mean) / sd;
    if (alpha < t4) {
	return mean + sd * nrs_a_inf(alpha);
    } else {
	return mean + sd * ers_a_inf(alpha);
    }
}

static R_INLINE double r_righttruncnorm(double b, double mean, double sd) {
    const double beta = (b - mean) / sd;
    /* Exploit symmetry: */
    return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0);
}

static R_INLINE double r_truncnorm(double a, double b, double mean, double sd) {
    const double alpha = (a - mean) / sd;
    const double beta = (b - mean) / sd;
    const double phi_a = dnorm(a, 0.0, 1.0, FALSE);
    const double phi_b = dnorm(b, 0.0, 1.0, FALSE);
    if (beta <= alpha) {
	return NA_REAL;
    } else if (alpha < 0 && 0 < beta) {
	if (phi_a < t1 && phi_b < t1) {
	    return mean + sd * nrs_a_b(alpha, beta);
	} else {
	    return mean + sd * urs_a_b(alpha, beta);
	}
    } else if (alpha >= 0) {
	if (phi_a / phi_b <= t2) {
	    return mean + sd * urs_a_b(alpha, beta);
	} else {
	    if (alpha < t3) {
		return mean + sd * hnrs_a_b(alpha, beta);
	    } else {
		return mean + sd * ers_a_b(alpha, beta);
	    }
	}
    } else { /* beta <= 0 */
	if (phi_a / phi_b <= t2) {
	    return mean - sd * urs_a_b(-beta, -alpha);
	} else {
	    if (alpha < t3) {
		return mean - sd * hnrs_a_b(-beta, -alpha);
	    } else {
		return mean - sd * ers_a_b(-beta, -alpha);
	    }
	}
    }
}

SEXP do_rtruncnorm(SEXP s_n, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, nn;
  UNPACK_INT(s_n, n);
  if (NA_INTEGER == n)
      error("n is NA - aborting.");
  UNPACK_REAL_VECTOR(s_a   , a   , n_a);
  UNPACK_REAL_VECTOR(s_b   , b   , n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);

  nn = MAX(n, MAX(MAX(n_a, n_b), MAX(n_mean, n_sd)));
  ALLOC_REAL_VECTOR(s_ret, ret, nn);

  GetRNGstate();
  for (i = 0; i < nn; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];

    if (R_FINITE(ca) && R_FINITE(cb)) {
	ret[i] = r_truncnorm(ca, cb, cmean, csd);
    } else if (R_NegInf == ca && R_FINITE(cb)) {
	ret[i] = r_righttruncnorm(cb, cmean, csd);
    } else if (R_FINITE(ca) && R_PosInf == cb) {
	ret[i] = r_lefttruncnorm(ca, cmean, csd);
    } else if (R_NegInf == ca && R_PosInf == cb) {
	ret[i] = rnorm(cmean, csd);
    } else {
	ret[i] = NA_REAL;
    }
  }
  PutRNGstate();
  UNPROTECT(1); /* s_ret */
  return s_ret;
}


