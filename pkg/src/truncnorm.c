/*
 * truncnorm.c - Implementation of truncated normal distribution
 *
 * Authors:
 *  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
 *  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
 *  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
 */

#include <R.h>
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

/*
 * These routines calculate the expected value and variance of the
 * left, right and doubly truncated normal distribution. The only
 * tricky bit is the calculation of the variance of the doubly
 * truncated normal distribution. We use a decompostion of the
 * variance of a mixture of distributions to here for numerical
 * reasons. For details see:
 * 
 *   Foulley JL. A completion simulator for the two-sided truncated
 *   normal distribution. Genetics, selection, evolution 2000
 *   Nov-Dec;32(6): p. 631-635.
 */
static R_INLINE double e_lefttruncnorm(double a, double mean, double sd) {
    const double alpha = (a - mean) / sd;
    const double phi_a = dnorm(alpha, 0.0, 1.0, FALSE);
    const double Phi_a = pnorm(alpha, 0.0, 1.0, TRUE, FALSE);
    double res = mean + sd * (phi_a / (1.0 - Phi_a));
    return res;	
}

static R_INLINE double e_truncnorm(double a, double b, double mean, double sd) {
    const double alpha = (a - mean) / sd;
    const double beta = (b - mean) / sd;

    const double phi_a = dnorm(alpha, 0.0, 1.0, FALSE);
    const double Phi_a = pnorm(alpha, 0.0, 1.0, TRUE, FALSE);
    const double phi_b = dnorm(beta, 0.0, 1.0, FALSE);
    const double Phi_b = pnorm(beta, 0.0, 1.0, TRUE, FALSE);
    return mean + sd * (phi_b - phi_a) / ( Phi_a - Phi_b);
}

static R_INLINE double e_righttruncnorm(double b, double mean, double sd) {
    const double beta = (b - mean) / sd;
    const double phi_b = dnorm(beta, 0.0, 1.0, FALSE);
    const double Phi_b = pnorm(beta, 0.0, 1.0, TRUE, FALSE);
    return mean + sd * (-phi_b / Phi_b);
}

static R_INLINE double v_lefttruncnorm(double a, double mean, double sd) {
    const double alpha = (a - mean) / sd;
    const double phi_a = dnorm(alpha, 0.0, 1.0, FALSE);
    const double Phi_a = pnorm(alpha, 0.0, 1.0, TRUE, FALSE);
    const double lambda = phi_a / (1.0 - Phi_a);
    return (sd*sd*(1.0 - lambda * (lambda - alpha)));
}

static R_INLINE double v_righttruncnorm(double b, double mean, double sd) {
    return (v_lefttruncnorm(-b, -mean, sd));
}

static R_INLINE double v_truncnorm(double a, double b, double mean, double sd) {
    const double pi1 = pnorm(a, mean, sd, TRUE, FALSE);
    const double pi2 = pnorm(b, mean, sd, TRUE, FALSE) - pnorm(a, mean, sd, TRUE, FALSE);
    const double pi3 = pnorm(b, mean, sd, FALSE, FALSE); /* 1 - F(b) */

    const double e1 = e_righttruncnorm(a, mean, sd);
    const double e2 = e_truncnorm(a, b, mean, sd);
    const double e3 = e_lefttruncnorm(b, mean, sd);

    const double v  = sd * sd;
    const double v1 = v_righttruncnorm(a, mean, sd);
    const double v3 = v_lefttruncnorm(b, mean, sd);    

    const double c1 = pi1 * (v1 + (e1 - mean)*(e1 - mean));
    const double c3 = pi3 * (v3 + (e3 - mean)*(e3 - mean));

    return (v - c1 - c3) / pi2 - (e2 - mean)*(e2 - mean);
}

SEXP do_dtruncnorm(SEXP s_x, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_x   , x   , n_x);
  UNPACK_REAL_VECTOR(s_a   , a   , n_a);
  UNPACK_REAL_VECTOR(s_b   , b   , n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  
  n = MAX(MAX(MAX(n_x, n_a), MAX(n_b, n_mean)), n_sd);  
  ALLOC_REAL_VECTOR(s_ret, ret, n);

  for (i = 0; i < n; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cx = x[i];
    if (ca <= cx && cx <= cb) { /* In range: */
      const double cmean = mean[i % n_mean];
      const double csd = sd[i % n_sd];

      const double c1 = pnorm(ca, cmean, csd, TRUE, FALSE);
      const double c2 = pnorm(cb, cmean, csd, TRUE, FALSE);
      const double c3 = csd * (c2 - c1);      
      const double c4 = dnorm((cx-cmean)/csd, 0.0, 1.0, FALSE);
      ret[i] = c4 / c3;
    } else { /* Truncated: */
      ret[i] = 0.0;
    }
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}

SEXP do_ptruncnorm(SEXP s_q, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_q   , q   , n_q);
  UNPACK_REAL_VECTOR(s_a   , a   , n_a);
  UNPACK_REAL_VECTOR(s_b   , b   , n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  
  n = MAX(MAX(MAX(n_q, n_a), MAX(n_b, n_mean)), n_sd);  
  ALLOC_REAL_VECTOR(s_ret, ret, n);

  for (i = 0; i < n; ++i) {
    const double cq = q[i % n_q];
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    if (cq < ca) {
      ret[i] = 0.0;
    } else if (cq > cb) {
      ret[i] = 1.0;
    } else {
      const double cmean = mean[i % n_mean];
      const double csd = sd[i % n_sd];
      const double c1 = pnorm(cq, cmean, csd, TRUE, FALSE);
      const double c2 = pnorm(ca, cmean, csd, TRUE, FALSE);
      const double c3 = pnorm(cb, cmean, csd, TRUE, FALSE);
      ret[i] = (c1 - c2) / (c3 - c2);
    }
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}

SEXP do_etruncnorm(SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_a   , a   , n_a);
  UNPACK_REAL_VECTOR(s_b   , b   , n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);

  n = MAX(MAX(n_a, n_b), MAX(n_mean, n_sd));
  ALLOC_REAL_VECTOR(s_ret, ret, n);
  
  for (i = 0; i < n; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];

    if (R_FINITE(ca) && R_FINITE(cb)) {
	ret[i] = e_truncnorm(ca, cb, cmean, csd);
    } else if (R_NegInf == ca && R_FINITE(cb)) {
	ret[i] = e_righttruncnorm(cb, cmean, csd);
    } else if (R_FINITE(ca) && R_PosInf == cb) {
	ret[i] = e_lefttruncnorm(ca, cmean, csd);
    } else if (R_NegInf == ca && R_PosInf == cb) {
	ret[i] = cmean;
    } else {
	ret[i] = NA_REAL;
    }
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;  
}

SEXP do_vtruncnorm(SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_a   , a   , n_a);
  UNPACK_REAL_VECTOR(s_b   , b   , n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);

  n = MAX(MAX(n_a, n_b), MAX(n_mean, n_sd));
  ALLOC_REAL_VECTOR(s_ret, ret, n);
  
  for (i = 0; i < n; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];

    if (R_FINITE(ca) && R_FINITE(cb)) {
	ret[i] = v_truncnorm(ca, cb, cmean, csd);
    } else if (R_NegInf == ca && R_FINITE(cb)) {
	ret[i] = v_righttruncnorm(cb, cmean, csd);
    } else if (R_FINITE(ca) && R_PosInf == cb) {
	ret[i] = v_lefttruncnorm(ca, cmean, csd);
    } else if (R_NegInf == ca && R_PosInf == cb) {
	ret[i] = csd * csd;
    } else {
	ret[i] = NA_REAL;
    }
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}
