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

SEXP dtruncnorm(SEXP s_x, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
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


SEXP ptruncnorm(SEXP s_q, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
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


SEXP rtruncnorm(SEXP s_n, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, nn;
  UNPACK_INT_VECTOR(s_n   , n   , n_n);
  UNPACK_REAL_VECTOR   (s_a   , a   , n_a);
  UNPACK_REAL_VECTOR   (s_b   , b   , n_b);
  UNPACK_REAL_VECTOR   (s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR   (s_sd  , sd  , n_sd);
  
  nn = MAX(MAX(MAX(n[0], n_n), MAX(n_a, n_b)), MAX(n_mean, n_sd));
  
  ALLOC_REAL_VECTOR(s_ret, ret, nn);
  GetRNGstate();
  for (i = 0; i < nn; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    double tmp;
    /* Rejection sampling until we generate an observation that is not
     * truncated. This can be slow for narrow and extreme intervals
     * [a,b].
     */
    while (1) { 
      tmp = rnorm(cmean, csd);
      if (ca <= tmp && tmp <= cb)
	break;
    }
    ret[i] = tmp;
  }
  PutRNGstate();
  UNPROTECT(1); /* s_ret */
  return s_ret;
}


SEXP etruncnorm(SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_a   , a   , n_a);
  UNPACK_REAL_VECTOR(s_b   , b   , n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);

  n = MAX(MAX(n_a, n_b), MAX(n_mean, n_sd));
  ALLOC_REAL_VECTOR(s_ret, ret, n);
  
  for (i = 0; i < n_a; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    
    const double c1 = dnorm((ca-cmean)/csd, 0.0, 1.0, FALSE);
    const double c2 = dnorm((cb-cmean)/csd, 0.0, 1.0, FALSE);
    const double C1 = pnorm(ca, cmean, csd, TRUE, FALSE);
    const double C2 = pnorm(cb, cmean, csd, TRUE, FALSE);   
    ret[i] = cmean + csd * ((c1 - c2) / (C2 - C1));
  } 
  UNPROTECT(1); /* s_ret */
  return s_ret;
}


SEXP vtruncnorm(SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_a   , a   , n_a);
  UNPACK_REAL_VECTOR(s_b   , b   , n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);

  n = MAX(MAX(n_a, n_b), MAX(n_mean, n_sd));
  ALLOC_REAL_VECTOR(s_ret, ret, n);
  
  for (i = 0; i < n_a; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    const double am = (ca - cmean)/csd;
    const double bm = (cb - cmean)/csd;
    const double c1 = dnorm(am, 0.0, 1.0, FALSE);
    const double c2 = dnorm(bm, 0.0, 1.0, FALSE);
    const double C1 = pnorm(ca, cmean, csd, TRUE, FALSE);
    const double C2 = pnorm(cb, cmean, csd, TRUE, FALSE);   
  
    const double v  = csd * csd;
    
    const double d = C2 - C1;
    const double m1 = (c1 - c2)/d;
    const double m2 = (am*c1 - bm*c2)/d;
    ret[i] = (1.0 + m2 + m1*m1)*v;
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}
