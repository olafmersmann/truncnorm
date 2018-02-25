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
#include "zeroin.h"

#define ALLOC_REAL_VECTOR(S, D, N)                                             \
  SEXP S;                                                                      \
  PROTECT(S = allocVector(REALSXP, N));                                        \
  double *D = REAL(S);

#ifndef MAX
#define MAX(A, B) ((A > B) ? (A) : (B))
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
  const double phi_a = dnorm(alpha, 0.0, 1.0, TRUE);
  const double Phi_a = pnorm(alpha, 0.0, 1.0, FALSE, TRUE);
  double res = mean + sd * exp(phi_a - Phi_a);
  return res;
}

static R_INLINE double e_truncnorm(double a, double b, double mean, double sd) {
  /* Special case numerically instable case when (a, b) is far away from the
   * center of mass. */
  if (b < mean - 6.0 * sd || a > mean + 6.0 * sd)
    return (a + b) / 2.0;

  double delta_phi = 0.0, delta_Phi = 0.0;

  const double alpha = (a - mean) / sd;
  const double beta = (b - mean) / sd;
  const double phi_a = dnorm(alpha, 0.0, 1.0, TRUE);
  const double Phi_a = pnorm(alpha, 0.0, 1.0, TRUE, TRUE);
  const double phi_b = dnorm(beta, 0.0, 1.0, TRUE);
  const double Phi_b = pnorm(beta, 0.0, 1.0, TRUE, TRUE);

  if (phi_b < phi_a) {
    delta_phi = logspace_sub(phi_a, phi_b);
  } else {
    sd = -sd;
    delta_phi = logspace_sub(phi_b, phi_a);
  }

  if (Phi_b > Phi_a) {
    sd = -sd;
    delta_Phi = logspace_sub(Phi_b, Phi_a);
  } else {
    delta_Phi = logspace_sub(Phi_a, Phi_b);
  }
  return mean + sd * -exp(delta_phi - delta_Phi);
}

static R_INLINE double e_righttruncnorm(double b, double mean, double sd) {
  const double beta = (b - mean) / sd;
  const double phi_b = dnorm(beta, 0.0, 1.0, TRUE);
  const double Phi_b = pnorm(beta, 0.0, 1.0, TRUE, TRUE);
  return mean + sd * -exp(phi_b - Phi_b);
}

static R_INLINE double v_lefttruncnorm(double a, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  const double phi_a = dnorm(alpha, 0.0, 1.0, FALSE);
  const double Phi_a = pnorm(alpha, 0.0, 1.0, TRUE, FALSE);
  const double lambda = phi_a / (1.0 - Phi_a);
  return (sd * sd * (1.0 - lambda * (lambda - alpha)));
}

static R_INLINE double v_righttruncnorm(double b, double mean, double sd) {
  return (v_lefttruncnorm(-b, -mean, sd));
}

static R_INLINE double v_truncnorm(double a, double b, double mean, double sd) {
  /* Special case numerically instable cases. These arise when (a, b) is far
   * away from mean +/- 6*sd */
  if (b < mean - 6.0 * sd || a > mean + 6.0 * sd)
    return 1.0 / 12 * (b - a) * (b - a);

  const double v = sd * sd;
  const double pi1 = pnorm(a, mean, sd, TRUE, FALSE);
  const double pi2 =
      pnorm(b, mean, sd, TRUE, FALSE) - pnorm(a, mean, sd, TRUE, FALSE);
  const double pi3 = pnorm(b, mean, sd, FALSE, FALSE); /* 1 - F(b) */
  const double e1 = e_righttruncnorm(a, mean, sd);
  const double e2 = e_truncnorm(a, b, mean, sd);
  const double e3 = e_lefttruncnorm(b, mean, sd);

  const double v1 = v_righttruncnorm(a, mean, sd);
  const double v3 = v_lefttruncnorm(b, mean, sd);

  const double c1 = pi1 * (v1 + (e1 - mean) * (e1 - mean));
  const double c3 = pi3 * (v3 + (e3 - mean) * (e3 - mean));
  const double cd = pi2 - (e2 - mean) * (e2 - mean);
  return (v - c1 - c3) / pi2 - (e2 - mean) * (e2 - mean);
}

static R_INLINE double ptruncnorm(const double q, const double a,
                                  const double b, const double mean,
                                  const double sd) {
  if (q < a) {
    return 0.0;
  } else if (q > b) {
    return 1.0;
  } else {
    const double c1 = pnorm(q, mean, sd, TRUE, FALSE);
    const double c2 = pnorm(a, mean, sd, TRUE, FALSE);
    const double c3 = pnorm(b, mean, sd, TRUE, FALSE);
    return (c1 - c2) / (c3 - c2);
  }
}

typedef struct { double a, b, mean, sd, p; } qtn;

/* qtmin - helper function to calculate quantiles of the truncated
 *   normal distribution.
 *
 * The root of this function is the desired quantile, given that *p
 * defines a truncated normal distribution and the desired
 * quantile. This function increases monotonically in x and is
 * positive for x=a and negative for x=b if 0 < p < 1.
 */
double qtmin(double x, void *p) {
  qtn *t = (qtn *)p;
  return ptruncnorm(x, t->a, t->b, t->mean, t->sd) - t->p;
}

SEXP do_dtruncnorm(SEXP s_x, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_x, x, n_x);
  UNPACK_REAL_VECTOR(s_a, a, n_a);
  UNPACK_REAL_VECTOR(s_b, b, n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd, sd, n_sd);

  n = MAX(MAX(MAX(n_x, n_a), MAX(n_b, n_mean)), n_sd);
  ALLOC_REAL_VECTOR(s_ret, ret, n);

  for (i = 0; i < n; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cx = x[i % n_x];
    if (ca <= cx && cx <= cb) { /* In range: */
      const double cmean = mean[i % n_mean];
      const double csd = sd[i % n_sd];

      const double c1 = pnorm(ca, cmean, csd, TRUE, FALSE);
      const double c2 = pnorm(cb, cmean, csd, TRUE, FALSE);
      const double c3 = csd * (c2 - c1);
      const double c4 = dnorm((cx - cmean) / csd, 0.0, 1.0, TRUE);
      if (!isfinite(log(c3))) {
        ret[i] = 1.0 / (cb - ca);
      } else {
        ret[i] = exp(c4 - log(c3));
      }
    } else { /* Truncated: */
      ret[i] = 0.0;
    }
    R_CheckUserInterrupt();
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}

SEXP do_ptruncnorm(SEXP s_q, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_q, q, n_q);
  UNPACK_REAL_VECTOR(s_a, a, n_a);
  UNPACK_REAL_VECTOR(s_b, b, n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd, sd, n_sd);

  n = MAX(MAX(MAX(n_q, n_a), MAX(n_b, n_mean)), n_sd);
  ALLOC_REAL_VECTOR(s_ret, ret, n);

  for (i = 0; i < n; ++i) {
    const double cq = q[i % n_q];
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    ret[i] = ptruncnorm(cq, ca, cb, cmean, csd);
    R_CheckUserInterrupt();
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}

SEXP do_qtruncnorm(SEXP s_p, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  qtn t;
  double tol;
  int maxit;
  UNPACK_REAL_VECTOR(s_p, p, n_p);
  UNPACK_REAL_VECTOR(s_a, a, n_a);
  UNPACK_REAL_VECTOR(s_b, b, n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd, sd, n_sd);

  n = MAX(MAX(MAX(n_p, n_a), MAX(n_b, n_mean)), n_sd);
  ALLOC_REAL_VECTOR(s_ret, ret, n);

  for (i = 0; i < n; ++i) {
    const double cp = p[i % n_p];
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    if (cp == 0.0) {
      ret[i] = ca;
    } else if (cp == 1.0) {
      ret[i] = cb;
    } else if (cp < 0.0 || cp > 1.0) {
      ret[i] = R_NaN;
    } else if (ca == R_NegInf && cb == R_PosInf) {
      ret[i] = qnorm(cp, cmean, csd, TRUE, FALSE);
    } else {
      /* We need to possible adjust ca and cb for R_zeroin(),
       * because R_zeroin() requires finite bounds and ca or cb (but
       * not both, see above) may be infinite. In that case, we use
       * a simple stepping out procedure to find a lower or upper
       * bound from which to begin the search.
       */
      double lower = ca, upper = cb;
      if (lower == R_NegInf) {
        lower = -1;
        while (ptruncnorm(lower, ca, cb, cmean, csd) - cp >= 0)
          lower *= 2.0;
      } else if (upper == R_PosInf) {
        upper = 1;
        while (ptruncnorm(upper, ca, cb, cmean, csd) - cp <= 0)
          upper *= 2.0;
      }
      t.a = ca;
      t.b = cb;
      t.mean = cmean;
      t.sd = csd;
      t.p = cp;
      maxit = 200;
      tol = 0.0; /* Set tolerance! */
      ret[i] = truncnorm_zeroin(lower, upper, qtmin(lower, &t),
                                qtmin(upper, &t), &qtmin, &t, &tol, &maxit);
    }
    R_CheckUserInterrupt();
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}

SEXP do_etruncnorm(SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_a, a, n_a);
  UNPACK_REAL_VECTOR(s_b, b, n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd, sd, n_sd);

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
    R_CheckUserInterrupt();
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}

SEXP do_vtruncnorm(SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, n;
  UNPACK_REAL_VECTOR(s_a, a, n_a);
  UNPACK_REAL_VECTOR(s_b, b, n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd, sd, n_sd);

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
    R_CheckUserInterrupt();
  }
  UNPROTECT(1); /* s_ret */
  return s_ret;
}
