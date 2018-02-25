#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>  // optional

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern SEXP do_dtruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP do_ptruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP do_qtruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP do_rtruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP do_etruncnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP do_vtruncnorm(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef R_CallDef[] = {
    CALLDEF(do_dtruncnorm, 5),
    CALLDEF(do_ptruncnorm, 5),
    CALLDEF(do_qtruncnorm, 5),
    CALLDEF(do_rtruncnorm, 5),
    CALLDEF(do_etruncnorm, 4),
    CALLDEF(do_vtruncnorm, 4),
    {NULL, NULL, 0}
};

void attribute_visible R_init_truncnorm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
