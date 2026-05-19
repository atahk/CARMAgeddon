#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call entry points */
extern SEXP monocarptr__new(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                            SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                            SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP monocarptr__print(SEXP, SEXP);
extern SEXP monocarptr__fun(SEXP, SEXP);
extern SEXP monocarptr__grad(SEXP, SEXP);
extern SEXP monocarptr__parallelgrad(SEXP, SEXP);
extern SEXP monocarptr__partial(SEXP, SEXP);
extern SEXP monocarptr__hess(SEXP, SEXP);
extern SEXP monocarptr__hist(SEXP, SEXP);
extern SEXP monocarptr__histrand(SEXP, SEXP, SEXP);
extern SEXP monocarptr__endopt(SEXP);
extern SEXP ousim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                  SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"monocarptr__new",          (DL_FUNC) &monocarptr__new,          28},
    {"monocarptr__print",        (DL_FUNC) &monocarptr__print,         2},
    {"monocarptr__fun",          (DL_FUNC) &monocarptr__fun,           2},
    {"monocarptr__grad",         (DL_FUNC) &monocarptr__grad,          2},
    {"monocarptr__parallelgrad", (DL_FUNC) &monocarptr__parallelgrad,  2},
    {"monocarptr__partial",      (DL_FUNC) &monocarptr__partial,       2},
    {"monocarptr__hess",         (DL_FUNC) &monocarptr__hess,          2},
    {"monocarptr__hist",         (DL_FUNC) &monocarptr__hist,          2},
    {"monocarptr__histrand",     (DL_FUNC) &monocarptr__histrand,      3},
    {"monocarptr__endopt",       (DL_FUNC) &monocarptr__endopt,        1},
    {"ousim",                    (DL_FUNC) &ousim,                    22},
    {NULL, NULL, 0}
};

void R_init_CARMAgeddon(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);  /* optional but recommended */
}
