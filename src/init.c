#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void glmm_boot(void *, void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *, void *,
		      void *, void *, void *, void *, void *, void *);
extern void glmm_boot0(void *, void *, void *, void *, void *, void *,
		       void *, void *, void *, void *, void *, void *,
		       void *, void *, void *, void *);
extern void glmm_ml(void *, void *, void *, void *, void *, void *,
		    void *, void *, void *, void *, void *, void *,
		    void *, void *, void *, void *, void *, void *,
		    void *, void *, void *, void *, void *, void *,
		    void *, void *, void *, void *, void *, void *,
		    void *, void *);

/* .Fortran calls */
extern void F77_NAME(ghq)(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"glmm_boot",  (DL_FUNC) &glmm_boot,  24},
    {"glmm_boot0", (DL_FUNC) &glmm_boot0, 16},
    {"glmm_ml",    (DL_FUNC) &glmm_ml,    32},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"ghq", (DL_FUNC) &F77_NAME(ghq), 4},
    {NULL, NULL, 0}
};

void R_init_glmmML(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
