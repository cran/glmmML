#ifndef GLMM_BOOT_H
#define GLMM_BOOT_H

#ifndef MATHLIB_STANDALONE
#include <R.h>
/* #include <R_ext/RS.h> */
#else
#include "testa.h"
#endif

typedef struct
{
    int family;
    int n;               /* = sum _0^(n_fam - 1) fam_size[i] */
    int p;               /* No. of covariates _excluding_    */ 
                         /* a constant                       */
    int *fam_size;       /* n_fam-vector.                    */
    int n_fam;           /* No. of families.                 */
    int *success;        /* n_fam-vector: Family totals      */
    double **x;          /* n x p 'matrix'                   */
    double *x_beta;      /* n: linear predictor              */
    double *pred;        /* n: logit(x_beta + gamma)         */           
    double *offset;      /* n                                */
    int *ki;             /* n                                */
    int *cluster;        /* n                                */
    int *fam_out;        /* n_fam: 0 == include,             */ 
			 /*       -1 all failures,           */
			 /*       +1 all successes           */
    double *gamma;       /* n_fam                            */
    double *gr;          /* bdim                             */
    double *hessian;    /* p x p                            */
    int *y;              /* n-vector: binary response (0-1). */
}
Extb;

void glmm_boot(int *family,
	       int *method,
	       int *p, 
	       double *start_beta,
	       int *cluster,
	       double *x, /* Now p x (\sum_1_{n_fam} fam_size[i]) */
	       int *y,
	       double *offset,
	       int *fam_size,
	       int *n_fam,
	       double *epsilon,
	       int *maxit,
	       int *trace,
	       int *boot,
	       double *beta,
	       double *loglik,
	       double *variance,
	       double *frail,
	       double *boot_p,
	       double *boot_log,
	       int *convergence);
    
#endif

