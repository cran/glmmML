#ifndef LOGISTIC_H
#define LOGISTIC_H

#ifndef MATHLIB_STANDALONE
#include <R.h>
#include <Rmath.h>
/* #include <R_ext/RS.h> */
#else
#include "testa.h"
#endif

typedef struct
{
    int family;          /* 0 = Bernoulli, logit link        */ 
                         /* 1 =            cloglog link      */
                         /* 2 = Poisson, log link            */

    int n;               /* = sum _0^(n_fam - 1) fam_size[i] */
    int p;               /* No. of covariates _including_    */ 
                         /* the constant (if any)            */
    int *cluster;        /* n                                */
    double **x;          /* n x p NOTE: check carefully **!! */
    double *offset;      /* n                                */
    int *ki;             /* n                                */
    double *x_beta;      /* <- x %*% beta                    */
    double *gr;          /* p + 1                            */
    double *hessian;     /* (p+1) x (p+1)                    */
    int *y;              /* n-vector: binary response (0-1). */
    int n_fam;           /* No. of families.                 */
    int *fam_size;       /* n_fam-vector.                    */
    int n_points;        /* No. of Gauss-Hermite points.     */
    double *weights;     /* n_points-vector.                 */
    double *zeros;       /* n_points-vector.                 */
}
Exts;

void glmm_ml(int *family,
             int *method,
             int *p, 
             double *start_beta,
	     int *cluster,
             double *start_sigma,
             double *x, /* Now p x (\sum_1_{n_fam} fam_size[i]) */
             int *y,
             double *offset,
             int *fam_size,
             int *n_fam,
             int *n_points, /* No. of pts in Gauss-Hermite quadrature */
             double *epsilon,
             int *maxit,
             int *trace,
	     int *boot,
	     double *predicted,
	     double *beta,
             double *sigma,
             double *loglik,
             double *variance,
             double *frail,
             double *mu,
	     double *boot_p,
	     double* boot_log,
             int *convergence,
	     int *info);

#endif
