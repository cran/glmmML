#ifndef LOGISTIC_H
#define LOGISTIC_H

#ifndef MATHLIB_STANDALONE
#include <R.h>
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
    double *x;           /* n x p                            */
    double *offset;      /* n                                */
    double *x_beta;      /* <- x %*% beta                    */
    double *gr;          /* p + 1                            */
    double *hessian;     /* (p+1) x (p+1)                    */
    int *y;              /* n-vector: binary response (0-1). */
    int *id;             /* n-vector.                        */
    int n_fam;           /* No. of families.                 */
    int *fam_size;       /* n_fam-vector.                    */
    int n_points;        /* No. of Gauss-Hermite points.     */
    double *weights;     /* n_points-vector.                 */
    double *zeros;       /* n_points-vector.                 */

}
Exts;

void logistic(int *method,
	      int *bdim, 
	      double *beta, /* beta[bdim - 1] == log(sigma) */
	      double *x,
	      int *y,
	      int *fam_size,
	      int *n_fam,
	      int *n_points, /* No. of pts in Gauss-Hermite quadrature */
	      double *Fmin,
	      double *hessian, /* bdim x bdim */
	      int *convergence);
    
#endif

