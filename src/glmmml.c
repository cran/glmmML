#include <stdio.h>
#include <R_ext/Applic.h>

#include "glmmml.h"
#include "fun.h"
#include "ghq.h"

P_fun *P;
G_fun *G;
Gprim_fun *Gprim;

void glmm_ml(int *family,
	     int *method,
	     int *p, 
	     double *start_beta,
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
	     double *beta,
	     double *sigma,
	     double *loglik,
	     double *variance,
	     double *frail,
	     double *mu,
	     int *convergence){
    
#ifndef MATHLIB_STANDALONE
    double abstol;
    double reltol;
    int nREPORT = 1;
    int fncount;
    int grcount;
    int fail;
    int *mask;
#endif
    Exts *ext;
    int i;
    
    double Fmin;
    double *b;
    double *gr;
    int bdim;
    int nr_maxit = 0;

    if (*family == 0){
	P = &P_logit;
	G = &G_logit;
	Gprim = &Gprim_logit;
    }else if (*family == 1){
	P = &P_cloglog;
	G = &G_cloglog;
	Gprim = &Gprim_cloglog;
    }else if (*family == 2){
	P = &P_poisson;
	G = &G_poisson;
	Gprim = &Gprim_poisson;
    }else{
	error("Unknown family\n");
    }
	
    abstol = *epsilon;
    reltol = abstol;

    bdim = *p + 1;

    b = Calloc(bdim, double);
    gr = Calloc(bdim, double);

    ext = Calloc(1, Exts);

    ext->family = *family; /* == 0 for binomial(logit) */

    for (i = 0; i < *p; i++){
	b[i] = start_beta[i];
    }

    if (*n_points == 1){/* Put in log(sigma) */
	b[*p] = 0.0;
    }else{
	b[*p] = log(*start_sigma);
    }
    ext->p = *p;

    ext->offset = offset;
    ext->weights = Calloc(*n_points, double);
    ext->zeros = Calloc(*n_points, double);

    mask = Calloc(bdim, int );
    for (i = 0; i < bdim; i++){
        mask[i] = 1;
    }

    F77_CALL(ghq)(n_points, ext->zeros, ext->weights);  

    ext->n = 0;
    for (i = 0; i < *n_fam; i++){
	ext->n += fam_size[i];
    }

    ext->x_beta = Calloc(ext->n, double); 
    ext->gr = gr;
    ext->hessian = variance;
    ext->x = x;
    ext->y = y;
    ext->n_fam = *n_fam;
    ext->fam_size = fam_size;
    ext->n_points = *n_points;

    if (ext->n_points == 1){
	mask[ext->p] = 0;
	b[ext->p] = 0.0;
    }
/* Note that this performs a minimum: (!!) */

    if (*method){

	vmmin(bdim, b, &Fmin,
	      fun, fun1, *maxit, *trace,
	      mask, abstol, reltol, nREPORT,
	      ext, &fncount, &grcount, &fail);
	*convergence = (fail == 0);

	fun1(bdim, b, gr, ext);
	if(*trace){
	    Rprintf("Max log likelihood after vmmin: %f\n", -Fmin);
	    printf("Gradients: ");
	    for (i = 0; i < bdim; i++){
		Rprintf(" %f, ", -ext->gr[i]);
	    }
	    Rprintf("\n");
	}

	nr_opt(bdim, b, &Fmin, mask, ext, *epsilon, nr_maxit, *trace);
	*loglik = Fmin;
	for (i = 0; i < *p; i++){
	    beta[i] = b[i];
	}
	*sigma = exp(b[*p]);
	if(*trace){
	    printf("Max log likelihood after Newton-Raphson: %f\n", -Fmin);
	    printf("Gradients: ");
	    for (i = 0; i < bdim; i++){
		printf(" %f, ", ext->gr[i]);
	    }
	    printf("\n");
	}
    }else{
	error("Nelder-Mead not implemented (yet)\n");
    }
    frail_fun(bdim, b, frail, ext);
    mu_fun(bdim, b, mu, ext);

    Free(ext->x_beta);

    Free(mask);

    Free(gr);
    Free(ext->zeros);
    Free(ext->weights);
    Free(ext);
    Free(gr);
    Free(b);
}
