#include <stdio.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

#include "glmmml.h"
#include "fun.h"
#include "ghq.h"

P_fun *P;
G_fun *G;
Gprim_fun *Gprim;

static void permute(int n, int *y, int *x)
{
/* Result in y; x is just a "work area" */

    int i, j, k;

    k = n; /* Eg, "wt sampling k-of-n" */

    for (i = 0; i < n; i++)
	x[i] = i;
    for (i = 0; i < k; i++) {
	j = n * unif_rand();
/*	y[i] = x[j] + 1; */
	y[i] = x[j];
	x[j] = x[--n];
    }
}

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
	     int * boot,
	     double *predicted,
	     double *beta,
	     double *sigma,
	     double *loglik,
	     double *variance,
	     double *frail,
	     double *mu,
	     double *boot_p,
	     double *boot_log,
	     int *convergence,
	     int *info){
    
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
    int i, j;
    
    double Fmin;
    double *b;
    double *gr;
    int bdim;
    int nr_maxit = 0;
/* New in 0.28; bootstrapping: */
    int *ki;
    int*ki_tmp;
    int upper;
    char *vmax;
    int *conditional;
    int condi = 1;
/* This is done to prepare for having conditional as an input parameter */     
    conditional = &condi;

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

    for (i = 0; i < *p; i++){
	b[i] = start_beta[i];
    }

    if (*n_points == 1){/* Put in log(sigma) */
	b[*p] = 0.0;
    }else{
	b[*p] = log(*start_sigma);
    }

    mask = Calloc(bdim, int );
    for (i = 0; i < bdim; i++){
        mask[i] = 1;
    }

    /******** Filling in 'ext': ***********/

    ext = Calloc(1, Exts);

    ext->family = *family; /* == 0 for binomial(logit) */

    ext->n = 0;
    for (i = 0; i < *n_fam; i++){
	ext->n += fam_size[i];
    }

    ext->p = *p;

    ext->cluster = Calloc(ext->n, int);
    for (i = 0; i < ext->n; i++)
	ext->cluster[i] = cluster[i];

    /* ext->x = x; */
    /* Changed 2006-06-18; may have catastrophic consequences!! */
    ext->x = Calloc(ext->n, double *);
    for (i = 0; i < ext->n; i++){
	ext->x[i] = x + i * (ext->p);
    }
    /*** Note that ext->x is not "filled"; ***/ 
    /*** only points to the right place    ***/

    ext->offset = Calloc(ext->n, double);
    for (i = 0; i < ext->n; i++){
	ext->offset[i] = offset[i];
    }

    ext->ki = Calloc(ext->n, int);

    ext->x_beta = Calloc(ext->n, double); 

    ext->gr = gr;

    ext->hessian = variance;

    ext->y = Calloc(ext->n, int);
    for (i = 0; i < ext->n; i++) ext->y[i] = y[i];

    ext->n_fam = *n_fam;
    ext->fam_size = fam_size;

    ext->n_points = *n_points;

    ext->weights = Calloc(*n_points, double);
    ext->zeros = Calloc(*n_points, double);

    F77_CALL(ghq)(n_points, ext->zeros, ext->weights);  

    if (ext->n_points == 1){
	mask[ext->p] = 0;
	b[ext->p] = 0.0;
    }

    /******** Done filling 'ext' *************/

    ki = ext->ki;
    for (i = 0; i < ext->n; i++){
	ki[i] = i;
    }
    ki_tmp = Calloc(ext->n, int);


/* Note that this performs a minimum: (!!) */

    if (*method){

	vmmin(bdim, b, &Fmin,
	      fun, fun1, *maxit, *trace,
	      mask, abstol, reltol, nREPORT,
	      ext, &fncount, &grcount, &fail);
	*convergence = (fail == 0);

	fun1(bdim, b, gr, ext);
	fun2(bdim, b, &Fmin, ext->gr, ext->hessian, ext);
	if(*trace){
	    Rprintf("Max log likelihood after vmmin: %f\n", -Fmin);
	    printf("beta: ");
	    for (i = 0; i < bdim; i++){
		Rprintf(" %f, ", b[i]);
	    }
	    Rprintf("\n");
	    printf("Gradients: ");
	    for (i = 0; i < bdim; i++){
		Rprintf(" %f, ", -ext->gr[i]);
	    }
	    Rprintf("\n");
	    printf("\n");
	    printf("hessian:\n");
	    for (i = 0; i < bdim; i++){
		for (j = 0; j < bdim; j++)
		    printf(" %f, ", ext->hessian[i * bdim + j]);
		printf("\n");
	    }
	}
	/* Let's avoid nr_opt! Just calculate the hessian and invert! */
	nr_opt(bdim, b, &Fmin, mask, ext, *epsilon, nr_maxit, info, *trace); 
	*loglik = Fmin;
	for (i = 0; i < *p; i++){
	    beta[i] = b[i];
	}
	*sigma = exp(b[*p]);
	if(*trace){
	    printf("Max log likelihood after Newton-Raphson: %f\n", -Fmin);
	    printf("Beta: ");
	    for (i = 0; i < bdim; i++){
		printf(" %f, ", b[i]);
	    }
	    printf("\n");
	    printf("Gradients: ");
	    for (i = 0; i < bdim; i++){
		printf(" %f, ", ext->gr[i]);
	    }
	    printf("\n");
	    printf("hessian:\n");
	    for (i = 0; i < bdim; i++){
		for (j = 0; j < bdim; j++)
		    printf(" %f, ", ext->hessian[i * bdim + j]);
		printf("\n");
	    }
	}
    }else{
	error("Nelder-Mead not implemented (yet)\n");
    }
    frail_fun(bdim, b, frail, ext);
    mu_fun(bdim, b, mu, ext);

    if (boot > 0){
/************** Bootstrapping starts *****************************/
	upper = 0;
	GetRNGstate();

	for (i = 0; i < *boot; i++){
	    if (*trace){
		if ((i / 10) * 10 == i)
		    printf("********************* Replicate No. %d\n", i);
	    }
	    if (*conditional){
/* Conditional bootstrap */
		permute(ext->n, ki, ki_tmp);
		for (j = 0; j < ext->n; j++){
		    ext->y[j] = y[ki[j]];
		    ext->x[j] = x + ki[j] * (ext->p);
		    ext->offset[j] = offset[ki[j]];
		    ext->cluster[j] = cluster[ki[j]];
		}
	    }else{
		if (*family <= 1){ /* Bernoulli */
		    for (j = 0; j < ext->n; j++)
			ext->y[j] = rbinom(1, predicted[j]);
		}else{
		    for (j = 0; j < ext->n; j++) /* Poisson */
			ext->y[j] = rpois(predicted[j]);
		}
	    }
/* Restore beta as start values: */
	    for ( j = 0; j < *p; j++) b[j] = beta[j];
	    
	    vmax = vmaxget();
	    vmmin(*p, b, &Fmin,
		  fun, fun1, *maxit, *trace,
		  mask, abstol, reltol, nREPORT,
		  ext, &fncount, &grcount, &fail);
	    vmaxset(vmax);
	    *convergence = (fail == 0);
	    boot_log[i] = -Fmin;
	    if (-Fmin >= *loglik) upper++;
	}
	
	if (*boot) *boot_p = (double)upper / (double)*boot;
	else *boot_p = 1.0;
	
	PutRNGstate();
	
    }
    
    free(ki_tmp);

    Free(ext->zeros);
    Free(ext->weights);
    Free(ext->y);
    Free(ext->x_beta);
    Free(ext->ki);
    Free(ext->offset);
    Free(ext->x);
    Free(ext->cluster);
    Free(ext);

    Free(mask);

    Free(gr);
    Free(b);
}
