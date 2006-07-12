#include <stdio.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

#include "glmmml.h"
#include "fun.h"
#include "ghq.h"

P_fun *P;
G_fun *G;
Gprim_fun *Gprim;
Hprim_fun *Hprim;

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
	     int *boot,
	     double *predicted,
	     double *beta,
	     double *sigma,
	     double *loglik,
	     double *variance,
	     double *post_mode,
	     double *post_mean,
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
    int i, j, k, m, indx;
    
    double Fmin;
    double *b;
    double *gr;
    double **hessian;
    double *hess_vec;

    double rcond;
    double *det;
    int lwork;
    double *work;
    int job = 11;
    char *vmax;

    int bdim;

/* New in 0.28; bootstrapping: */
    int upper;

    int *conditional;
    int condi = 1;
/* This is done to prepare for having conditional as an input parameter */     
    conditional = &condi;

    if (*trace){
	Rprintf("************* Entering [glmmml] **************** \n");
	Rprintf(" p = %d\n\n", *p);
    }

    det = Calloc(2, double);

    if (*family == 0){
	P = &P_logit;
	G = &G_logit;
	Gprim = &Gprim_logit;
	Hprim = &Hprim_logit;
    }else if (*family == 1){
	P = &P_cloglog;
	G = &G_cloglog;
	Gprim = &Gprim_cloglog;
	Hprim = &Hprim_cloglog;
    }else if (*family == 2){
	P = &P_poisson;
	G = &G_poisson;
	Gprim = &Gprim_poisson;
	Hprim = &Hprim_poisson;
    }else{
	error("Unknown family\n");
    }
	
    abstol = *epsilon;
    reltol = abstol;

    bdim = *p + 1;
    lwork = 11 * bdim;
    work = Calloc(lwork, double);
    b = Calloc(bdim, double);
    gr = Calloc(bdim, double);

    hessian = Calloc(bdim, double *);
    hess_vec = Calloc(bdim * bdim, double);
    for (j = 0; j < bdim; j++) hessian[j] = hess_vec + j * bdim;

    /**** Build up 'ext' ********************/
    ext = Calloc(1, Exts);

    ext->family = *family; /* == 0 for binomial(logit) */

    ext->n = 0;
    for (i = 0; i < *n_fam; i++){
	ext->n += fam_size[i];
    }
    ext->p = *p;
    ext->cluster = cluster;
    /* Changed 2006-06-18; may have catastrophic consequences!! */
    ext->x = Calloc(ext->p, double *);
    for (i = 0; i < ext->p; i++){
	ext->x[i] = x + i * (ext->n);
    }
    /*** Note that ext->x is not "filled"; ***/ 
    /*** only points to the right place    ***/
    ext->offset = offset;
    ext->x_beta = Calloc(ext->n, double);
    ext->y = Calloc(ext->n, int); /* We cannot copy if bootstrapping! */
    for (i = 0; i < ext->n; i++)
	ext->y[i] = y[i];
    ext->n_fam = *n_fam;
    ext->fam_size = fam_size;
    ext->post_mode = Calloc(*n_fam, double); 
    ext->post_mean = Calloc(*n_fam, double); 
    ext->n_points = *n_points;
    ext->weights = Calloc(*n_points, double);
    ext->zeros = Calloc(*n_points, double);
    F77_CALL(ghq)(n_points, ext->zeros, ext->weights);  

/******* Done with 'ext' ***************/

    for (j = 0; i < *p; j++){
	b[j] = start_beta[j];
    }
    b[*p] = log(*start_sigma);

    mask = Calloc(bdim, int );
    for (i = 0; i < bdim; i++){
        mask[i] = 1;
    }

/* Note that this performs a minimum: (!!) */

    if (*method){
     /* *trace = 1;  Should be REMOVED!!!!!!!!!!!!!!!!!!!!!! */
	vmmin(bdim, b, &Fmin,
	      fun, fun1, *maxit, *trace,
	      mask, abstol, reltol, nREPORT,
	      ext, &fncount, &grcount, &fail);
	*convergence = (fail == 0);

	fun1(bdim, b, gr, ext);
	fun2(bdim, b, &Fmin, gr, hess_vec, ext);
	if(*trace){
	    Rprintf("Max log likelihood after vmmin: %f\n", -Fmin);
	    printf("beta: ");
	    for (i = 0; i < bdim; i++){
		Rprintf(" %f, ", b[i]);
	    }
	    Rprintf("\n");
	    printf("Gradients: ");
	    for (i = 0; i < bdim; i++){
		Rprintf(" %f, ", -gr[i]);
	    }
	    Rprintf("\n");
	    Rprintf("\n");
	    Rprintf("hessian:\n");
	    for (i = 0; i < bdim; i++){
		for (j = 0; j < bdim; j++)
		    Rprintf(" %f, ", hessian[i][j]);
		Rprintf("\n");
	    }
	}
	/* Let's avoid nr_opt! Just calculate the hessian and invert! */
/*
	nr_opt(bdim, b, &Fmin, mask, ext, *epsilon, nr_maxit, info, *trace); 
*/
	*loglik = Fmin;
	for (i = 0; i < *p; i++){
	    beta[i] = b[i];
	}
	*sigma = exp(b[*p]);
        F77_CALL(dpoco)(*hessian, &bdim, &bdim, &rcond, work, info);
        if (*info == 0){
            F77_CALL(dpodi)(*hessian, &bdim, &bdim, det, &job);
            for (m = 0; m < bdim; m++){
                for (k = 0; k < m; k++){
                    hessian[k][m] = hessian[m][k];
                }
            }
            indx = 0;
            for (m = 0; m < bdim; m++){
                for (k = 0; k < bdim; k++){
                    variance[indx] = hessian[m][k];
                    indx++;
                }
            }

        }else{
            Rprintf("info = %d\n", *info);
            warning("Hessian non-positive definite. No variance!");
        }

	if(*trace){
	    printf("Max log likelihood: %f\n", -Fmin);
	    printf("Beta: ");
	    for (i = 0; i < bdim; i++){
		printf(" %f, ", b[i]);
	    }
	    printf("\n");
	    printf("Gradients: ");
	    for (i = 0; i < bdim; i++){
		printf(" %f, ", gr[i]);
	    }
	    printf("\n");
	    printf("hessian:\n");
	    for (i = 0; i < bdim; i++){
		for (j = 0; j < bdim; j++)
		    printf(" %f, ", hessian[i][j]);
		printf("\n");
	    }
	}
    }else{
	error("Nelder-Mead not implemented (yet)\n");
    }

    frail_fun(bdim, b, ext);
    mu_fun(bdim, b, mu, ext);

    for (i = 0; i < ext->n_fam; i++){
	post_mode[i] = ext->post_mode[i];
	post_mean[i] = ext->post_mean[i];
    }

    if (*boot > 0){
/************** Bootstrapping starts *****************************/
	upper = 0;
	GetRNGstate();
	for (i = 0; i < *boot; i++){
	    if (*trace){
		if ((i / 10) * 10 == i)
		    printf("********************* Replicate No. %d\n", i);
	    }
	 
	    if (*family <= 1){ /* Bernoulli */
		for (j = 0; j < ext->n; j++)
		    ext->y[j] = rbinom(1, predicted[j]);
	    }else{
		for (j = 0; j < ext->n; j++) /* Poisson */
		    ext->y[j] = rpois(predicted[j]);
	    }
	
/* Restore beta as start values: */
	    for ( j = 0; j < *p; j++) b[j] = beta[j];
	    if (*trace){
/*
		Rprintf("Sampled values by cluster:\n");
		for (i = 0; i < ext->n; i++){
		    Rprintf("y[%d] = %d, cluster[%d] = %d\n", 
			    i, ext->y[i], i, ext->cluster[i]);
		}
*/
		Rprintf("Start value to vmmin: %f\n", fun(*p, b, ext)); 
	    }
	    vmax = vmaxget();
	    vmmin(*p, b, &Fmin,
		  fun, fun1, *maxit, *trace,
		  mask, abstol, reltol, nREPORT,
		  ext, &fncount, &grcount, &fail);
	    vmaxset(vmax);
	    *convergence = (fail == 0);
	    if (!convergence){
		Rprintf("No converhǵence...\n");
	    }
	    boot_log[i] = -Fmin;
	    if (*trace){
		Rprintf("boot_log[%d] = %f; loglik = %f\n", i, -Fmin, *loglik);
		Rprintf("beta[0] = %f\n", b[0]);
	    }
	    if (-Fmin >= *loglik) upper++;
	}
	
	if (*boot) *boot_p = (double)upper / (double)*boot;
	else *boot_p = 1.0;
	
	PutRNGstate();
	
    }
    
    Free(mask);

    Free(ext->zeros);
    Free(ext->weights);

    Free(ext->post_mean);
    Free(ext->post_mode);
    Free(ext->x_beta);
    Free(ext->y);
    Free(ext->x);
    Free(ext);

    Free(hessian);
    Free(hess_vec);
    Free(gr);
    Free(b);
    Free(work);
    Free(det);
}
