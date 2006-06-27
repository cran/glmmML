#include <stdio.h>
#include <R_ext/Linpack.h>

#include "glmmboot.h"
#include "bfun.h"
#include "fun.h"

extern P_fun *P;
extern G_fun *G;
extern Gprim_fun *Gprim;

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

void glmm_boot0(int *family,
		int *method,
		int *cluster,
		int *y,
		double *offset,
		int *fam_size,
		int *n_fam,
		int *conditional,
		int *trace,
		int *boot,
		double *predicted,
		double *loglik,
		double *frail,
		double *boot_p,
		double *boot_log,
		int *convergence){
    
#ifndef MATHLIB_STANDALONE
    double abstol;
    double reltol;
#endif
    Extb *ext;
    int i;
    int j;

    int *ki;
    int *ki_tmp;

    double Fmin;
    double *b = NULL;
    int upper;

    GetRNGstate(); /* For random number generation */

/*    vmax1 = vmaxget(); */

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

    abstol = 0.00000001;
    reltol = abstol;

    ext = Calloc(1, Extb);
/************************ Fill in ext: *****************/
    ext->family = *family; /* == 0 for binomial(logit) */

    ext->n = 0;
    for (i = 0; i < *n_fam; i++){
	ext->n += fam_size[i];
    }
    ext->p = 0;
    ext->fam_size = fam_size;
    ext->n_fam = *n_fam;
    ext->success = Calloc(*n_fam, int);
    ext->fam_out = Calloc(ext->n_fam, int);
    ext->x_beta = Calloc(ext->n, double);

    ext->x = NULL;

    ext->pred = Calloc(ext->n, double);
    ext->offset = Calloc(ext->n, double);
    for (i = 0; i < ext->n; i++){
	ext->offset[i] = offset[i];
    }
    
    ext->ki = Calloc(ext->n, int);
    ext->cluster = Calloc(ext->n, int);
    for (i = 0; i < ext->n; i++)
	ext->cluster[i] = cluster[i];
    ext->gamma = Calloc(ext->n_fam, double);
    ext->y = Calloc(ext->n, int);
    for (i = 0; i < ext->n; i++) ext->y[i] = y[i];
/**************** Filled in ext  *************************/

    ki = ext->ki;
    ki_tmp = Calloc(ext->n, int);

    for (i = 0; i < ext->n; i++){
	ki[i] = i;
    }

    
/* Note that this searches for a minimum: (!!) */
/*
    vmax = vmaxget();

    vmmin(*p, b, &Fmin,
	  bfun, bfun_gr, *maxit, *trace,
	  mask, abstol, reltol, nREPORT,
	  ext, &fncount, &grcount, &fail);
    *convergence = (fail == 0);
    vmaxset(vmax);
    bfun_gr(*p, b, gr, ext); 
*/
    Fmin = bfun(ext->p, b, ext);

    *loglik = -Fmin;

    for (i = 0; i < ext->n_fam; i++)
	frail[i] = ext->gamma[i];

/* Gone for now; predicted comes from calling R function
    if (ext->family <= 1){ 
	for (j = 0; j < ext->n; j++)
	    predicted[j] = ext->pred[j];
*/

    upper = 0;

/************** Bootstrapping starts *****************************/

    for (i = 0; i < *boot; i++){
	if ((i / 10) * 10 == i)
	    if (*trace)
		printf("****************************** Replicate No. %d\n", i);
	if (*conditional){
	    permute(ext->n, ki, ki_tmp);
	    for (j = 0; j < ext->n; j++){
		ext->y[j] = y[ki[j]];
		ext->offset[j] = offset[ki[j]];
		ext->cluster[j] = cluster[ki[j]];
	    }
	}else{
	    if (*family <= 1){
		for (j = 0; j < ext->n; j++)
		    ext->y[j] = rbinom(1, predicted[j]);
	    }else{
		for (j = 0; j < ext->n; j++)
		    ext->y[j] = rpois(predicted[j]);
	    }		
	}

	Fmin = bfun(ext->p, b, ext);
	boot_log[i] = -Fmin;
	if (-Fmin >= *loglik) upper++;
    }

    if (*boot) *boot_p = (double)upper / (double)*boot;
    else *boot_p = 1.0;

    PutRNGstate();

/*    vmaxset(vmax1); */

    Free(ext->y);
    Free(ext->gamma);
    Free(ext->cluster);
    Free(ext->ki);
    Free(ext->offset);
    Free(ext->pred);
    Free(ext->x_beta);
    Free(ext->fam_out);
    Free(ext->success);

    Free(ext);

    Free(ki_tmp);
}
