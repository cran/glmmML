#include <stdio.h>
#include <R_ext/Applic.h>

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
	       double *hessian,
	       double *frail,
	       double *boot_p,
	       double *boot_log,
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
    Extb *ext;
    int i;
    int j;

    int *ki;
    int *ki_tmp;

    double Fmin;
    double *b;
    double *gr;
    int upper;
    char *vmax;

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
    ext->p = *p;
    ext->fam_size = fam_size;
    ext->n_fam = *n_fam;
    ext->success = Calloc(*n_fam, int);
    ext->fam_out = Calloc(ext->n_fam, int);

    ext->x = Calloc(ext->n, double *);
    for (i = 0; i < ext->n; i++){
	ext->x[i] = x + i * (ext->p);
    }
    /*** Note that ext->x is not "filled"; ***/ 
    /*** only points to the right place    ***/

    ext->x_beta = Calloc(ext->n, double);
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
    ext->gr = Calloc(ext->p, double);
    ext->hessian = Calloc(ext->p * ext->p, double);
    ext->y = Calloc(ext->n, int);
    for (i = 0; i < ext->n; i++) ext->y[i] = y[i];
/**************** Filled in ext  *************************/
    mask = Calloc(ext->p, int);    

    b = Calloc(ext->p, double);

    gr = ext->gr;

    for (i = 0; i < *p; i++){
	b[i] = start_beta[i];
    }

    for (i = 0; i < ext->p; i++){
        mask[i] = 1;
    }

    ki = ext->ki;
    ki_tmp = Calloc(ext->n, int);

    for (i = 0; i < ext->n; i++){
	ki[i] = i;
    }

/* Note that this searches for a minimum: (!!) */
    vmax = vmaxget();

    vmmin(*p, b, &Fmin,
	  bfun, bfun_gr, *maxit, *trace,
	  mask, abstol, reltol, nREPORT,
	  ext, &fncount, &grcount, &fail);
    *convergence = (fail == 0);
    vmaxset(vmax);
    bfun_gr(*p, b, gr, ext);
    if(*trace){
	Rprintf("Max log likelihood after vmmin: %f\n", -Fmin);
	printf("Gradients: ");
	for (i = 0; i < *p; i++){
	    Rprintf(" %f, ", -ext->gr[i]);
	}
	Rprintf("\n");
    }
    
    *loglik = -Fmin;
    for (i = 0; i < *p; i++){
	beta[i] = b[i];
    }
    for (i = 0; i < ext->n_fam; i++){
	frail[i] = ext->gamma[i];
    }

    bfun_hess(*p, beta, hessian, ext);

    upper = 0;

/************** Bootstrapping starts *****************************/

    for (i = 0; i < *boot; i++){
	if (*trace){
	    if ((i / 10) * 10 == i)
		printf("********************* Replicate No. No. %d\n", i);
	}
	permute(ext->n, ki, ki_tmp);
	for (j = 0; j < ext->n; j++){
	    ext->y[j] = y[ki[j]];
	    ext->x[j] = x + ki[j] * (ext->p);
	    ext->offset[j] = offset[ki[j]];
	    ext->cluster[j] = cluster[ki[j]];
	}
/* Restore beta as start values: */
	for ( j = 0; j < *p; j++) b[j] = beta[j];
	
	vmax = vmaxget();
	vmmin(*p, b, &Fmin,
	      bfun, bfun_gr, *maxit, *trace,
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

/*    vmaxset(vmax1); */

    Free(ext->y);
    Free(ext->hessian);
    Free(ext->gr);
    Free(ext->gamma);
    Free(ext->cluster);
    Free(ext->fam_out);
    Free(ext->ki);
    Free(ext->offset);
    Free(ext->x_beta);
    Free(ext->pred);
    Free(ext->x);
    Free(ext->success);

    Free(ext);

    Free(mask);
    Free(ki_tmp);
    Free(b);

}
