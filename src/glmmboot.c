#include <stdio.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

#include "glmmboot.h"
#include "bfun.h"
#include "fun.h"

extern P_fun *P;
extern G_fun *G;
extern H_fun *H;

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
	       int *p, 
	       double *start_beta,
	       int *cluster,
	       double *weights,
	       double *x, /* Now p x (\sum_1_{n_fam} fam_size[i]) */
	       double *y,
	       double *offset,
	       int *fam_size,
	       int *n_fam,
	       double *epsilon,
	       int *maxit,
	       int *trace,
	       int *boot,
	       double *beta,
	       double *predicted,
	       double *loglik,
	       double *variance,
	       int *info, 
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
    int bdim, m, k;

    double ** hessian;
    double *hess_vec;

    double *det;
    int lwork;
    double *work;
    double rcond;
    int job = 11;

    bdim = *p;
    lwork = 11 * (*p);
    work = Calloc(lwork, double);
    det = Calloc(2, double);

    hessian = Calloc(bdim, double *);
    hess_vec = Calloc(bdim * bdim, double);
    for (j = 0; j < bdim; j++) hessian[j] = hess_vec + j * bdim;

    GetRNGstate(); /* For random number generation */

/*    vmax1 = vmaxget(); */

    if (*family == 0){
	P = &P_logit;
	G = &G_logit;
	H = &H_logit;
    }else if (*family == 1){
	P = &P_cloglog;
	G = &G_cloglog;
	H = &H_cloglog;
    }else if (*family == 2){
	P = &P_poisson;
	G = &G_poisson;
	H = &H_poisson;
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
    ext->yw = Calloc(ext->n, double);
    for (i = 0; i < ext->n; i++) 
	ext->yw[i] = y[i] * weights[i];
    ext->weights = weights;
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

/* Done in calling R function.... 
    if (*family <= 1){
	for (j = 0; j < ext->n; j++)
	    predicted[j] = P(ext->x_beta[j], 1);

    }else{
	for (j = 0; j < ext->n; j++)
	    predicted[j] = exp(ext->x_beta[j]);
    }
*/
    bfun_hess(*p, beta, hess_vec, ext);


    Rprintf("Hessian...\n\n");
    for (i = 0; i < *p; i++){
	for (j = 0; j < *p; j++){
	    Rprintf("%f  ", hessian[i][j]);
	}
	Rprintf("\n");
    }


    F77_CALL(dpoco)(*hessian, &bdim, &bdim, &rcond, work, info);
    if (*info == 0){
	F77_CALL(dpodi)(*hessian, p, p, det, &job);
	for (m = 0; m < bdim; m++){
	    for (k = 0; k < m; k++){
		hessian[k][m] = hessian[m][k];
	    }
	}
	for (m = 0; m < bdim * bdim; m++) variance[m] = hess_vec[m];
	
    }else{
	Rprintf("info[dpoco] = %d\n", *info);
	warning("[glmmboot:] Information non-positive definite. No variance!");
    }

    upper = 0;

/************** Bootstrapping starts *****************************/
 
    for (i = 0; i < *boot; i++){
	if (*trace){
	    if ((i / 10) * 10 == i)
		printf("********************* Replicate No. No. %d\n", i);
	}
	if (*family <= 1){ /* Bernoulli */
	    for (j = 0; j < ext->n; j++)
		ext->yw[j] = rbinom((int)weights[j], predicted[j]);
	}else{
	    for (j = 0; j < ext->n; j++) /* Poisson */
		ext->yw[j] = rpois(weights[j] * predicted[j]);
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


    Free(hessian);
    Free(hess_vec);
    Free(det);
    Free(work);

    Free(ext->yw);
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
