#include <stdio.h>
#include <math.h>

#include "glmmboot.h"
#include "fun.h"
#include "bfun.h"
#include <Rmath.h>
#include "GB_zeroin.h"
#include <R_ext/Linpack.h>

extern P_fun *P;
extern G_fun *G;
extern H_fun *H;

static double get2_gam(int n, /* family size */
		       double *weights,
		       double *lin, /* x %*% beta for this family */ 
		       double *yw,
		       double gam_start);

static double get3_gam(int n, /* family size */
		       double *weights,
		       double *lin, /* x %*% beta for this family */ 
		       double ytot);
	
static double gam_fun(double gam, void *extra){
    int i;
    double res;
    Ext_gam *ex;

    ex = extra;

    /* res = ex->ytot; */
    res = 0.0;
    for (i = 0; i < ex->n; i++){
	/* res -= P(ex->lin[i] + gam, 1.0, 1.0);  check !!!!! */
	res += G(ex->lin[i] + gam, ex->yw[i], ex->weights[i]);
    }
    return(res);
}

static double get_gam(int n, /* family size */
	       double *lin, /* x %*% beta for this family */ 
	       int ytot){
    int i, iter;

    double gam, delta, eps;
    double *p;
    double tal, nam, upper, lower;

    p = Calloc(n, double);
    eps = 0.000001;

    if ((ytot == 0)) error("ytot must be strictly positive.");
    if ((ytot == n)) error("ytot must be strictly less than n.");

    upper = 0.0;
    lower = 0.0;
    for (i = 0; i < n; i++){
	if (lower < lin[i]) lower = lin[i];
	if (upper > lin[i]) upper = lin[i];
    }
    gam = log((double)ytot / (double)(n - ytot));
    lower = gam - lower;
    upper = gam - upper;
    gam = (upper + lower) * 0.5;
    if ((upper - lower) < eps) return(gam);

    delta = 1.0;
    iter = 0;
    while ((fabs(delta) > eps) && (iter < 100)){
	iter++;
	tal = ytot;
	nam = 0.0;
	for (i = 0; i < n; i++){
	    p[i] = P(lin[i] + gam, 1.0, 1.0); /* Change !!!!!! */
	    tal -= p[i];
	    nam += p[i] * (1.0 - p[i]);
	}
	delta = tal / nam;
	gam += delta;
    }
    tal = 0.0;
    for (i = 0; i < n; i++) tal += p[i];

    if (iter > 99){
	Rprintf("gamma = %f\n", gam);
	Rprintf("ytot = %d\n", ytot); 
	Rprintf("n = %d\n", n);
	Rprintf("upper = %f\n", upper);
	Rprintf("lower = %f\n", lower);
	for (i = 0; i < n; i++){
	    Rprintf("lin[%d] = %f\n", i, lin[i]);
	}
	error("Too many iterations in [get_gam]");
    }

    Free(p);
    return(gam);
}

static double get2_gam(int n, /* family size */
		       double *weights,
		       double *lin, /* x %*% beta for this family */ 
		       double *yw,
		       double gam_start){

/* For Binomial families */
    double ax;
    double bx;
    Ext_gam *extra;
    double Tol;
    int Maxit;
    double res;
    double gam;
    int i;

    extra = Calloc(1, Ext_gam);
    /*
    if (abs(ytot) < 0.0001) error("ytot must be strictly positive.");
    if (abs(ytot - wtot) < 0.0001) 
	error("ytot must be strictly less than wtot.");
    */
    ax = 0.0;
    bx = 0.0;
    for (i = 0; i < n; i++){
	if (ax < lin[i]) ax = lin[i];
	if (bx > lin[i]) bx = lin[i];
    }
    gam = gam_start;
    ax = gam - ax;
    bx = gam - bx; /* Check if not bx = gam + bx!? */

    extra->yw = yw;
    extra->weights = weights;
    extra->n = n;
    extra->lin = lin;

    Tol = 0.0;
    Maxit = 25;
    res = GB_zeroin(ax, bx, &gam_fun, extra, &Tol, &Maxit);
    Free(extra);
    return(res);
}

static double get3_gam(int n, /* family size */
		       double *weights,
		       double *lin, /* x %*% beta for this family */ 
		       double ytot){

/* For Poisson family */
    int j;
    double denom;

    denom = 0.0;
    for (j = 0; j < n; j++){
	denom += weights[j] * exp(lin[j]);
    }

    return ( log( ytot / denom ) );
}


double bfun(int p, double *b, void *ex){
    int i;
    int j;
    int q; /* No. of families */
    int indx;
    double ytot, wtot;
    double loglik;

    double *lin;

    int *cluster;
    int *ki;
    double **x;
    double *y;
    double *offset;

    Extb *ext;

    ext = ex;

    cluster = ext->cluster;
    ki = ext->ki;
    x = ext->x;
    y = ext->yw;
    offset = ext->offset;

    lin = ext->x_beta;

    q = ext->n_fam;

/* Get the "linear predictor": */

    if (p > 0){
	indx = -1;
	for (i = 0; i < ext->n; i++){
	    lin[i] = offset[i];
	    for (j = 0; j < p; j++){
		indx++;
		lin[i] += b[j] * x[i][j];
	    }
	}
    }else{ /* Null model */
	for (i = 0; i < ext->n; i++) lin[i] = offset[i];
    }

/* Now get the gamma's: */
    
    indx = 0;
    if (ext->family <= 1){ /* binomial family */
	for (i = 0; i < ext->n_fam; i++){  /* NOT Excluding first family!! */
	    ytot = 0.0;
	    wtot = 0.0;
	    for (j = 0; j < ext->fam_size[i]; j++){
		ytot += ext->yw[indx + j];
		wtot += ext->weights[indx + j];
	    }
	    if (abs(ytot) < 0.001){
		ext->fam_out[i] = -1;
		ext->gamma[i] = -1000; /* -inf */
	    }else if (abs(ytot - wtot) < 0.001){
		ext->fam_out[i] = 1;
		ext->gamma[i] = 1000;  /* +inf */
	    }else{
		ext->fam_out[i] = 0;
		ext->gamma[i] = get2_gam(ext->fam_size[i],
					 (ext->weights + indx), 
					 (lin + indx),
					 (ext->yw + indx),
					 log( ytot / (wtot - ytot) ));
	    }
	    indx += ext->fam_size[i];
	}
    }else{ /* Poisson; ext->family == 2 */
	for (i = 1; i < ext->n_fam; i++){ /* Excluding first family!! */
	    ytot = 0.0;
	    wtot = 0.0;
	    for (j = 0; j < ext->fam_size[i]; j++){
		ytot += ext->yw[indx + j];
		wtot += ext->weights[indx + j];
	    }
	    if (abs(ytot) < 0.001){
		ext->fam_out[i] = -1;
		ext->gamma[i] = -1000; /* -inf */
	    }else{
		ext->fam_out[i] = 0;
		ext->gamma[i] = get3_gam(ext->fam_size[i],
					 ext->weights + indx,
					 (lin + indx),
					 ytot);
	    }
	    indx += ext->fam_size[i];
	}
    }

/* Now get the log likelihood: */
    loglik = 0.0;
    indx = -1;
    /* printf("beta[%d] = %f\n", 0, b[0]);  */
    for (i = 0; i < ext->n_fam; i++){
	if (ext->fam_out[i] == 0){
	    for (j = 0; j < ext->fam_size[i]; j++){
		indx++;
       
		loglik += log( P(lin[indx] + ext->gamma[i], 
				 y[indx], ext->weights[indx]) );
	    }
	}else{
	    indx += ext->fam_size[i];
	}
    }

    return(-loglik); /* Return: -loglikelihood */
    }

void bfun_gr(int n, double *b, double *gr, void *ex){
    int i;
    int j;
    int s;
    int indx;
    double tmp;

    double *pe, *lin;
    double **x;

    Extb *ext;

    ext = ex;
    pe = ext->pred;
    x = ext->x;

    lin = ext->x_beta;

/* Calculate the 'predicted values' at (beta, gamma): */
    indx = -1;

    if (ext->family <= 1){ /* Bernoulli case */
	for (i = 0; i < ext->n_fam; i++){
	    if (ext->fam_out[i] == 0){
		for (j = 0; j < ext->fam_size[i]; j++){
		    indx++;
		    ext->pred[indx] = P(ext->gamma[i] + lin[indx], 1.0, 
					ext->weights[indx]); /* ????? */
		}
	    }else{
		tmp = (double)(ext->fam_out[i] == 1);
		for(j = 0; j < ext->fam_size[i]; j++){
		    indx++;
		    ext->pred[indx] = tmp;
		}
	    }
	}
	 
    }else{ /* Poisson case */
	for (i = 0; i < ext->n_fam; i++){
	    if (ext->fam_out[i] == 0){
		for (j = 0; j < ext->fam_size[i]; j++){
		    indx++;
		    ext->pred[indx] = exp(ext->gamma[i] + lin[indx]);
		}
	    }else{
		for(j = 0; j < ext->fam_size[i]; j++){
		    indx++;
		    ext->pred[indx] = 0.0;
		}
	    }
	}
    }

    for (s = 0; s < ext->p; s++){
	gr[s] = 0.0;
	indx = -1;
	for (i = 0; i < ext->n_fam; i++){
	    if (ext->fam_out[i]) indx += ext->fam_size[i];
	    else{
/* Calculate gr[s]: */
		for (j = 0; j < ext->fam_size[i]; j++){
		    indx++;
		    gr[s] += (ext->yw[indx] - pe[indx]) *
			x[indx][s];
		}
	    }
	}
    }


    /* Return minus the gradient! */
    for (s = 0; s < n; s++) gr[s] = -gr[s];
}

void bfun_hess(int p, double *b, double *hessian, Extb *ext){

    int m, s, i, j, indx;
    
    double *h, *h_fam;
    double **hess;
    double gam, t1, t2;
    
    h = Calloc(ext->n, double);
    h_fam = Calloc(ext->n_fam, double);
    hess = Calloc(p, double *);
    for (m = 0; m < p; m++){
	hess[m] = hessian + m * p;
    }

    for (i = 0; i < ext->n; i++) h[i] = 0.0;

    indx = -1;
    for (i = 0; i < ext->n_fam; i++){
	h_fam[i] = 0.0;
	if (ext->fam_out[i] == 0){
	    gam = ext->gamma[i];
	    for (j = 0; j < ext->fam_size[i]; j++){
		indx++;
		h[indx] = H(ext->x_beta[indx] + gam, ext->yw[indx],
				ext->weights[indx]);
		h_fam[i] += h[indx];
	    }
	}else{
	    indx += ext->fam_size[i];
	}
    }
    
    for (m = 0; m < p; m++){
	for (s = 0; s <= m; s++){
	    hess[m][s] = 0.0;
	}
    }
    
    for (m = 0; m < p; m++){
	for (s = 0; s <= m; s++){
	    for (i = 0; i < ext->n; i++){
		hess[m][s] += ext->x[i][m] * ext->x[i][s] * h[i];
	    }
	    indx = -1;
	    for (i = 0; i < ext->n_fam; i++){
		if (ext->fam_out[i] == 0){
		    t1 = 0.0;
		    t2 = 0.0;
		    for (j = 0; j < ext->fam_size[i]; j++){
			indx++;
			t1 += ext->x[indx][m] * h[indx];
			t2 += ext->x[indx][s] * h[indx];
		    }
		    hess[m][s] -= t1 * t2 / h_fam[i];
		}else{
		    indx += ext->fam_size[i];
		}
	    }
	}
    }
    
    for (m = 0; m < p; m++){
	for (s = m + 1; s < p; s++){
	    hess[m][s] = hess[s][m];
	}
    }
    Free(h_fam);
    Free(h);
}
   
