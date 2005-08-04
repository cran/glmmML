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
extern Gprim_fun *Gprim;

double bfun(int p, double *b, void *ex){
    int i;
    int j;
    int q; /* No. of families */
    int indx;
    int ytot;
    double loglik;

    double *lin;

    int *cluster;
    int *ki;
    double **x;
    int *y;
    double *offset;

    Extb *ext;

    ext = ex;

    cluster = ext->cluster;
    ki = ext->ki;
    x = ext->x;
    y = ext->y;
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
	    ytot = 0;
	    for (j = 0; j < ext->fam_size[i]; j++)
		ytot += ext->y[indx + j];
	    if (ytot == 0){
		ext->fam_out[i] = -1;
		ext->gamma[i] = -1000; /* -inf */
	    }else if (ytot == ext->fam_size[i]){
		ext->fam_out[i] = 1;
		ext->gamma[i] = 1000;  /* +inf */
	    }else{
		ext->fam_out[i] = 0;
		ext->gamma[i] = get2_gam(ext->fam_size[i], 
					 (lin + indx),
					 ytot);
	    }
	    indx += ext->fam_size[i];
	}
    }else{ /* Poisson; ext->family == 2 */
	for (i = 1; i < ext->n_fam; i++){ /* Excluding first family!! */
	    ytot = 0;
	    for (j = 0; j < ext->fam_size[i]; j++)
		ytot += ext->y[indx + j];
	    if (ytot == 0){
		ext->fam_out[i] = -1;
		ext->gamma[i] = -1000; /* -inf */
	    }else{
		ext->fam_out[i] = 0;
		ext->gamma[i] = get3_gam(ext->fam_size[i],
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
		loglik += log( P(lin[indx] + ext->gamma[i], y[indx]) );
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
		    ext->pred[indx] = P(ext->gamma[i] + lin[indx], 1);
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
		    gr[s] += (ext->y[indx] - pe[indx]) *
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
    
    double *H, *H_fam;
    double **hess;
    double gam, t1, t2;
    
    H = Calloc(ext->n, double);
    H_fam = Calloc(ext->n_fam, double);
    hess = Calloc(p, double *);
    for (m = 0; m < p; m++){
	hess[m] = hessian + m * p;
    }

    for (i = 0; i < ext->n; i++) H[i] = 0.0;

    indx = -1;
    for (i = 0; i < ext->n_fam; i++){
	H_fam[i] = 0.0;
	if (ext->fam_out[i] == 0){
	    gam = ext->gamma[i];
	    for (j = 0; j < ext->fam_size[i]; j++){
		indx++;
		H[indx] = Gprim(ext->x_beta[indx] + gam, ext->y[indx]);
		H_fam[i] += H[indx];
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
		hess[m][s] += ext->x[i][m] * ext->x[i][s] * H[i];
	    }
	    indx = -1;
	    for (i = 0; i < ext->n_fam; i++){
		if (ext->fam_out[i] == 0){
		    t1 = 0.0;
		    t2 = 0.0;
		    for (j = 0; j < ext->fam_size[i]; j++){
			indx++;
			t1 += ext->x[indx][m] * H[indx];
			t2 += ext->x[indx][s] * H[indx];
		    }
		    hess[m][s] -= t1 * t2 / H_fam[i];
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
    Free(H_fam);
    Free(H);
}
   
	
double gam_fun(double gam, void *extra){
    int i;
    double res;
    Ext_gam *ex;

    ex = extra;

    res = ex->ytot;
    for (i = 0; i < ex->n; i++){
	res -= P(ex->lin[i] + gam, 1);
    }
    return(res);
}

double get_gam(int n, /* family size */
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
	    p[i] = P(lin[i] + gam, 1);
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

double get2_gam(int n, /* family size */
	       double *lin, /* x %*% beta for this family */ 
	       int ytot){
    double ax;
    double bx;
    Ext_gam *extra;
    double Tol;
    int Maxit;
    double res;
    double gam;
    int i;

    extra = Calloc(1, Ext_gam);
    if ((ytot == 0)) error("ytot must be strictly positive.");
    if ((ytot == n)) error("ytot must be strictly less than n.");

    ax = 0.0;
    bx = 0.0;
    for (i = 0; i < n; i++){
	if (ax < lin[i]) ax = lin[i];
	if (bx > lin[i]) bx = lin[i];
    }
    gam = log((double)ytot / (double)(n - ytot));
    ax = gam - ax;
    bx = gam - bx;

    extra->ytot = ytot;
    extra->n = n;
    extra->lin = lin;

    Tol = 0.0;
    Maxit = 25;
    res = GB_zeroin(ax, bx, &gam_fun, extra, &Tol, &Maxit);
    Free(extra);
    return(res);
}

double get3_gam(int n, /* family size */
	       double *lin, /* x %*% beta for this family */ 
	       int ytot){
    int j;
    double denom;

    denom = 0.0;
    for (j = 0; j < n; j++){
	denom += exp(lin[j]);
    }

    return ( log( (double)ytot / denom ) );
}



#include "fun.h"

void bnr_opt(int bdim, double *beta, double *loglik, int *mask, 
	     Extb *ext, double epsilon, int maxit, int trace){

    /* Start a Newton-Raphson thing: */
    
    double fstart, fprev;
    int iter;

    double *db = NULL;
    double *gr = NULL;
    double *hess = NULL;
    void *ex;

    int one = 1;
    int conver = 0;

    int info = 1;
    double L1;
    int i;

    int p, k, m;
    int true_bdim;

    int *ipiv;
    double *work;
    int lwork;

    double rcond;
    double *det;
    int job = 11;
    int pos_def;
    double *c_hess;

    c_hess = Calloc(bdim * bdim, double);
    det = Calloc(2, double);

    p = ext->p;

    true_bdim = p; /* p = bdim here!!!!!!!! */

    db = Calloc(bdim, double);
    ipiv = Calloc(bdim, int);
    lwork = 11 * bdim;
    work = Calloc(lwork, double);

    gr = ext->gr;
    /* hess = ext->hessian; */
    ex = ext;

    fun2(bdim, beta, loglik, gr, hess, ex);
    fstart = *loglik;
    fprev = fstart;
  
    for (iter = 0; iter < maxit; iter++){
	F77_CALL(dcopy)(&true_bdim, gr, &one, db, &one);

	pos_def = 0;
	while (!pos_def){
	    F77_CALL(dpoco)(hess, &bdim, &true_bdim, &rcond, work, &info);
	    
	    if (info){
		printf("Hessian not positive definite.\n");
		printf("info = %d\n", info);
		if (true_bdim == bdim){
		    fun2(bdim, beta, loglik, gr, hess, ex);
		    printf("We try fixing sigma at %f\n", exp(beta[bdim-1]));
		    true_bdim--;
		}else{
		    printf("sigma currently = %f", exp(beta[bdim -1]));
		    error("Try another start value for sigma.\n");
		}
		F77_CALL(dpoco)(hess, &bdim, &true_bdim, 
				&rcond, work, &info);
		if (info) error("Try another start value for sigma.\n");
		pos_def = 1;
	    }else{
		pos_def = 1;
	    }
	}
	F77_CALL(dposl)(hess, &bdim, &true_bdim, db);
	
	L1 = 0.0; /*F77_CALL(dnrm2)(&bdim, db, &one);*/
	for (i = 0; i < true_bdim; i++){
	    L1 += fabs(db[i]);
	    beta[i] += db[i];
	}
	if (trace)
	    printf("*** Iteration %d: L1 = %f, loglik = %f\n", iter, 
		   L1, *loglik);
	conver = (L1 < epsilon);
	if (!conver){ 
	    conver = 
		fabs(*loglik - fprev) < epsilon;/* * (fabs(fprev) + eps); */
	}

	if (conver){
	    if (iter > 0){ /* Always take at least one step! */
		if (trace){
		    printf("Newton-Raphson CONVERGENCE in %d step(s)!!\n", 
			   iter);
		}
		break;
	    }
	}
	fprev = *loglik;
	fun2(bdim, beta, loglik, gr, hess, ex);
	if (*loglik < fprev){
	    Rprintf("Warning: Decreasing loglik!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
    }

    /* Note: __packed__ form (should use pos. definite versions?) */
    /* F77_CALL(dpptrf)(&u, &(bdim), hess, &info); */

    for ( m = 0; m < bdim * bdim; m++){
	c_hess[m] = hess[m];
    }
    F77_CALL(dpoco)(hess, &bdim, &true_bdim, &rcond, work, &info);
    if (info == 0){
	F77_CALL(dpodi)(hess, &bdim, &true_bdim, det, &job);

	for (m = 0; m < bdim; m++){
	    for (k = 0; k < m; k++){
		hess[m + k * bdim] = hess[k + m * bdim];
	    }
	}
    }else{
	Rprintf("info = %d\n", info);
	for ( m = 0; m < bdim * bdim; m++){
	    hess[m] = c_hess[m];
	}
	warning("Hessian non-positive definite. No variance!");
    }

    Free(c_hess);
    Free(work);
    Free(ipiv);
    Free(db);
    Free(det);
}
