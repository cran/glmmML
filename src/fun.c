#include <stdio.h>
#include <math.h>

#include "glmmml.h"
#include "fun.h"
#include <Rmath.h>
#include <R_ext/Linpack.h>

extern P_fun *P;
extern G_fun *G;
extern Gprim_fun *Gprim;

/***********************************************************/
/*         Bernoulli distribution, logit link:             */

double P_logit(double x, int y){ /* logit link */
    
    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;
    double res;

    res = plogis(x, location, scale, y, give_log);
    return ( res );
}

double G_logit(double x, int y){
    
    /* Calculates G = P'/P */
    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;
    double res;

    if (y) {
	res = plogis(x, location, scale, 0, give_log);
    }else{
	res = -plogis(x, location, scale, 1, give_log);
    }
    return ( res );
}

double Gprim_logit(double x, int y){ /* Note: Independent of y */
 
    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;

    return ( -dlogis(x, location, scale, give_log) );
}

/******************************************************************/
/*         Bernoulli distribution, cloglog link:                  */

double P_cloglog(double x, int y){
    double q;

    q = exp(-exp(x));

    if (y) {
	return ( 1.0 - q );
    }else{
	return ( q );
    }
}

double G_cloglog(double x, int y){

    double exp_x;

    if (y){
	if (x < -6.0){
	    return (1.0);
	}else if (x > 6.0){
	    return (0.0);
	}else{
	    exp_x = exp(x);
	    return ( exp_x / expm1(exp_x) );
	}
    }else{
	return(-exp(x));
    }
}

double Gprim_cloglog(double x, int y){
    double exp_x;
    double rat;

    if (y){
	if (x < -6.0){
	    return (0.0);
	}else if (x > 6.0){
	    return (0.0);
	}else{
	    exp_x = exp(x);
/*
	    return ( exp_x * (expm1(exp_x) - exp_x * exp(exp_x)) / 
		    expm1(exp_x) / expm1(exp_x) );
*/
	    rat = -exp_x * exp(-exp_x) / expm1(-exp_x);
	    return (-rat * (expm1(x) + rat) );
	}
    }else{
	return(-exp(x));
    }
}

/*****************************************************************/
/*       Poisson distribution, log link:                         */

double P_poisson(double x, int y){

    int give_log = 0;
    return ( dpois(y, exp(x), give_log) );
}

double G_poisson(double x, int y){

    return (y - exp(x));
}

double Gprim_poisson(double x, int y){ /* Note: Independent of y */

    return (-exp(x));
}

/*****************************************************************/
/*           Normal distribution, identity link:                 */
/* Needs modification!!! ('int y' no good!)                      */

/*********************
double P_normal(double x, int y){

}

double G_normal(double x, int y){

}

double Gprim_normal(double x, int y){

}
********************/
static void update(int level,
		   int p, 
		   double *beta,
		   double *loglik, 
		   double *score,
		   double *hessian,
		   int n,
		   double *x,
		   double *x_beta,
		   int *y,
		   Exts *ext){
    
    double h;
    double *hb = NULL;
    double *hbb = NULL; /* (p+1) x (p+1) */
    double *pip = NULL; /* ext->n_points */
    double *xG = NULL;  /* ext->n_points * (p+1) */
    double sigma;
    double tmp, tmp2;
    int i, j, m, k;
    double factor;
    int count;
    
    if (level < 0) return;

    factor = 1.0;

    sigma = exp(beta[p]);

    pip = Calloc(ext->n_points, double);
    if (level > 0){
	xG = Calloc((p + 1) * ext->n_points, double);
	hb = Calloc((p + 1), double);
	hbb = Calloc((p + 1) * (p + 1), double);
    }


/**********************************************************************/    
/* Calculate  h  for the loglik, and pip ("weights"): */
    h = 0.0;
    for (i = 0; i < ext->n_points; i++){
	tmp = 1.0;
	for (j = 0; j < n; j++){
	    tmp *= P(x_beta[j] + ext->zeros[i] * sigma, y[j]);
	}
	pip[i] = tmp * ext->weights[i];
	h += pip[i];
    }

    /* Add into the loglik; We allow *loglik = -Inf!! */
    *loglik += log(h);

    if (level == 0) {
	Free(pip);
	return;
    }

    if (h <= 0.0) { /* We don't allow h == 0.0 in the following! */
	printf("h = 0.0; trying to fix...\n");
	factor = 2.0 * factor;
	h = 0.0;
	count = 0;
	while ( (h <= 0.0) && (count < 10) ){
	    count++;
	    for (i = 0; i < ext->n_points; i++){
		tmp = 1.0;
		for (j = 0; j < n; j++){
		    tmp *= factor * P(x_beta[j] + ext->zeros[i] * sigma, y[j]);
		}
		pip[i] = tmp * ext->weights[i];
		h += pip[i];
	    }
	}
	if (h <= 0){
	    error("Unable do get likelihood function POSITIVE!!!!!!!!!\n");
	}
    }


/* Fill xG: */
    for (m = 0; m < p; m++){ /* First without w = log(sigma) */
	for (i = 0; i < ext->n_points; i++){
	    tmp = 0.0;
	    for (j = 0; j < n; j++){
/*
		tmp2 = G(x_beta[j] + sigma * ext->zeros[i], y[j]);
	        printf("G = %f,  X = %f, Y = %d\n", 
		       tmp2, x_beta[j] + sigma * ext->zeros[i], y[j]);
*/
		tmp += factor * x[m + j * p] * 
		    G(x_beta[j] + sigma * ext->zeros[i], y[j]);
	    }
	    xG[i + m * ext->n_points] = tmp;
	}
    }
    for (i = 0; i < ext->n_points; i++){ /* Then for w = log(sigma) */
	tmp = 0.0;
	for (j = 0; j < n; j++){
	    tmp += factor * sigma * ext->zeros[i] * 
		G(x_beta[j] + sigma * ext->zeros[i], y[j]);
	}
	xG[i + p * ext->n_points] = tmp;
    }
    /* Done with xG */

/***********************************************************************/    
/* First derivatives, hb[]. */
    for (m = 0; m <= p; m++){ /* hb[m]; note w = log(sigma) INCLUDED! */
	tmp = 0.0;
	for (i = 0; i < ext->n_points; i++){
	    tmp += pip[i] * xG[i + m * ext->n_points];
	}
	hb[m] = tmp;
    }

    /* Add into first derivatives: */
    for (m = 0; m <= p; m++){
	score[m] += hb[m] / h;
    }

    if (level == 1){
	Free(xG);
	Free(pip);
	Free(hbb);
	Free(hb);
	return;
    }
	
/***********************************************************************/
/* Done with first derivatives. On to the hessian: */

    /* First the pxp matrix of 'coefficients': */
    for (m = 0; m < p; m++){
	for (k = 0; k <= m; k++){
	    tmp = 0.0;
	    for (i = 0; i < ext->n_points; i++){
		tmp += pip[i] * xG[i + m * ext->n_points] *
		    xG[i + k * ext->n_points];
	    }
	    hbb[m + k * (p + 1)] = tmp;
	    tmp = 0.0;
	    for (i = 0; i < ext->n_points; i++){
		tmp2 = 0.0;
		for (j = 0; j < n; j++){
		    tmp2 += factor * x[m + j * p] * x[k + j * p] *
			Gprim(x_beta[j] + sigma * ext->zeros[i], y[j]);
		}
		tmp += tmp2 * pip[i];
	    } /* And add in: */
	    hbb[m + k * (p + 1)] += tmp;
	}
    }

    /* Then the "beta[m] with beta[p] = log(sigma), m = 0,...,(p-1) */
    m = p;
    for (k = 0; k < m; k++){ /* Changed k <= m --> k < m */ 
	tmp = 0.0;
	for (i = 0; i < ext->n_points; i++){
	    tmp += pip[i] * xG[i + m * ext->n_points] *
		xG[i + k * ext->n_points];
	}
	hbb[m + k * (p + 1)] = tmp;
	tmp = 0.0; 
	for (i = 0; i < ext->n_points; i++){
	    tmp2 = 0.0;
	    for (j = 0; j < n; j++){ /* Here is the difference: */
		tmp2 += factor * x[k + j * p] * sigma * ext->zeros[i] * 
		    /* Changed n*m --> n*k above! */
		    Gprim(x_beta[j] + sigma * ext->zeros[i], y[j]);
	    }
	    tmp += tmp2 * pip[i];
	} /* And add in: */
	hbb[m + k * (p + 1)] += tmp;
/*
	if (k < p)
	    printf("hbb[%d, %d] = %f\n", m, k, hbb[m+k*(p+1)]);
*/
    }

    /* And finally beta[p] with beta[p] (AKA log(sigma) = w): */
    k = p;
    tmp = 0.0;
    for (i = 0; i < ext->n_points; i++){ /* Here is a difference: */
	tmp += pip[i] * xG[i + m * ext->n_points] *
	    (1.0 + xG[i + k * ext->n_points]);
    }
    hbb[m + k * (p + 1)] = tmp;
    tmp = 0.0; 
    for (i = 0; i < ext->n_points; i++){
	tmp2 = 0.0;
	for (j = 0; j < n; j++){ /* Here is the difference: */
	    tmp2 += factor * 
		(sigma * ext->zeros[i]) * (sigma * ext->zeros[i]) *
		Gprim(x_beta[j] + sigma * ext->zeros[i], y[j]);
	}
	tmp += tmp2 * pip[i];
    } /* And add in: */

    hbb[m + k * (p + 1)] += tmp;

   /* Now, add it into the hessian (lower triangle): */
    for (m = 0; m <= p; m++){
	for (k = 0; k <= m; k++){
	    hessian[m + k * (p+1)] += hbb[m + k * (p+1)] / h - 
		(hb[m] / h) * (hb[k] / h);
	    /*
	      if ((m == p) && (k == (p-1)))
		printf("hbb[%d, %d] = %f\n", m, k, hbb[m + k * (p+1)]);
	    */
	}
    }
    /* Fill in the upper triangle (symmetry): Necessary?? */
    for (m = 0; m <= p; m++){
	for (k = (m + 1); k <= p; k++){
	    hessian[m + k * (p+1)] = hessian[k + m * (p+1)];
	}
    }

/*********************************************************************/
/* We are done! Clean up 'the mess'! */

    Free(xG);
    Free(pip);
    Free(hbb);
    Free(hb);
}

static double frail_mean(int level,
		       int p, 
		       double *beta,
		       double *loglik, 
		       double *score,
		       double *hessian,
		       int n,
		       double *x,
		       double *x_beta,
		       int *y,
		       Exts *ext){
    
    double h, h_mean;
    double sigma;
    double tmp;
    int i, j;
    
    sigma = exp(beta[p]);


/**********************************************************************/    
/* Calculate  h  for the loglik, and pip ("weights"): */
    h = 0.0;
    h_mean =0.0;
    for (i = 0; i < ext->n_points; i++){
	tmp = 1.0;
	for (j = 0; j < n; j++){
	    tmp *= P(x_beta[j] + ext->zeros[i] * sigma, y[j]);
	}
	h += tmp * ext->weights[i];
	h_mean += tmp * ext->weights[i] * ext->zeros[i];
    }

    return ( h_mean / h);
}

void frail_fun(int pp1, 
		 double *beta,
		 double *frail,
		 void *ex){

    int start;
    double *gr = NULL;

    double tmp;
    int i, j;

    Exts *ext;
    double loglik;
    double *x;
    double *hessian = NULL;

    int level = 0;

    ext = ex;

    loglik = 0.0;

    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j + i * ext->p];
	}
	ext->x_beta[i] = tmp;
    }

/*
    F77_CALL(dgemm)(&trans, &trans, &(ex->n), 
    &p, &one, &alpha, x, &n, beta, &p,
		    &alpha, x_beta, &p);

    w = beta[ext->p];
    sigma = exp(w);
*/
    start = 0;
    for (i = 0; i < ext->n_fam; i++){
	x = ext->x + start * ext->p;
	frail[i] = frail_mean(level,
			      ext->p,
			      beta,
			      &loglik,
			      gr,
			      hessian,
			      ext->fam_size[i],
			      x,
			      ext->x_beta + start,
			      ext->y + start,
			      ext);
	start += ext->fam_size[i];
    }
}

void mu_fun(int bdim, double *b, double *mu, void *ex){

    int i;
    Exts *ext;

    ext = ex;

    for (i = 0; i < ext->n; i++){
	mu[i] = 0.0; /* More to come later */
    }
}

double fun(int pp1, 
	   double *beta, 
	   void *ex){

/* Dimensions:
   +++++++++++
   beta[p + 1] (0, ... p-1 = regression coefficients; p = log(sigma) = w
   x[p][n]
   y[n]
   fam_size[n_fam] { sum(fam_size) == n! }
   points[n_points][2] : first col: abscissas, second: weights.

   It is assumed that equal values in  'id'  comes together in groups.
*/

    int start;
    double *gr = NULL;

/*
    char trans = 'N';
    double alpha = 1.0;
    int one = 1;
*/
    double tmp;
    int i, j;

    Exts *ext;
    double loglik;
    double *x;
    double *hessian = NULL;

    int level = 0;

    ext = ex;

    loglik = 0.0;

    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j + i * ext->p];
	}
	ext->x_beta[i] = tmp;
    }

/*
    F77_CALL(dgemm)(&trans, &trans, &(ex->n), 
    &p, &one, &alpha, x, &n, beta, &p,
		    &alpha, x_beta, &p);

    w = beta[ext->p];
    sigma = exp(w);
*/
    start = 0;
    for (i = 0; i < ext->n_fam; i++){
	x = ext->x + start * ext->p;
	update(level,
	       ext->p,
	       beta,
	       &loglik,
	       gr,
	       hessian,
	       ext->fam_size[i],
	       x,
	       ext->x_beta + start,
	       ext->y + start,
	       ext);
	start += ext->fam_size[i];
    }
/*    printf("[fun]; loglik = %f\n", loglik); */
    return ( -loglik ); /* Note: minimizing!!! */
}


void fun1(int pp1, 
	  double *beta,
	  double *gr,
	  void *ex){

    int i, j, k;
    int start;

    double loglik;
    double *hessian = NULL;
    double *x;
    Exts *ext;
    double tmp;

    int level = 1;

/* Note that we here trust 'ext-x_beta' to be properly updated !! */
/* In 'fun' */
 
    ext = ex;

    loglik = 0.0;

    for (k = 0; k < pp1; k++){
	gr[k] = 0.0;
    }

    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j + i * ext->p];
	}
	ext->x_beta[i] = tmp;
    }
    
    start = 0;
    for (i = 0; i < ext->n_fam; i++){
	x = ext->x + start * ext->p;
	update(level,
	       ext->p,
	       beta,
	       &loglik,
	       gr,
	       hessian,
	       ext->fam_size[i],
	       x,
	       ext->x_beta + start,
	       ext->y + start,
	       ext);
	start += ext->fam_size[i];
    }
/*    printf("Gradient: "); */
    for (i = 0; i < pp1; i++){
	/* printf(" %f, ", gr[i]); */ 
	gr[i] = -gr[i]; /* Minimization! */
	ext->gr[i] = gr[i];
    }
    /* printf("\n"); */
}

void fun2(int pp1, 
	  double *beta,
	  double *loglik,
	  double *gr,
	  double *hessian,
	  void *ex){

    int i, j, k;
    int start;
    double tmp;
    double *x;
    Exts *ext;

    int level = 2;

    ext = ex;

    *loglik = 0.0;

    for (k = 0; k < pp1; k++){
	gr[k] = 0.0;
    }
    
    for (i = 0; i < pp1 * pp1; i++){
	hessian[i] = 0.0;
    }

    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j + i * ext->p];
	}
	ext->x_beta[i] = tmp;
    }
    
    start = 0;
    for (i = 0; i < ext->n_fam; i++){
	x = ext->x + start * ext->p;
	update(level,
	       ext->p,
	       beta,
	       loglik,
	       gr,
	       hessian,
	       ext->fam_size[i],
	       x,
	       ext->x_beta + start,
	       ext->y + start,
	       ext);
/*
	Free(x);
*/
        start += ext->fam_size[i];
	}
/* } */


   for (i = 0; i < pp1 * pp1; i++){
	hessian[i] = -hessian[i];
    }

}


void nr_opt(int bdim, double *beta, double *loglik, int *mask, 
	    Exts *ext, double epsilon, int maxit, int trace){
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

    true_bdim = 0;
    for (i = 0; i < bdim; i++){
	true_bdim += mask[i];
    }
    if ( (true_bdim < (bdim - 1)) || (true_bdim > bdim) ) 
	error("Error in [nr_opt]: true dimension wrong.");

    db = Calloc(bdim, double);
    ipiv = Calloc(bdim, int);
    lwork = 11 * bdim;
    work = Calloc(lwork, double);

    gr = ext->gr;
    hess = ext->hessian;
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
