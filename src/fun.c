#include <stdio.h>
#include <math.h>

#include "glmmml.h"
#include "fun.h"
#include <Rmath.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

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

double Hprim_logit(double x, int y){

    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;
    int lower_tail = 0;

    /* exp(x) * (exp(x) - 1.0) / (1 + exp(x))^3 */

    return(dlogis(x, location, scale, give_log) * 
	   expm1(x) * plogis(x, location, scale, lower_tail, give_log));
    
}

double Hbis_logit(double x, int y){

    return((4.0 * exp(2.0 * x) - exp(3.0 * x) - exp(x)) /
	   exp(4.0 * log(1 + exp(x))));

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

double Hprim_cloglog(double x, int y){
    double px, p2x, p3x, expx, exp2x;
    double res;

    expx = exp(x);
    res = -expx;

    if (y){
	px = 1.0 - exp(-exp(x));
	p2x = px * px;
	p3x = p2x * px;
	exp2x = expx * expx;
	res += expx * (p2x * (exp2x + 3.0 * expx + 1.0) -
		       3.0 * px * expx * (1.0 + expx) -
		       2.0 * exp2x) / p3x;
    }
    return(res);
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

double Hprim_poisson(double x, int y){

    return(-exp(x));
}

double Hbis_poisson(double x, int y){

    return(-exp(x));
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

static double fun_fam(int one, double *u, void *ex){
/* Calculates minus the log of the 'family' likelihood */
    double res;
    int j;
    double zero = 0.0;
    double ett = 1.0;
    int give_log = 1; /* NOTE! */

    Family *fam;

    fam = ex;

    res = dnorm(*u, zero, ett, give_log);

    for (j = 0; j < fam->n; j++){
/* This should be fixed via a parameter 'give_log' to P! */
 	    res += log(P(fam->x_beta[j] + *u * fam->sigma, fam->y[j]));
    }

    return(-res);
}

static void du_fun_fam(int one, double *u, double *gr, void *ex){
/* Calculates minus the derivative of the log of fun_fam at u */
    
    int j;
    double tmp;
    Family *fam;

    fam = ex;

    tmp = 0.0;
    for (j = 0; j < fam->n; j++) 
	tmp += G(fam->x_beta[j] + *u * fam->sigma, fam->y[j]);

    *gr = *u - fam->sigma * tmp;
}

static void du2_fun_fam(int one, double *u, double *grr, void *ex){
/* Calculates minus the second derivative of the log of fun_fam at u */
    
    int j;
    double tmp;
    Family *fam;

    fam = ex;

    tmp = 0.0;
    for (j = 0; j < fam->n; j++) 
	tmp += Gprim(fam->x_beta[j] + *u * fam->sigma, fam->y[j]);

    *grr = 1.0 - fam->sigma *fam->sigma * tmp;
}

static double fun1_fam(int one, double *u, void *ex){
    /* Calculates the first derivative map beta[m] of */
    /* minus the log of the 'family' likelihood       */
    double res, tmp;
    int j;
    double zero = 0.0;
    double ett = 1.0;
    int give_log = 1; /* NOTE! */
    int m;

    Family *fam;

    fam = ex;

    m = fam->m;

    res = dnorm(*u, zero, ett, give_log);

    for (j = 0; j < fam->n; j++){
/* This should be fixed via a parameter 'give_log' to P! */
	res += log(P(fam->x_beta[j] + *u * fam->sigma, fam->y[j]));
/*	    Rprintf("[fun_fam:] x_beta[%d] = %f\n", j, fam->x_beta[j]); */
    }

    tmp = 0.0;
    
    return(tmp);
/* This will be continued when we know whether it is worthwhile!!! */
/* Not used at the moment */
}

    
static void int_fun(double *x, 
		    int n, 
		    void *ex){
    /* Calculates one integrand (for one "family") in the expression     */
    /* of the log likelihood function.                                   */
    /* it is vectorizing, x is of length n. On input, x is the vector    */
    /* of points where int_fun should be evaluated, on output x contains */
    /* the n results.                                                    */

    /* NOTE: x[i] and fam->x[i] are different things here!! */

    double *res;
    int i, j;
    double zero = 0.0;
    double one = 1.0;
    int give_log = 0;

    Family *fam;

    fam = ex;

    res = Calloc(n, double);
    for (i = 0; i < n; i++) res[i] = dnorm(x[i], zero, one, give_log);

    for (i = 0; i < n; i++){
	for (j = 0; j < fam->n; j++){
 	    res[i] *= P(fam->x_beta[j] + x[i] * fam->sigma,fam->y[j]);
	}
	/* x[i] = res[i]; Improve by 'memcpy' later */
    }
    memcpy(x, res, n * sizeof(double));

    Free(res);
}

static void int_fun1(double *x, 
		     int n, 
		     void *ex){
    /* Calculates one integrand (for one "family") in the expression     */
    /* of the partial derivatives wrt                                    */
    /* 'beta' of the log likelihood function.                            */
    /* it is vectorizing, x is of length n. On input, x is the vector    */
    /* of points where int_fun should be evaluated, on output x contains */
    /* the n results.                                                    */
    /* NOTE: x[i] and fam->x[i] are different things here!! */

    double *res;
    int i, j;
    double zero = 0.0;
    double one = 1.0;
    int give_log = 0;

    double tsum;
    int m;

    Family *fam;

    fam = ex;

    m = fam->m;

    res = Calloc(n, double);
    for (i = 0; i < n; i++) res[i] = dnorm(x[i], zero, one, give_log);

    for (i = 0; i < n; i++){
	tsum = 0.0;
	for (j = 0; j < fam->n; j++){
	    tsum += fam->x[m][j] * 
		G(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	res[i] *= tsum;
    }
	
    for (i = 0; i < n; i++){
	for (j = 0; j < fam->n; j++){ 
	    res[i] *= P(fam->x_beta[j] + x[i] * fam->sigma,fam->y[j]);
	}
	/* x[i] = res[i]; Improve by 'memcpy' later */
    }
    memcpy(x, res, n * sizeof(double));

    Free(res);
}

static void int_fun1s(double *x, 
		      int n, 
		      void *ex){
    /* NOTE: x[i] and fam->x[i] are different things here!! */

    /* Calculates one integrand (for one "family") in the expression     */
    /* of the partial derivatives wrt                                    */
    /* 'log(sigma)' of the log likelihood function.                      */
    /* it is vectorizing, x is of length n. On input, x is the vector    */
    /* of points where int_fun should be evaluated, on output x contains */
    /* the n results.                                                    */
    double *res;
    int i, j;
    double zero = 0.0;
    double one = 1.0;
    int give_log = 0;

    double tsum;

    Family *fam;

    fam = ex;

    res = Calloc(n, double);
    for (i = 0; i < n; i++) res[i] = dnorm(x[i], zero, one, give_log);

    for (i = 0; i < n; i++){
	tsum = 0.0;
	for (j = 0; j < fam->n; j++){
	    tsum += G(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	res[i] *= tsum * fam->sigma * x[i];
	
    }
	
    for (i = 0; i < n; i++){
	for (j = 0; j < fam->n; j++){ 
	    res[i] *= P(fam->x_beta[j] + x[i] * fam->sigma,fam->y[j]);
	}
	/* x[i] = res[i]; Improve by 'memcpy' later */
    }
    memcpy(x, res, n * sizeof(double));

    Free(res);
}

static void int_fun2(double *x, 
		     int n, 
		     void *ex){
    /* Calculates one integrand (for one "family") in the expression     */
    /* of the partial second derivatives (hessian) wrt                   */
    /* 'beta[m], beta[k]' of the log likelihood function.                */
    /* it is vectorizing, x is of length n. On input, x is the vector    */
    /* of points where int_fun should be evaluated, on output x contains */
    /* the n results.                                                    */
    /* NOTE: x[i] and fam->x[i] are different things here!! */

    double *res;
    int i, j;
    double zero = 0.0;
    double one = 1.0;
    int give_log = 0;

    double *g1sum, *g2sum, *hsum;
    double tmp;
    int m, k;

    Family *fam;

    fam = ex;

    m = fam->m;
    k = fam->k;

    res = Calloc(n, double);
    g1sum = Calloc(n, double);
    g2sum = Calloc(n, double);
    hsum = Calloc(n, double);

    for (i = 0; i < n; i++) res[i] = dnorm(x[i], zero, one, give_log);

    /* First G-sum (m) */
    for (i = 0; i < n; i++){
	tmp = 0.0;
	for (j = 0; j < fam->n; j++){
	    tmp += fam->x[m][j] * 
		G(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	g1sum[i] = tmp;
    }
	
    /* Second G-sum (k) */
    for (i = 0; i < n; i++){
	tmp = 0.0;
	for (j = 0; j < fam->n; j++){
	    tmp += fam->x[k][j] * 
		G(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	g2sum[i] = tmp;
    }
	
    /* H-sum (m, k) */
    for (i = 0; i < n; i++){
	tmp = 0.0;
	for (j = 0; j < fam->n; j++){
	    tmp += 
		fam->x[m][j] * 
		fam->x[k][j] * 
		Gprim(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	hsum[i] = tmp;
    }
	
    for (i = 0; i < n; i++){
	tmp = 1.0;
	for (j = 0; j < fam->n; j++){
	    tmp *= P(fam->x_beta[j] + x[i] * fam->sigma,fam->y[j]);
	}
	res[i] *= tmp * (g1sum[i] * g2sum[i] + hsum[i]);
	/* x[i] = res[i];  Improve by 'memcpy' later */
    }

    memcpy(x, res, n * sizeof(double));

    Free(res);
    Free(g1sum);
    Free(g2sum);
    Free(hsum);
}

static void int_fun2s(double *x, 
		     int n, 
		     void *ex){
    /* Calculates one integrand (for one "family") in the expression     */
    /* of the partial second derivatives (hessian) wrt                   */
    /* 'beta[m], log(sigma)' of the log likelihood function.             */
    /* it is vectorizing, x is of length n. On input, x is the vector    */
    /* of points where int_fun should be evaluated, on output x contains */
    /* the n results.                                                    */
    /* NOTE: x[i] and fam->x[i] are different things here!! */

    double *res;
    int i, j;
    double zero = 0.0;
    double one = 1.0;
    int give_log = 0;

    double *g1sum, *g2sum, *hsum;
    double tmp;

    int m;

    Family *fam;

    fam = ex;

    m = fam->m;

    res = Calloc(n, double);
    g1sum = Calloc(n, double);
    g2sum = Calloc(n, double);
    hsum = Calloc(n, double);

    for (i = 0; i < n; i++) res[i] = dnorm(x[i], zero, one, give_log);

    /* First G-sum (m) */
    for (i = 0; i < n; i++){
	tmp = 0.0;
	for (j = 0; j < fam->n; j++){
	    tmp += fam->x[m][j] * 
		G(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	g1sum[i] = tmp;
}
	
    /* Second G-sum (sigma) */
    for (i = 0; i < n; i++){
	tmp = 0.0;
	for (j = 0; j < fam->n; j++){
	    tmp += G(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	g2sum[i] = tmp;
    }
	
    /* H-sum (m, k) */
    for (i = 0; i < n; i++){
	tmp = 0.0;
	for (j = 0; j < fam->n; j++){
	    tmp += 
		fam->x[m][j] * 
		Gprim(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	hsum[i] = tmp;
    }
	
    for (i = 0; i < n; i++){
	tmp = 1.0;
	for (j = 0; j < fam->n; j++){
	    tmp *= P(fam->x_beta[j] + x[i] * fam->sigma,fam->y[j]);
	}
	res[i] *= tmp * x[i] * fam->sigma *
	    (g1sum[i] * g2sum[i] + hsum[i]);
	/* x[i] = res[i];  Improve by 'memcpy' later */
    }

    memcpy(x, res, n * sizeof(double));

    Free(res);
    Free(g1sum);
    Free(g2sum);
    Free(hsum);
}

static void int_funss(double *x, 
		      int n, 
		      void *ex){
    /* Calculates one integrand (for one "family") in the expression     */
    /* of the partial second derivatives (hessian) wrt                   */
    /* 'log(sigma)' of the log likelihood function.                      */
    /* it is vectorizing, x is of length n. On input, x is the vector    */
    /* of points where int_fun should be evaluated, on output x contains */
    /* the n results.                                                    */
    /* NOTE: x[i] and fam->x[i] are different things here!! */

    double *res;
    int i, j;
    double zero = 0.0;
    double one = 1.0;
    int give_log = 0;

    double *gsum, *hsum;
    double tmp;

    Family *fam;

    fam = ex;

    res = Calloc(n, double);
    gsum = Calloc(n, double);
    hsum = Calloc(n, double);

    for (i = 0; i < n; i++) res[i] = dnorm(x[i], zero, one, give_log);

    /* G-sum (sigma) */
    for (i = 0; i < n; i++){
	tmp = 0.0;
	for (j = 0; j < fam->n; j++){
	    tmp += G(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	gsum[i] = tmp;
    }
	
    /* H-sum (m, k) */
    for (i = 0; i < n; i++){
	tmp = 0.0;
	for (j = 0; j < fam->n; j++){
	    tmp += 
		Gprim(fam->x_beta[j] + fam->sigma * x[i], fam->y[j]);
	}
	hsum[i] = tmp;
    }
	
    for (i = 0; i < n; i++){
	tmp = 1.0;
	for (j = 0; j < fam->n; j++){
	    tmp *= P(fam->x_beta[j] + x[i] * fam->sigma,fam->y[j]);
	}
	res[i] *= fam->sigma * x[i] *  tmp *
	    (gsum[i] * (1.0 + fam->sigma * x[i] * gsum[i]) +
	     fam->sigma * x[i] * hsum[i]);
    	/* x[i] = res[i];  Improve by 'memcpy' later */
    }

    memcpy(x, res, n * sizeof(double));

    Free(res);
    Free(gsum);
    Free(hsum);
}

static void update(int level,
		   int p, 
		   double *beta,
		   double *loglik, 
		   double *score,
		   double *hessian,
		   double *post_mode,
		   Family *fam,
		   int n_points,
		   double *weights,
		   double *zeros){
    
    double h;
    double *hb = NULL;
    double *hbb = NULL; /* (p+1) x (p+1) */
    double sigma;
    double tmp, tmp2;
    int i, m, k;

    /* For integrate */
    int inf = 2;
    double bound = 0.0;
    double epsabs = 0.000001; /* Check this!!!!!!!!!!!! */
    double epsrel = 0.000001; /* !!!!!!!!!!!!!!!!!!!!!! */
    double abserr;
    int neval;
    int ier;
    int limit = 100;
    int lenw;
    int *iwork;
    double *work;
    int last;
    double res;
    /* End For integrate */
    
    void *ex;

    int n;
    double **x;
    double *x_beta;
    int *y;

    int one = 1;

    double u, fu, sigma_hat;
    double abstol = 0.000001;
    double reltol = 0.000001;
    int maxit = 100;
    int trace = 0;
    int nREPORT = 1;
    int fncount, grcount;
    int fail = 0;
    int mask = 1;

    /* temporary gig: */
    double *x_val;
    double *wght;

    char *vmax;

    if (level < 0) return;
/*
    x_val = Calloc(n_points, double);
    wght = Calloc(n_points, double);
*/
    x_val = zeros;
    wght = weights;
/*
    for (i = 0; i < n_points; i++){
	x_val[i] = zeros[i] / 1.41421356237;
	wght[i] = weights[i] * 1.77245385091;
    }
*/
    ex = fam;

    n = fam->n;
    x = fam->x;
    x_beta = fam->x_beta;
    y = fam->y;

    lenw = 4 * limit;
    iwork = Calloc(limit, int);
    work = Calloc(lenw, double);

    sigma = exp(beta[p]);

    if (level > 0){
	hb = Calloc((p + 1), double);
	hbb = Calloc((p + 1) * (p + 1), double);
    }


    tmp = 0.0;
    for (i = 0; i < fam->n; i++){
	tmp += fam->x_beta[i];
    }

    tmp /= (double)fam->n;
    u = -tmp / fam->sigma; /* Start value for u */
    tmp = fun_fam(one, &u, ex);
/*    Rprintf("[fun:] Start value to vmmin: %f\n", tmp); */
    vmax = vmaxget();
    vmmin(one, &u, &fu, fun_fam, du_fun_fam, maxit, trace,
	  &mask, abstol, reltol, nREPORT, 
	  ex, &fncount, &grcount, &fail);
    vmaxset(vmax);
    
    du2_fun_fam(one, &u, &tmp, ex);

    sigma_hat = sqrt(1.0 / tmp);
    *post_mode = u;

/**********************************************************************/    
/* Calculate  h  for the loglik: */
    
/* We choose GHQ for the loglik....

    Rdqagi(int_fun, fam, &bound, &inf, &epsabs, &epsrel, &h,
	   &abserr, &neval, &ier, &limit, &lenw, &last,
	   iwork, work);
*/

/* GHQ: */

    h = 0.0;
    for (i = 0; i < n_points; i++){
	tmp = sqrt(2.0) * x_val[i] * sigma_hat + u;
	tmp2 = tmp;
	int_fun(&tmp2, one, ex);
	tmp2 *= wght[i] * exp(x_val[i] * x_val[i]);
	h += tmp2;
    }

    h = sigma_hat * sqrt(2) * h;

    /* Add into the loglik; We allow *loglik = -Inf!! */
    *loglik += log(h);

    if (level == 0) {
	Free(iwork);
	Free(work);
	return;
    }

    if (h <= 0){
	error("Unable do get likelihood function POSITIVE!!!!!!!!!\n");
    }


/***********************************************************************/    
/* First order derivatives, hb[]. */
    for (m = 0; m <= p; m++){ /* hb[m]; note w = log(sigma) INCLUDED! */
	if (m < p){

/* We don't use GHQ here...

	    hb[m] = 0.0;
	    fam->m = m;
	    for(i = 0; i < n_points; i++){
		tmp = sqrt(2.0) * x_val[i] * sigma_hat + u;
		tmp2 = tmp;
		int_fun1(&tmp2, one, ex);
		hb[m] += tmp2 * wght[i] 
		    * exp(x_val[i] * x_val[i]);
	    }
	    hb[m] = sigma_hat * sqrt(2) * hb[m];
*/
	    /*
	    Rprintf("hb[%d] med GHQ = %f\n", m, hb[m]);
	    */
	    fam->m = m;
	    Rdqagi(int_fun1, fam, &bound, &inf, &epsabs, &epsrel, &res,
		   &abserr, &neval, &ier, &limit, &lenw, &last,
		   iwork, work);
	    if (ier) Rprintf("[update] (Warning) ier = %d, hb[%d] = %f\n", 
			     ier, m, res);
	    hb[m] = res;
            /*
	    Rprintf("hb[%d] med int = %f\n\n", m, hb[m]);
	    */
	}else{ /* p == m */
/*
	    hb[m] = 0.0;
	    for(i = 0; i < n_points; i++){
		tmp = sqrt(2.0) * x_val[i] * sigma_hat + u;
		tmp2 = tmp;
		int_fun1s(&tmp2, one, ex);
		hb[m] += tmp2 * wght[i] 
		    * exp(x_val[i] * x_val[i]);
	    }
	    hb[m] = sigma_hat * sqrt(2) * hb[m];
*/
	    /*
	    Rprintf("hb[%d] med GHQ = %f\n", m, hb[m]);
	    */
	    Rdqagi(int_fun1s, fam, &bound, &inf, 
		   &epsabs, &epsrel, 
		   &res, &abserr, &neval, &ier, 
		   &limit, &lenw, &last,
		   iwork, work);
	    if (ier) Rprintf("[update] (Warning) ier = %d, hb[%d] = %f\n", 
			     ier, m, res);
	    hb[p] = res;
	    /*
	    Rprintf("hb[%d] med int = %f\n\n", p, hb[p]);
	    */
	}
    }

    
    /* Add into first derivatives: */
    for (m = 0; m <= p; m++){
	score[m] += hb[m] / h;
    }

    if (level == 1){
	Free(hbb);
	Free(hb);
	Free(iwork);
	Free(work);
	return;
    }
	
/***********************************************************************/
/* Done with first derivatives. On to the hessian: */

    /* First the pxp matrix of 'coefficients': */
    for (m = 0; m < p; m++){
	for (k = 0; k <= m; k++){
/*
	    hbb[m + k * (p + 1)] = 0.0;
	    fam->m = m;
	    fam->k = k;
	    for(i = 0; i < n_points; i++){
		tmp = sqrt(2.0) * x_val[i] * sigma_hat + u;
		tmp2 = tmp;
		int_fun2(&tmp2, one, ex);
		hbb[m + k * (p + 1)] += tmp2 * wght[i] 
		    * exp(x_val[i] * x_val[i]);
	    }
	    hbb[m + k * (p + 1)] = sigma_hat * sqrt(2) * 
		hbb[m + k * (p + 1)];
*/
	    /*
	    Rprintf("hbb[%d][%d] med GHQ = %f\n", m, k, hbb[m + k * (p + 1)]);
	    */
	    fam->m = m;
	    fam->k = k;
	    Rdqagi(int_fun2, fam, &bound, &inf, &epsabs, &epsrel, &res,
		   &abserr, &neval, &ier, &limit, &lenw, &last,
		   iwork, work);
	    if (ier) Rprintf("[update] (Warning) ier = %d, hbb[m][k] = %f\n", 
			     ier, m, k, res);

	    hbb[m + k * (p + 1)] = res;
	    /*
	    Rprintf("hbb[%d][%d] med int = %f\n", m, k, hbb[m + k * (p + 1)]);
	    */
	}
    }

    /* Then the "beta[m] with beta[p] = log(sigma), m = 0,...,(p-1) */
    m = p;
    for (k = 0; k < m; k++){ /* Changed k <= m --> k < m */ 
/*
	hbb[m + k * (p + 1)] = 0.0;
	for(i = 0; i < n_points; i++){
	    tmp = sqrt(2.0) * x_val[i] * sigma_hat + u;
	    tmp2 = tmp;
	    int_fun2s(&tmp2, one, ex);
	    hbb[m + k * (p + 1)] += tmp2 * wght[i] 
		* exp(x_val[i] * x_val[i]);
	}
	hbb[m + k * (p + 1)] = sigma_hat * sqrt(2) * 
	    hbb[m + k * (p + 1)];
*/
	/*
	Rprintf("hbb[%d][%d] med GHQ = %f\n", m, k, hbb[m + k * (p + 1)]);
	*/
	fam->m = k;
	fam->k = k;
	Rdqagi(int_fun2s, fam, &bound, &inf, &epsabs, &epsrel, &res,
	       &abserr, &neval, &ier, &limit, &lenw, &last,
	       iwork, work);
	if (ier) Rprintf("[update] (Warning) ier = %d, hb[%d][%d] = %f\n", 
			 ier, m, k, res);

	hbb[m + k * (p + 1)] = res;
	/*
	Rprintf("hbb[%d][%d] med int = %f\n", m, k, hbb[m + k * (p + 1)]);
	*/
    }

    /* And finally beta[p] with beta[p] (AKA log(sigma) = w): */
    k = p;
/*
    hbb[m + k * (p + 1)] = 0.0;

    for(i = 0; i < n_points; i++){
	tmp = sqrt(2.0) * x_val[i] * sigma_hat + u;
	tmp2 = tmp;
	int_funss(&tmp2, one, ex);
	hbb[m + k * (p + 1)] += tmp2 * wght[i] 
	    * exp(x_val[i] * x_val[i]);
    }
    hbb[m + k * (p + 1)] = sigma_hat * sqrt(2) * 
	hbb[m + k * (p + 1)];
*/
    /*
    Rprintf("hbb[%d][%d] med GHQ = %f\n", m, k, hbb[m + k * (p + 1)]);
    */
    Rdqagi(int_funss, fam, &bound, &inf, &epsabs, &epsrel, &res,
	   &abserr, &neval, &ier, &limit, &lenw, &last,
	   iwork, work);
	    if (ier)
		Rprintf("[update] (Warning) ier = %d, hbb[%d][%d] = %f\n", 
			ier, m, k, hb[m]);

    hbb[m + k * (p + 1)] = res;
/*
    Rprintf("hbb[%d][%d] med int = %f\n", m, k, hbb[m + k * (p + 1)]);    
*/
   /* Now, add it into the hessian (lower triangle): */
    for (m = 0; m <= p; m++){
	for (k = 0; k <= m; k++){
	    hessian[m + k * (p+1)] += hbb[m + k * (p+1)] / h  -   
	    	(hb[m] / h) * (hb[k] / h); /* 0 at the solution? No!!! */
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

    Free(hbb);
    Free(hb);
    Free(iwork);
    Free(work);
}

static double frail_mean(int level,
			 int p, 
			 double *beta,
			 double *loglik, 
			 double *score,
			 double *hessian,
			 Family *fam,
			 int n_points,
			 double *weights,
			 double *zeros){
    
    double h, h_mean;
    double sigma;
    double tmp;
    int i, j;
    
    sigma = exp(beta[p]);

/**********************************************************************/    
/* Calculate  h  for the loglik, and pip ("weights"): */
    h = 0.0;
    h_mean = 0.0;
    for (i = 0; i < n_points; i++){
	tmp = 1.0;
	for (j = 0; j < fam->n; j++){
	    tmp *= P(fam->x_beta[j] + zeros[i] * sigma, fam->y[j]);
	}
	h += tmp * weights[i];
	h_mean += tmp * weights[i] * zeros[i];
    }

    return ( h_mean / h);
}

void frail_fun(int pp1, 
		 double *beta,
		 void *ex){

    int start;
    double *gr = NULL;

    double tmp;
    int i, j;

    Exts *ext;
    Family *fam;

    double loglik;

    double *hessian = NULL;

    int level = 0;

    ext = ex;

    fam = Calloc(1, Family);
    fam->x = Calloc(ext->p, double *);
    fam->p = ext->p;
    fam->m = 0;
    fam->k = 0;
    fam->sigma = exp(beta[ext->p]);

    loglik = 0.0;

    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j][i];
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
	fam->n = ext->fam_size[i];
	fam->x_beta = ext->x_beta + start;
	fam->y = ext->y + start;
	for (j = 0; j < ext->p; j++)	
	    fam->x[j] = ext->x[j] + start;  /* check this! */

	ext->post_mean[i] = frail_mean(level,
				       ext->p,
				       beta,
				       &loglik,
				       gr,
				       hessian,
				       fam,
				       ext->n_points,
				       ext->weights,
				       ext->zeros);
	
	start += ext->fam_size[i];
    }
    Free(fam->x);
    Free(fam);
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
    double post_mode;
/*
    char trans = 'N';
    double alpha = 1.0;
    int one = 1;
*/
    double tmp;
    int i, j;

    Exts *ext;
    Family *fam;

    double loglik;

    double *hessian = NULL;

    int level = 0;

    ext = ex;

    fam = Calloc(1, Family);
    fam->x = Calloc(ext->p, double *);
    fam->p = ext->p;
    fam->m = 0;
    fam->k = 0;
    fam->sigma = exp(beta[ext->p]);

    loglik = 0.0;

    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j][i]; /* check this! */
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
	fam->n = ext->fam_size[i];
	fam->x_beta = ext->x_beta + start;
	fam->y = ext->y + start;
	for (j = 0; j < ext->p; j++)	
	    fam->x[j] = ext->x[j] + start;  /* check this! */
	update(level,
	       ext->p,
	       beta,
	       &loglik,
	       gr,
	       hessian,
	       &post_mode,
	       fam,
	       ext->n_points,
	       ext->weights,
	       ext->zeros);
	ext->post_mode[i] = post_mode;
	start += ext->fam_size[i];
    }
    Free(fam->x);
    Free(fam);
    /* printf("[fun]; loglik = %f\n", loglik); */ 
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

    Exts *ext;
    Family *fam;

    double tmp;
    double post_mode;

    int level = 1;

/* Note that we here trust 'ext->x_beta' to be properly updated !! */
/* In 'fun' */
 
    ext = ex;

    fam = Calloc(1, Family);
    fam->x = Calloc(ext->p, double *);
    fam->p = ext->p;
    fam->m = 0;
    fam->k = 0;
    fam->sigma = exp(beta[ext->p]);

    loglik = 0.0;

    for (k = 0; k < pp1; k++){
	gr[k] = 0.0;
    }

    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j][i]; 
	}
	ext->x_beta[i] = tmp;
    }
    
    start = 0;

    for (i = 0; i < ext->n_fam; i++){
	fam->n = ext->fam_size[i];
	fam->x_beta = ext->x_beta + start;
	fam->y = ext->y + start;
	for (j = 0; j < ext->p; j++)
	    fam->x[j] = ext->x[j] + start; 
	update(level,
	       ext->p,
	       beta,
	       &loglik,
	       gr,
	       hessian,
	       &post_mode,
	       fam,
	       ext->n_points,
	       ext->weights,
	       ext->zeros);
	ext->post_mode[i] = post_mode;
	start += ext->fam_size[i];
    }

    for (i = 0; i < pp1; i++){
	gr[i] = -gr[i]; /* Minimization! */
    }

    Free(fam->x);
    Free(fam);
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
    double post_mode;

    Exts *ext;
    Family *fam;

    int level = 2;

    ext = ex;

    fam = Calloc(1, Family);
    fam->x = Calloc(ext->p, double *);
    fam->p = ext->p;
    fam->m = 0;
    fam->k = 0;
    fam->sigma = exp(beta[ext->p]);

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
	    tmp += beta[j] * ext->x[j][i]; /* check this! */
	}
	ext->x_beta[i] = tmp;
    }
    
    start = 0;
    for (i = 0; i < ext->n_fam; i++){
	fam->n = ext->fam_size[i];
	fam->x_beta = ext->x_beta + start;
	fam->y = ext->y + start;
	for (j = 0; j < ext->p; j++)
	    fam->x[j] = ext->x[j] + start; /* check this! */

	update(level,
	       ext->p,
	       beta,
	       loglik,
	       gr,
	       hessian,
	       &post_mode,
	       fam,
	       ext->n_points,
	       ext->weights,
	       ext->zeros);
	ext->post_mode[i] = post_mode;
        start += ext->fam_size[i];
	}

   for (i = 0; i < pp1 * pp1; i++){
	hessian[i] = -hessian[i];
    }
   Free(fam->x);
   Free(fam);
}
