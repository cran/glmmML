#ifndef FUN_H
#define FUN_H

#ifdef MATHLIB_STANDALONE
#include "testa.h"
#endif

#include "glmmml.h"



void frail_fun(int pp1, 
		 double *beta,
	       void *ex);

void mu_fun(int bdim, 
	    double *b, 
	    double *mu, 
	    void *ex);

double fun(int pp1, 
	   double *beta, 
	   void *ex);

void fun1(int pp1, 
	  double *beta,
	  double *gr,
	  void *ex);

void fun2(int pp1, 
	  double *beta,
	  double *loglik,
	  double *gr,
	  double *hessian,
	  void *ex);

void nr_opt(int bdim, double *beta, double *loglik, int *mask, 
	    Exts *ext, double epsilon, int maxit, int *info, int trace);

typedef double P_fun(double, int);

typedef double G_fun(double, int);

typedef double Gprim_fun(double, int);

typedef double Hprim_fun(double, int);

double P_logit(double x, int y); /* logit link */
    
double G_logit(double x, int y);

double Gprim_logit(double x, int y);

double Hprim_logit(double x, int y);

double Hbis_logit(double x, int y);

double P_cloglog(double x, int y);

double G_cloglog(double x, int y);

double Gprim_cloglog(double x, int y);

double Hprim_cloglog(double x, int y);

double P_poisson(double x, int y);

double G_poisson(double x, int y);

double Gprim_poisson(double x, int y);

double Hprim_poisson(double x, int y);

double Hbis_poisson(double x, int y);




#endif
