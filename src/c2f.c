/* Makes some C routines (mainly from R) available in Fortran */

#include <R.h>

void F77_SUB(randstart)(void) { GetRNGstate();}
void F77_SUB(randend)  (void) { PutRNGstate();}
void F77_SUB(ranf)(double *x) { *x = unif_rand();}
void F77_SUB(rvsort)(double *x, int *index, int *n) {
	revsort(x, index, *n);
}
