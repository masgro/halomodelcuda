#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "polifitgsl.c"

double polim(double x, double *a, int orden) {
  int i,j;
  double tmp, tmp1;
  tmp  = a[0];
  for(i = 1; i <= orden; i++) {
    tmp1 = a[i];
    for(j = 1; j <= i; j++) tmp1 *= x;
    tmp += tmp1;
  }
  return(tmp);
}

void fiteador(int nlines, float *x, float *y, float *e) {
  double *xx, *yy, *coeff;
  int DEGREE = 7;
  int k,kk;

  coeff = (double *) malloc(DEGREE*sizeof(double));
  xx = (double *) malloc(nlines*sizeof(double));
  yy = (double *) malloc(nlines*sizeof(double));

  kk = 0;
  for(k = 0; k < nlines; k++) {
    if(!isnormal(y[k])) {
      continue;
    }

    if(isnormal(e[k]) && e[k]/y[k] > 0.5) {
      continue;
    }

    xx[kk] = log10((double)x[k]);
    yy[kk] = log10((double)y[k]);
    kk++; 
  }

  xx = (double *) realloc(xx,kk*sizeof(double));
  yy = (double *) realloc(yy,kk*sizeof(double));
  
  polynomialfit(kk, DEGREE, xx, yy, coeff);
  for(k = 0; k < nlines; k++) {
    y[k] = (float)pow(10.0,polim(log10((double)x[k]),coeff,DEGREE-1));
  }
  free(xx);
  free(yy);
}
