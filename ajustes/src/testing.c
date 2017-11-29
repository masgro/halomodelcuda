#include <stdio.h>
#include <stdlib.h>
 
#include "polifitgsl.h"
 
#define NP 100
double *x;
double *y;
 
#define DEGREE 7
double coeff[DEGREE];
 
int main()
{
  int i, j;
 
  x = (double *) malloc(NP*sizeof(double));
  y = (double *) malloc(NP*sizeof(double));

  srand48(192837465L);

  for(i = 0; i < DEGREE; i++){
    coeff[i] = drand48()*2.0 - 1.0;
    fprintf(stdout,"%lf\n",coeff[i]);
  }
  fprintf(stdout,"\n");

  for(i = 0; i < NP; i++){
    x[i] = 0.01*(double)i;
    y[i] = 0.0;
    for(j = 0; j < DEGREE; j++)
      y[i] += coeff[j]*pow(x[i],j);
  }

  for(i = 0; i < DEGREE; i++) coeff[i] = 0.0;

  polynomialfit(NP, DEGREE, x, y, coeff);
  for(i = 0; i < DEGREE; i++) {
    printf("%lf\n", coeff[i]);
  }
  return 0;
}
