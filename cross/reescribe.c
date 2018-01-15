#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#ifndef ANGULO
#define ANGULO 0
#endif

#include "polifitgsl.c"

#define NDIR 3

double polim(double x, double *a, int orden)
{
  int i,j;
  double tmp, tmp1;

  tmp  = a[0];
  for(i = 1; i <= orden; i++)
  {
    tmp1 = a[i];
    for(j = 1; j <= i; j++) tmp1 *= x;
    tmp += tmp1;
  }

  return(tmp);
}

struct function
{
  float *x;
  float *y;
  float *e;
} f[NDIR];

int main(int argc, char **argv)
{
  char filename[200];
  FILE *pfin, *pfout;
  float _tmp0, _tmp1, _tmp2;
  int nrun = atoi(argv[1]);
  int i,j,k;
  int nlines; 

  char cmd[200], term[200], kind[200];

#ifdef GG
  sprintf(kind,"GG");
#endif
#ifdef CG
  sprintf(kind,"CG");
#endif
#ifdef MM
  sprintf(kind,"MM");
#endif
#ifdef HM
  sprintf(kind,"HM");
#endif

#ifdef DOSH
  sprintf(term,"2h");
#else
  sprintf(term,"1h");
#endif

  sprintf(filename,"funcorr_%1d_%s.%02d",0,term,0);

  sprintf(cmd,"wc -l %s\n",filename);
  pfin = popen(cmd,"r");
  assert(pfin != NULL);
  fscanf(pfin,"%d",&nlines);
  pclose(pfin);

  for(j = 0; j < NDIR; j++)
  {
    f[j].x = (float *) malloc(nlines*sizeof(float));
    f[j].y = (float *) malloc(nlines*sizeof(float));
    f[j].e = (float *) malloc(nlines*sizeof(float));

    for(k = 0; k < nlines; k++){
      f[j].y[k] = 1.E26;
      f[j].e[k] = -1.E26;
    }

    for(i = 0; i < nrun; i++)
    {
      sprintf(filename,"funcorr_%1d_%2s.%02d",j,term,i);
      pfin = fopen(filename,"r");
      assert(pfin!=NULL);
      for(k = 0; k < nlines; k++)
      {
        fscanf(pfin,"%e %e %e\n",&_tmp0,&_tmp1,&_tmp2);
        f[j].x[k] = _tmp0;
        //if(_tmp1 < f[j].y[k]) f[j].y[k] = _tmp1;
        f[j].y[k] = _tmp1;
        if(_tmp2 > f[j].e[k]) f[j].e[k] = _tmp2;
      }
      fclose(pfin);
    }
  }

#ifdef DOSH
  double *xx, *yy, *coeff;
  int DEGREE = 7;
  int kk;

  coeff = (double *) malloc(DEGREE*sizeof(double));

  for(j = 0; j < NDIR; j++)
  {
    xx = (double *) malloc(nlines*sizeof(double));
    yy = (double *) malloc(nlines*sizeof(double));

    kk = 0;
    for(k = 0; k < nlines; k++)
    {
      if(!isnormal(f[j].y[k])) continue;
      if(isnormal(f[j].e[k]) && f[j].e[k]/f[j].y[k] > 0.5) continue;

      xx[kk] = log10((double)f[j].x[k]);
      yy[kk] = log10((double)f[j].y[k]);
      kk++; 
    }

    xx = (double *) realloc(xx,kk*sizeof(double));
    yy = (double *) realloc(yy,kk*sizeof(double));
    
    polynomialfit(kk,DEGREE,xx,yy,coeff);
    for(k = 0; k < nlines; k++)
    {
      f[j].y[k] = (float)pow(10.0,polim(log10((double)f[j].x[k]),coeff,DEGREE-1));
    }

    free(xx);
    free(yy);
  }
#endif

  //sprintf(filename,"outputs/funcorr_%s_%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d.%s.%02d",
  //                  term,CENTROS_MASA_MIN,CENTROS_MASA_MAX,ANGULO,
  //                  BCMEDIO,ABMEDIO,ALIGN_B,ALIGN_C,MAG,kind,NUM);
  sprintf(filename,"outputs/funcorr_%s_%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d.%s.dat",
                    term,CENTROS_MASA_MIN,CENTROS_MASA_MAX,ANGULO,
                    BCMEDIO,ABMEDIO,ALIGN_B,ALIGN_C,MAG,kind);
  pfout = fopen(filename,"w");
  assert(pfout != NULL);
  for(k = 0; k < nlines; k++)
  {
    fprintf(pfout,"%e ",f[0].x[k]);
    for(j = 0; j < NDIR; j++)
      fprintf(pfout,"%e %e ",f[j].y[k],f[j].e[k]);
    fprintf(pfout,"\n");
  }
  fclose(pfout);
}
