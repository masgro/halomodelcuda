#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#ifndef ANGULO
#define ANGULO 0
#endif

#include "polifitgsl.c"

#define NDIR 3

//#define NPASOS   nlines

//#define NPASOS   20
//float h_radio_vector[NPASOS] = {0.119379,0.165646,0.229845,0.318927,0.442533,0.614046,0.852031,1.182253,
//																1.640459,2.276252,3.158458,4.382582,6.081139,8.438007,11.708325,16.246120,
//																22.542629,31.279476,43.402466,60.223957};

#define NPASOS   25
float h_radio_vector[NPASOS] = {1.1399923e-01,1.4815143e-01,1.9253504e-01,2.5021520e-01,3.2517537e-01,4.2259231e-01,5.4919362e-01,7.1372253e-01,9.2754149e-01,1.2054168e+00,
1.5665389e+00,2.0358467e+00,2.6457512e+00,3.4383726e+00,4.4684496e+00,5.8071198e+00,7.5468326e+00,9.8077326e+00,1.2745960e+01,1.6564426e+01,
2.1526846e+01,2.7975912e+01,3.6357018e+01,4.7248943e+01,6.1403919e+01};

double polim(double x, double *a, int orden){
  int i,j;
  double tmp, tmp1;

  tmp  = a[0];
  for(i = 1; i <= orden; i++){
    tmp1 = a[i];
    for(j = 1; j <= i; j++) tmp1 *= x;
    tmp += tmp1;
  }

  return(tmp);
}

struct function{
  float *x;
  float *y;
  float *e;
} f[NDIR];

int main(int argc, char **argv){
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

  for(j = 0; j < NDIR; j++){
    f[j].x = (float *) malloc(nlines*sizeof(float));
    f[j].y = (float *) malloc(nlines*sizeof(float));
    f[j].e = (float *) malloc(nlines*sizeof(float));

    for(k = 0; k < nlines; k++){
      f[j].y[k] =  0.0;
      f[j].e[k] = -1.E26;
    }

    for(i = 0; i < nrun; i++){
      sprintf(filename,"funcorr_%1d_%2s.%02d",j,term,i);
      pfin = fopen(filename,"r");
      assert(pfin!=NULL);
      for(k = 0; k < nlines; k++)
      {
        fscanf(pfin,"%e %e %e\n",&_tmp0,&_tmp1,&_tmp2);
        f[j].x[k] = _tmp0;
        f[j].y[k] += _tmp1;
        //f[j].y[k] = _tmp1;
        if(_tmp2 > f[j].e[k]) f[j].e[k] = _tmp2;
      }
      fclose(pfin);
    }

    for(k = 0; k < nlines; k++)
      f[j].y[k] /= (float)nrun;
  }

/*
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

    //free(f[j].x);
    //free(f[j].y);
    //f[j].x = (float *) malloc(NPASOS*sizeof(float));
    //f[j].y = (float *) malloc(NPASOS*sizeof(float));
    
    polynomialfit(kk,DEGREE,xx,yy,coeff);
    for(k = 0; k < NPASOS; k++)
    {
      //f[j].x[k] = pasos[k];
      f[j].y[k] = (float)pow(10.0,polim(log10((double)f[j].x[k]),coeff,DEGREE-1));
    }

    free(xx);
    free(yy);
  }
#endif
*/

	float mmin = atof(argv[2]);
	float mmax = atof(argv[3]);
	float abmedio = atof(argv[4]);
	float bcmedio = atof(argv[5]);
	float align_b = atof(argv[6]);
	float align_c = atof(argv[7]);
  sprintf(filename,"outputs/NEW_ALIGN/funcorr_%s_%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d.%s.dat",
                    term,mmin,mmax,ANGULO,bcmedio,abmedio,align_b,align_c,MAG,kind);
  pfout = fopen(filename,"w");
  assert(pfout != NULL);
  for(k = 0; k < NPASOS; k++){
    fprintf(pfout,"%e ",f[0].x[k]);
    for(j = 0; j < NDIR; j++)
      fprintf(pfout,"%e %e ",f[j].y[k],f[j].e[k]);
    fprintf(pfout,"\n");
  }
  fclose(pfout);
}
