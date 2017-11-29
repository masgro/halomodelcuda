#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"
#include "deltas.h"
#include "polifitgsl.h"

int main(int argc, char **argv){
  read_parameterfile(argv[1]);

  set_units();

  initialize_powerspectrum();

  print_sigma();

  print_spec();

  xilin();

  exit(EXIT_SUCCESS);
}

void set_units(void){                /* ... set some units */
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;
  DeltaC = deltacri(Redshift);
  fprintf(stdout,"DELTA_C: %lf\n",DeltaC);
}

static double A, B, alpha, beta, V, gf;

double fnl(double x){                /* Peacock & Dodds formula */
  return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) /
                 (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)), 1 / beta);
}

void print_sigma(void){
  double chi2;
  double dm;
  double factor, norma;
  double dc2 = DeltaC*DeltaC;
  double *coeff;
  int    i, j;
  int    degree;
  char   filename[200];
  FILE   *pfout;
  
  factor = 3.0/4.0/M_PI/RHOCRIT/Omega;
  norma = Sigma8*Sigma8/TopHatSigma2(8.0 * (3.085678e24 / UnitLength_in_cm));

  sprintf(filename,"salida.dat");
  pfout = fopen(filename, "w");

  dm = (LMMAX - LMMIN)/(double)NTAB;
  for(i = 0; i < NTAB; i++){
    Mass[i]  = (double)i*dm + LMMIN;
    Scale[i] = pow(factor*pow(10.0,Mass[i]),1.0/3.0)*(3.085678e24 / UnitLength_in_cm);
    Sigma[i] = TopHatSigma2(Scale[i])*norma;
    Nu[i]    = sqrt(dc2/Sigma[i]);
    Nu[i]    = log10(Nu[i]);

    fprintf(pfout,"%E %E %E\n",Mass[i],Scale[i]/1000.0,Nu[i]);
  }
  fclose(pfout);
  
  double chi, chimin;
  int mydegree;
  chimin = 1.E26;
  for(i=1;i<30;i++){
    degree = i;
    coeff = (double *) malloc(degree*sizeof(double));
    polynomialfit(NTAB, degree, Mass, Nu, coeff);
    chi = 0.0;
    for(j=0;j<NTAB;j++){
      chi += pow((Nu[j]-polim(degree,coeff,Mass[j]))/Nu[j],2);
    }
    chi /= (double)NTAB;
    if((chi < chimin) && (chi > 1.E-8)){
      chimin = chi;
      mydegree = degree;
    }
    free(coeff);
  }
  fprintf(stdout,"\n%.7lf %d\n",chimin,mydegree);

  degree = 7;
  coeff = (double *) malloc(degree*sizeof(double));
  fprintf(stdout,"\nNu(mass)\n");
  polynomialfit(NTAB, degree, Mass, Nu, coeff);
  for(i = 0; i < degree; i++)
    fprintf(stdout,"%E\n",coeff[i]);

  pfout = fopen("Nu_M.dat","w");
  for(j = 0; j < NTAB; j++){
    fprintf(pfout,"%.7E %.7E %.7E\n",Mass[j],Nu[j],polim(degree,coeff,Mass[j]));
  }
  fclose(pfout);

  free(coeff);

  chimin = 1.E26;
  for(i = 1; i < 30; i++){
    degree = i;
    coeff = (double *) malloc(degree*sizeof(double));
    polynomialfit(NTAB, degree, Nu, Mass, coeff);
    chi = 0.0;
    for(j=0;j<NTAB;j++){
      chi += pow((Mass[j]-polim(degree,coeff,Nu[j]))/Mass[j],2);
    }
    chi /= (double)NTAB;
    if(chi<chimin && chi > 1.E-8){
      chimin = chi;
      mydegree = degree;
    }
    free(coeff);
  }
  fprintf(stdout,"\n%lf %d\n",chimin,mydegree);

  degree = 7;
  coeff = (double *) malloc(degree*sizeof(double));
  fprintf(stdout,"\nMass(nu)\n");
  polynomialfit(NTAB, degree, Nu, Mass, coeff);
  for(i = 0; i < degree; i++)
    fprintf(stdout,"%lf\n",coeff[i]);

  pfout = fopen("M_Nu.dat","w");
  for(j = 0; j < NTAB; j++){
    fprintf(pfout,"%.7E %.7E %.7E\n",Nu[j],Mass[j],polim(degree,coeff,Nu[j]));
  }
  fclose(pfout);

  free(coeff);
}

double polim(int degree, double *coeff, double x){
  double tmp = coeff[0];
  int    i;
  for(i = 1; i < degree; i++)
    tmp += coeff[i]*pow(x,i);
  return(tmp);
}

void xilin(void){
  double k, po, dl, kstart, kend, po2, po1;
  double k1, k2, dk = 1.001;
  double rmin, rmax, dr, rtmp, xitmp;
  double *xi, *r;

  double result, abserr;
  gsl_integration_workspace *w, *cycle_w;
  gsl_integration_qawo_table *wf;
  gsl_function F;

  int    i, j;
  FILE *pfout;

  r = (double *) malloc(NTAB*sizeof(double));
  xi = (double *) malloc(NTAB*sizeof(double));

  rmin = 0.01;
  rmax = 100.0;

  kstart = 2 * PI / (1000.0 * (3.085678e24 / UnitLength_in_cm));     /* 1000 Mpc/h */
  kend = 2 * PI / (0.001 * (3.085678e24 / UnitLength_in_cm));        /* 0.001 Mpc/h */

  w = gsl_integration_workspace_alloc(100000);
  cycle_w = gsl_integration_workspace_alloc(100000);

  rtmp = rmin*(3.085678e24 / UnitLength_in_cm);
  wf = gsl_integration_qawo_table_alloc(rtmp, 1.E5, GSL_INTEG_SINE, 100);

  F.function = &integrando;

  dr = (log10(rmax) - log10(rmin))/(double)NTAB;
  for(i = 0, j = 0; i < NTAB; i++){
    rtmp = (double)i*dr + log10(rmin);
    rtmp = pow(10.0, rtmp)*(3.085678e24 / UnitLength_in_cm);

    gsl_integration_qawo_table_set(wf, rtmp, 1.E5, GSL_INTEG_SINE);

    F.params = &rtmp;
    gsl_integration_qawf(&F, kstart, 1.0e-8, 100000, w, cycle_w, wf, &result, &abserr);

    if(result < 0.0) continue;
    r[j]  = log10(rtmp/1000.0);
    xi[j] = log10(result);
    j++;
  }

  int nlines = j; 

  r  = (double *) realloc(r, nlines*sizeof(double));
  xi = (double *) realloc(xi, nlines*sizeof(double));

  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(w);
  gsl_integration_workspace_free(cycle_w);

  double *coeff;
  int degree;
  double chi, chimin;
  int mydegree;
  chimin = 1.E26;
  for(i = 1; i < 30; i++){
    degree = i;

    coeff = (double *) malloc(degree*sizeof(double));

    polynomialfit(j, degree, r, xi, coeff);

    chi = 0.0;
    for(j = 0; j < nlines; j++){
      chi += pow((xi[j]-polim(degree,coeff,r[j]))/xi[j],2.);
    }
    chi /= (double)nlines;

    if((chi < chimin) && (chi > 1.E-8)){
      chimin = chi;
      mydegree = degree;
    }

    free(coeff);
  }
  fprintf(stdout,"\n%lf %d\n",chimin,mydegree);

  degree = 7;
  coeff = (double *) malloc(degree*sizeof(double));
  polynomialfit(j, degree, r, xi, coeff);

  pfout = fopen("xi.dat","w");
  for(j=0;j<nlines;j++)
    fprintf(pfout,"%lf %lf %lf\n",r[j],xi[j],polim(degree,coeff,r[j]));
  fclose(pfout);


  fprintf(stdout,"\nXiLin(r)\n");
  for(i = 0; i < degree; i++)
    fprintf(stdout,"%lf\n",coeff[i]);
}

double integrando(double k, void *param)
{
  double tmp;
  double r = *(double *) param;

  tmp  = 4.0*M_PI*PowerSpec(k);
  tmp *= k*k;
  tmp /= (r*k);

  return(tmp);
}

void print_spec(void)
{
  double k, knl, po, dl, dnl, neff, kf, kstart, kend, po2, po1, DDD;
  char buf[1000];
  FILE *fd;

  sprintf(buf, "%s/inputspec_%s.txt", OutputDir, FileBase);

  fd = fopen(buf, "w");

  gf = GrowthFactor(0.001, 1.0) / (1.0 / 0.001);

  DDD = GrowthFactor(1.0 / (Redshift + 1), 1.0);

  fprintf(fd, "%12g %12g\n", Redshift, DDD); /* print actual starting redshift and 
                                             linear growth factor for this cosmology */

  kstart = 2 * PI / (1000.0 * (3.085678e24 / UnitLength_in_cm));        /* 1000 Mpc/h */
  kend = 2 * PI / (0.001 * (3.085678e24 / UnitLength_in_cm));        /* 0.001 Mpc/h */

  for(k = kstart; k < kend; k *= 1.025)
    {
      po = PowerSpec(k);
      dl = 4.0 * PI * k * k * k * po;

      kf = 0.5;

      po2 = PowerSpec(1.001 * k * kf);
      po1 = PowerSpec(k * kf);

      if(po != 0 && po1 != 0 && po2 != 0)
        {
          neff = (log(po2) - log(po1)) / (log(1.001 * k * kf) - log(k * kf));

          if(1 + neff / 3 > 0)
            {
              A = 0.482 * pow(1 + neff / 3, -0.947);
              B = 0.226 * pow(1 + neff / 3, -1.778);
              alpha = 3.310 * pow(1 + neff / 3, -0.244);
              beta = 0.862 * pow(1 + neff / 3, -0.287);
              V = 11.55 * pow(1 + neff / 3, -0.423) * 1.2;

              dnl = fnl(dl);
              knl = k * pow(1 + dnl, 1.0 / 3);
            }
          else
            {
              dnl = 0;
              knl = 0;
            }
        }
      else
        {
          dnl = 0;
          knl = 0;
        }

      fprintf(fd, "%12g %12g    %12g %12g\n", k, dl, knl, dnl);
    }
  fclose(fd);
}
