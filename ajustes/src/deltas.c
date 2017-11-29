#include <math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "deltas.h"


double deltacri( double z )
{
  double q ;   
  if((Omega - 1.0) < 1.0E-8 & OmegaLambda < 1.0E-8)
  {
    q = (3./20.)*pow(12.*M_PI,2./3.)/lingrowth(z) ;
  }
  if(Omega < 1. & OmegaLambda < 1.0E-8)
  {
    q = (3./20.)*pow(12.*M_PI,2./3.)*pow(Omega_z(z),0.0185)/lingrowth(z) ;
  }
  if(Omega < 1. & (Omega + OmegaLambda - 1.0) < 1.0E-8)
  {
    q = (3./20.)*pow(12.*M_PI,2./3.)*pow(Omega_z(z),0.0055)/lingrowth(z) ;
  }
  return q;
}

double Omega_z( double z )
{
  double q ;
  if((Omega - 1.) < 1.0E-8)
  {
    q = Omega*(1.+z);
  }
  if(Omega < 1. & OmegaLambda < 1.0E-8 )
  {
    q = Omega*(1.+z)/(1+Omega*z) ;
  }
  if(Omega < 1. & (Omega + OmegaLambda - 1.) < 1.0E-8)
  {
    q = Omega*pow(1.+z,3.)/(1.-Omega+Omega*pow(1.+z,3.)) ;
  }
  return(q) ;
}

double lingrowth( double z )
{
  double w, w0, y, y0;
  double b;
  w0 = 1./Omega - 1 ;

  if((Omega - 1.0) < 1.0E-8)
  {
    b = 1./(1.+z) ;
  }
  if(Omega < 1. & (OmegaLambda < 1.0E-8))
  {
    w  = w0/(1. + z) ;
    b  = f1(w) / f1(w0) ;
  }
  if(Omega < 1. & (Omega + OmegaLambda - 1.0) < 1.0E-8)
  {
    y0 = pow(2.*w0,1./3.) ;
    y  = y0/(1. + z) ;
    b  = f2(y)*f3(y)/f2(y0)/f3(y0) ;
  }
  return b;
}
  
double f1(double u)
{   
  double q ;
  q = 1+3./u+3*pow(1+u,0.5)/pow(u,1.5)*log(pow(1.+u,0.5)-pow(u,0.5)) ;
  return q ;
}

double f2(double u)
{
  double q ;
  q = pow(pow(u,3.)+2.,0.5)/pow(u,1.5) ;
  return q ;
}

double f3(double u)
{
  double resultado, error ;
  size_t neval ;
  gsl_function F ;
  F.function=&intf3 ;  
  gsl_integration_qng(&F,0.0,u,1.0E-7,1.0E-7,&resultado,&error,&neval) ;
  return(resultado);
}  

double intf3(double u, void *par)
{
  double q ; 
  q = pow(u/(pow(u,3.)+2.),1.5) ;
  return q ;
}
