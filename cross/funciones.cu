__device__ float T1h(float *r, float *x){
  float tmp;
  float lgm      = x[0];
  float ab       = x[1];
  float bc       = x[2];
  float m,nu;

  nu = Nu_M(lgm);
  nu = powf(10.0f,nu);
  m  = powf(10.0f,lgm);

  //tmp  = u_esferico(r,lgm,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
  tmp  = u(r,bc,ab,lgm,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
#ifdef CG
	tmp *= RHOMEDIO/m;
  tmp *= mom1(lgm);
#endif
  tmp *= f_nu(nu)*nu*dlogNu_dlogM(lgm)*FMNORM;
  /*FUNCION DE MASA FIT*/
  //tmp *= n(lgm)*m*m/RHOMEDIO;
  tmp *= forma(bc,ab);
  tmp *= LOGE10; //integracion logaritmica en m1
  return(tmp);
}

__device__ float T2h(float *r, float *x){
  float tmp;
  float lgm1, lgm2;
  float m1, m2;
  float nu1, nu2;
  float ab2, bc2;
  float ab1, bc1;
  float xc, yc, zc;
  float y[3],yy[3];

  /***********
  x[0]       = M1
  x[2,1]     = bc1,ab1 //bugy
  x[3]       = M2
  x[5,4]     = bc2,ab2 //bugy
  x[6,7,8]   = x,y,z
  x[9,10,11] = alpha,cos(beta),gamma -> angulos de euler
  **********/

  lgm1 = x[0];
   ab1 = x[1];
   bc1 = x[2];

  lgm2 = x[3];
   ab2 = x[4];
   bc2 = x[5];

  y[0] = x[6]; 
  y[1] = x[7]; 
  y[2] = x[8];

   m1 = powf(10.0f,lgm1);
  nu1 = Nu_M(lgm1);
  nu1 = powf(10.0f,nu1);

   m2 = powf(10.0f,lgm2);
  nu2 = Nu_M(lgm2);
  nu2 = powf(10.0f,nu2);

  xc = y[0] + r[0];
  yc = y[1] + r[1];
  zc = y[2] + r[2];

  rotate(y,x[9],x[10],x[11],yy);

  //tmp  = u_esferico(yy,lgm2,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
  tmp  = u(yy,bc2,ab2,lgm2,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
  tmp *= (RHOMEDIO/m1);
  #ifdef CG
  tmp *= (RHOMEDIO/m2);
  tmp *= mom1(lgm2);
  #endif
  tmp *= f_nu(nu2)*nu2*dlogNu_dlogM(lgm2)*FMNORM;
  /*FUNCION DE MASA FIT*/
  //tmp *= n(lgm2)*m2*m2/RHOMEDIO;
  tmp *= forma(bc2,ab2);
  tmp *= NORMA_ORIENTACION_CROSS; //factor para normalizar la func. de probabilidad de
                                    //orientacion del halo vecino.
  tmp *= f_nu(nu1)*nu1*dlogNu_dlogM(lgm1)*FMNORM;
  /*FUNCION DE MASA FIT*/
  //tmp *= n(lgm1)*m1*m1/RHOMEDIO;
  tmp *= forma(bc1,ab1);
  tmp *= bias(nu1);
  tmp *= bias(nu2);
  tmp *= align_masa(zc,yc,xc);
  tmp *= LOGE10;  //integracion logaritmica en m2
  tmp *= LOGE10;  //integracion logaritmica en m1
  tmp *= xilin(xc,yc,zc);
  return(tmp);
}
