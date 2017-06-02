__host__ __device__ float cambio_de_coordenadas(float ab, float bc, float costheta, 
                                               float phi, float *r){
  float psi,cp,sp;
  float ct,st,cf,sf,st2,ct2,sf2,cf2;
  float j,k,l;
  float xt,yt,zt;
  float t,u;
  float theta;
  float f;

  ///************/
  //t = ab*ab;
  //t = 1.0f/t;
  //u = bc*bc;
  //u = 1.0f/u;

  //ct = costheta;
  //ct2 = ct*ct;
  //st2 = 1.0f - ct2;
  //st = sqrtf(st2);

  //cf = cosf(phi);
  //sf = sinf(phi);
  //sf2 = sf*sf;
  //cf2 = cf*cf;

  //A = ct2*(sf2*t*u + cf2*u) + st2*t*u*u;
  //B = ct*sinf(2.0f*phi)*(t*u - u);
  //C = (sf2*u + cf2*t*u);

  //psi = 0.5f*atan2f(B,(A-C));

  //f = st*st*(cf2*t*u+u*sf2) + ct2;
  ///************/


  /**** NEW: siguiendo STARK********/
  //theta = acosf(costheta);
  ct = costheta;
  ct2 = ct*ct;
  st = 1.0f-ct2;
  st2 = st*st;

  cf = cosf(phi);
  sf = sinf(phi);
  sf2 = sf*sf;
  cf2 = cf*cf;

  t = 1./(ab*bc);
  u = 1./(bc);

  j = st2*t*t*u*u + t*t*cf2*ct2 + u*u*sf2*ct2;
  k = (u*u - t*t)*sf*cf*ct;
  l = t*t*sf2 + u*u*cf2;

  f = st2*(t*t*sf2 + u*u*cf2) + ct2;

  psi = 0.5f*atanf(2.0*k/(j-l));
  /************/

  //if(((A - C)*cosf(2.0f*psi) + B*sinf(2.0f*psi)) > 0.0f)
  //  psi = psi - PI_CUDA*0.5f;

  cp = cosf(psi);
  sp = sinf(psi);
  
  xt = r[0];
  yt = r[1];
  zt = r[2];

  /** X = A0^{-1} * A1^{-1} * A2^{-1} * X" **/
  r[0] = xt*(cf*cp-sf*ct*sp)-yt*(cf*sp+sf*ct*cp)+sf*st*zt;
  r[1] = xt*(sf*cp+cf*ct*sp)-yt*(sf*sp-cf*ct*cp)-cf*st*zt;
  r[2] = st*sp*xt+st*cp*yt+ct*zt;

  return(f);
}

__device__ float T1h(float *r, float *x){
  float tmp;
  float lgm      = x[0];
  float ab       = x[1];
  float bc       = x[2];
  float costheta = x[3];
  float phi      = x[4];
  float m,nu;
  float f;

  nu = Nu_M(lgm);
  nu = powf(10.0f,nu);
  m  = powf(10.0f,lgm);

  f = cambio_de_coordenadas(ab,bc,costheta,phi,r);
  
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
  tmp *= NORMA_ORIENTACION_PIZZA;
  tmp *= LOGE10; //integracion logaritmica en m1
  return(tmp);
}

__device__ float T2h(float *r, float *x){
  float tmp;
  float lgm1, lgm2;
  float m1, m2;
  float nu1, nu2;
  //float ab2, bc2;
  float ab1, bc1;
  float costheta1, phi1;
  //float costheta2, phi2;
  float xc, yc, zc;
  float y[3];
  float f;

  /***********
  x[0]       = M1
  x[2,1]     = bc1,ab1 //bugy
  x[3,4]     = e11,e12
  x[5]       = M2
  x[6,7]     = bc2,ab2 //bugy
  x[8,9]     = e21,e22
  x[10,11,12]   = x,y,z
  **********/

  lgm1 = x[0];
   ab1 = x[1];
   bc1 = x[2];

  costheta1 = x[3];
  phi1      = x[4];

  lgm2      = x[5];
  //ab2       = x[6];
  //bc2       = x[7];
  //costheta2 = x[8];
  //phi2      = x[9];
  y[0]      = x[10];
  y[1]      = x[11];
  y[2]      = x[12];

   m1 = powf(10.0f,lgm1);
  nu1 = Nu_M(lgm1);
  nu1 = powf(10.0f,nu1);

   m2 = powf(10.0f,lgm2);
  nu2 = Nu_M(lgm2);
  nu2 = powf(10.0f,nu2);

  f = cambio_de_coordenadas(ab1,bc1,costheta1,phi1,r);

  xc = y[0] + r[0];
  yc = y[1] + r[1];
  zc = y[2] + r[2];

  //tmp  = u_esferico(y,lgm2,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
  tmp  = u(y,1.f,1.f,lgm2,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
  tmp *= (RHOMEDIO/m1);
  #ifdef CG
  tmp *= (RHOMEDIO/m2);
  tmp *= mom1(lgm2);
  #endif
  tmp *= f_nu(nu2)*nu2*dlogNu_dlogM(lgm2)*FMNORM;
  /*FUNCION DE MASA FIT*/
  //tmp *= n(lgm2)*m2*m2/RHOMEDIO;
  //tmp *= forma(bc2,ab2);
  //tmp *= NORMA_ORIENTACION; //factor para normalizar la func. de probabilidad de
                            //orientacion del halo vecino.
  tmp *= f_nu(nu1)*nu1*dlogNu_dlogM(lgm1)*FMNORM;
  /*FUNCION DE MASA FIT*/
  //tmp *= n(lgm1)*m1*m1/RHOMEDIO;
  tmp *= forma(bc1,ab1);
  tmp *= NORMA_ORIENTACION_PIZZA;
  tmp *= bias(nu1);
  tmp *= bias(nu2);
  tmp *= align_masa(zc,yc,xc);
  tmp *= LOGE10;  //integracion logaritmica en m2
  tmp *= LOGE10;  //integracion logaritmica en m1
  tmp *= xilin(xc,yc,zc);
  return(tmp);
}

__global__ void proyecta_xilin(curandState *state, float r){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int it = threadIdx.x;
  int   i,j,k;
  float value, sigma;
  float p[3], tmp;

  __shared__ float s_value[THREADS_PER_BLOCK];
  __shared__ float s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState xx = state[tid];

  /*Esperan hasta que todos terminen*/
  __syncthreads();

  /*Tira un numero random x de dimension NDIM*/
  value = 0.0f; sigma = 0.0f;
  for(j = 0; j < LAZOS; j++){
    k = 1;
    while(k != 0){
      p[0] = r;
      p[1] = 0.0f;

      p[2] = curand_normal(&xx); /*Linea de la visual*/

      tmp  = SQRT_TWOPI_CUDA*RANGO/exp(-p[2]*p[2]*0.5f);
      p[2] *= RANGO;

      tmp *= xilin(p[0],p[1],p[2]);
      tmp *= bias(powf(10.0f,Nu_M(14.0)));

      value += tmp;
      sigma += tmp*tmp;

      k = isnan(tmp);
    }
  }

  /*Guarda en la memoria compartida*/
  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();

  if(it == 0){
    value = 0.0f; sigma = 0.0f;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }
    j = blockIdx.x;
    /*Suma a la memoria global*/
    d_int[j] = value;
    d_sig[j] = sigma;
  }

  /*Guarda el estado del RNG*/
  state[tid] = xx;
}
