__host__ __device__ float polim(float x, float *a, int n){
  int i,j;
  float tmp, tmp1;
  tmp  = a[0];
  for(i = 1; i <= n; i++){
    tmp1 = a[i];
    for(j = 1; j <= i; j++) tmp1 *= x;
    tmp += tmp1;
  }
  return(tmp);
}

/**Log10(nu) en funcion de Log10(M)**/
__host__ __device__ float Nu_M(float m){
#ifndef __CUDA_ARCH__
  return(polim(m,h_coef_nu_m,h_coef_nu_m_n-1));
#else
  return(polim(m,d_coef_nu_m,d_coef_nu_m_n-1));
#endif
}

__host__ __device__ float dlogNu_dlogM(float m){
#ifndef __CUDA_ARCH__
  return(polim(m,h_coef_dlognu_dlogm,h_coef_dlognu_dlogm_n-1));
#else
  return(polim(m,d_coef_dlognu_dlogm,d_coef_dlognu_dlogm_n-1));
#endif
}

__host__ __device__ float M_Nu(float nu){
#ifndef __CUDA_ARCH__
  return(polim(nu,h_coef_m_nu,h_coef_m_nu_n-1));
#else
  return(polim(nu,d_coef_m_nu,d_coef_m_nu_n-1));
#endif
}

__host__ __device__ float xilin(float x, float y, float z){
  float tmp;
  tmp = sqrtf(x*x + y*y + z*z);
  tmp = log10f(tmp);
#ifndef __CUDA_ARCH__
  tmp = polim(tmp,h_coef_xilin,h_coef_xilin_n-1);
#else
  tmp = polim(tmp,d_coef_xilin,d_coef_xilin_n-1);
#endif
  tmp = powf(10.0f,tmp);
  return(tmp);
}

__device__ float db_gauss(float f, float m1, float m2){
  float fi;
  fi  = exp(-(f - m1)*(f - m1)*0.5f/m2);
  fi /= SQRT_TWOPI_CUDA;
  fi /= sqrtf(m2);
  return(fi);
}

__device__ float forma(float bc, float ab){
  float prob;
  float dbc,dab;
  float m1bc,m2bc;
  float m1ab,m2ab;

  m1bc = BCMEDIO; m2bc = 0.1f*0.1f; //sigma^2
  dbc  = db_gauss(bc,m1bc,m2bc)*(1.0f-bc)*bc;

  m1ab = ABMEDIO; m2ab = 0.1f*0.1f;  //sigma^2
  dab  = db_gauss(ab,m1ab,m2ab)*(1.0f-ab)*ab;

  prob = dab*dbc;
  if(prob < 0.0f) prob = 0.0f;

  return(prob);
}

__host__ __device__ float bias(float nu){
  float q;
  float tmp;

  tmp = BIAS_ALPHA*nu*nu;
  q  = 1.0f;
  q += (tmp - 1.0f)/DELTAC;
  q += 2.0f*BIAS_POTEN/(DELTAC*(1.0f + powf(tmp,BIAS_POTEN)));

  //tmp = BIAS_ALPHA*nu*nu;
  //q = BIAS_SQRT_ALPHA*tmp;
  //q += (BIAS_SQRT_ALPHA*BIAS_BETA*powf(tmp,1.0f - BIAS_GAMMA));
  //q -= (powf(tmp,BIAS_GAMMA)/(powf(tmp,BIAS_GAMMA)+BIAS_BETA*(1.0f - BIAS_GAMMA)*(1.0f - BIAS_GAMMA*0.5f)));
  //q /= (BIAS_SQRT_ALPHA*DELTAC);
  //q += 1.0f;

  return(q);
}

__host__ __device__ float f_nu(float nu){
  float q;
  float x;
  x  = ALPHA*nu*nu;
  q  = 2.0f;
  q *= (1.0f + 1.0f/powf(x,POTEN));
  q *= sqrtf(x*0.5f/PI_CUDA);
  q *= expf(-0.5f*x);
  q /= nu;
  return(q);
}

__device__ float align(float theta, float phi){
  float p;
  p = ALIGN_C*expf(-(theta - PI_CUDA*0.5f)*(theta - PI_CUDA*0.5f)*0.5f/(ALIGN_B*ALIGN_B));
  p *= expf(-phi*phi*0.5f/(ALIGN_B*ALIGN_B));
  p += expf(-theta*theta*0.5f/(ALIGN_B*ALIGN_B));
  return(p);
}

__device__ float align_masa(float x, float y, float z){
  float costheta;
  float r;

  r = x*x + y*y + z*z;
  r = sqrtf(r);

  if(r < 1.E-7)
    costheta = 1.0f;
  else
    costheta = fabsf(x)/r;

  return(align(acosf(costheta),atan2f(fabsf(z),fabsf(y))));
}

__host__ __device__ float ncen(float lgm){
  float q;
  q = 0.5f;
  q *= (1.0f + (float)erf((lgm - NCEN_LOGMMIN)/NCEN_SIGMA));
  return(q);
}

__host__ __device__ float nsat(float lgm){
  float q;
  q = (powf(10.0f,lgm) - NSAT_LOGM0);
  q /= NSAT_LOGM1;
  q = (q > 0.0f)? powf(q,NSAT_ALPHA) : 0.0f;

  return(q);
}

__host__ __device__ float f_projected(float r)
{
  float q, p;

  p = 1.0f - r*r;
  if(p > 0.0f)
  {
    q = sqrtf((1.0f - r)/(1.0f + r));
    q = -1.0f + 2.0f/sqrtf(p)*atanh(q);
    q = q/p;
  }
  else
  {
    if(p < 0.0f)
    {
      q = sqrtf((r - 1.0f)/(r + 1.0f));
      q = 1.0f - 2.0f/sqrtf(-p)*atan(q);
      q = q/(-p);
    }
    else
      q = 1.0f/3.0f;
  }

  return q;
}

__host__ __device__ float profile_1h(float *x, float lgm, float ab, float bc)
{
  float cvir,rvir,ce,Deltae,dc,rs,r;
  float rho;

  if(bc < 0.1f) bc = 0.1f;
  if(ab < 0.1f) ab = 0.1f;

  /** parametro de concentracion **/
  cvir = 1.020f - 0.109f*(lgm - 12.0f); 
  cvir = powf(10.0f,cvir);
  /********************************/

  /** radio virial esferico **/
  rvir = log10(4.0f/3.0f*PI_CUDA*DELTAV*RHOCRIT);
  rvir = (lgm - rvir)/3.0f;
  rvir = powf(10.0f,rvir);
  /***************************/

  /** ajuste de JS **/
  Deltae  = (5.0f*DELTAV);
  Deltae *= powf(1.0f/((bc*bc)*ab),0.75f);

  ce = FACTOR_JS*cvir;

  dc  = (Deltae*ce*ce*ce/3.0f);
  dc /= (log(1.0f + ce) - ce/(1.0f + ce));

  rs = rvir/cvir;
  /******************/

  r = x[2]*x[2] + x[1]*x[1]/(bc*bc) + x[0]*x[0]/(ab*ab*bc*bc);
  r = sqrtf(r);
  if(r < 1.e-7) r = 1.e-7;
  r /= rs;

  rho  = (dc*RHOCRIT);
  rho /= r;
  rho /= ((r + 1.0f)*(r + 1.0f));

  return rho;
}
__host__ __device__ float profile_nfw(float *x, float lgm)
{
  float cvir,rvir,dc,rs,r;
  float rho;

  /** parametro de concentracion **/
  cvir = 1.020f - 0.109f*(lgm - 12.0f); 
  cvir = exp10f(cvir);
  /********************************/

  /** radio virial esferico **/
  rvir = log10(4.0f/3.0f*PI_CUDA*DELTAV*RHOCRIT);
  rvir = (lgm - rvir)/3.0f;
  rvir = exp10f(rvir);
  /***************************/

  /** ajuste de JS **/
  dc  = (DELTAV*cvir*cvir*cvir/3.0f);
  dc /= (log(1.0f + cvir) - cvir/(1.0f + cvir));

  rs = rvir/cvir;
  /******************/

  r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  r = sqrtf(r);
  if(r < 1.e-4) r = 1.e-4;
  r /= rs;

  rho  = (dc*RHOCRIT);
  rho /= r;
  rho /= ((r + 1.0f)*(r + 1.0f));

  return rho;
}

__host__ __device__ float projected_profile(float *x, float lgm, float ab, float bc, 
                                   float costheta, float phi)
{
  float f,A,B,C;
  float st2,ct2,sp2,cp2;
  float st,ct,sp,cp;
  float cvir,rvir,ce,Deltae,dc,rs,r;
  float rho;
#ifndef DOSH
  float SQRTABC,qx2,qy2;
#endif

  if(bc < 0.1f) bc = 0.1f;
  if(ab < 0.1f) ab = 0.1f;

  /** parametro de concentracion **/
  cvir = 1.020f - 0.109f*(lgm - 12.0f); 
  cvir = exp10f(cvir);
  /********************************/

  /** radio virial esferico **/
  rvir = log10(4.0f/3.0f*PI_CUDA*DELTAV*RHOCRIT);
  rvir = (lgm - rvir)/3.0f;
  rvir = exp10f(rvir);
  /***************************/

  /** ajuste de JS **/
  Deltae  = (5.0f*DELTAV);
  Deltae *= powf(1.0f/((bc*bc)*ab),0.75f);

  ce = FACTOR_JS*cvir;
  dc  = (Deltae*ce*ce*ce/3.0f);
  dc /= (log(1.0f + ce) - ce/(1.0f + ce));

  rs = rvir/cvir;
  /******************/

  /** A, B, C **/
  ct = costheta;
  st = sqrtf(1.0f - ct*ct);
  sp = sin(phi);
  cp = cos(phi);
  
  st2 = st*st;
  ct2 = ct*ct;
  sp2 = sp*sp;
  cp2 = cp*cp;

  f = st2*(cp2/(ab*ab) + sp2)/(bc*bc) + ct2;
  
  A = ct2*(sp2/(ab*ab) + cp2)/(bc*bc) + st2/(bc*bc)/(bc*bc)/(ab*ab);
  B = ct*sin(2.0f*phi)*(1.0f/(ab*ab) - 1.0f)/(bc*bc);
  C = (sp2 + cp2/(ab*ab))/(bc*bc);
  /************/

#ifdef DOSH
  r = (A*x[0]*x[0]+B*x[0]*x[1]+C*x[1]*x[1])/f;
#else
  /** Q's **/
  SQRTABC = sqrt((A-C)*(A-C)+B*B);
  qx2 = 2.0f*f/(A + C - SQRTABC);
  qy2 = 2.0f*f/(A + C + SQRTABC);
  /*********/
  r = x[0]*x[0]/qx2 + x[1]*x[1]/qy2;
#endif
  r = sqrtf(r);
  if(r < 1.e-4) r = 1.e-4;
  r /= rs;

  rho  = (2.0f*dc*RHOCRIT);
  rho /= sqrtf(f);
  rho *= f_projected(r);

  return rho;
}

__host__ __device__ void cambio_de_coordenadas(float ab, float bc, float costheta, float phi, 
                                      float *x, float *y, float *z)
{
  float psi,cp,sp;
  float ct,st,cf,sf,st2,ct2,sf2,cf2;
  float A,B,C;
  float xt, yt, zt;

  ct = costheta;
  st = sqrtf(1.0f - ct*ct);
  cf = cosf(phi);
  sf = sinf(phi);

  st2 = st*st;
  ct2 = ct*ct;
  sf2 = sf*sf;
  cf2 = cf*cf;

  A = ct2*(sf2/(ab*ab)/(bc*bc) + cf2/(bc*bc)) + st2/(ab*ab*bc*bc*bc*bc);
  B = ct*sinf(2.0f*phi)*(1.0f/(ab*ab)/(bc*bc) - 1.0f/(bc*bc));
  C = (sf2/(bc*bc) + cf2/(ab*ab)/(bc*bc));
  /************/

  psi = 0.5f*atan(B/(A-C));

  if(((A - C)*cosf(2.0f*psi) + B*sinf(2.0f*psi)) > 0.0f) psi = psi - PI_CUDA*0.5f;

  cp = cos(psi);
  sp = sin(psi);
  
  xt = *x;
  yt = *y;
  zt = *z;

  *x = (-sf*cp - cf*ct*sp)*xt + (sf*sp - cf*ct*cp)*yt + (cf*st)*zt;
  *y = (cf*cp - sf*ct*sp)*xt + (-cf*sp - sf*ct*cp)*yt + (sf*st)*zt;
  *z = (st*sp)*xt + (st*cp)*yt + ct*zt;
}

__device__ float T1h(float *r, float *x)
{
  float tmp;
  float lgm      = x[0];
  float ab       = x[1];
  float bc       = x[2];
  float costheta = x[3];
  float phi      = x[4];
  float m,nu;
  float xp, yp, zp;

  nu = Nu_M(lgm);
  nu = exp10f(nu);
  m  = exp10f(lgm);

  xp = r[0];
  yp = r[1];
  zp = r[2];

  cambio_de_coordenadas(ab,bc,costheta,phi,&xp,&yp,&zp);

  r[0] = xp;
  r[1] = yp;  
  r[2] = zp;

  tmp  = profile_1h(r,lgm,ab,bc)/m;
#ifdef CG
	tmp *= RHOMEDIO/m;
  tmp *= (ncen(lgm) + nsat(lgm));
#endif
  tmp *= f_nu(nu)*nu*dlogNu_dlogM(lgm)*FMNORM;
  tmp *= forma(bc,ab);
  tmp *= LOGE10; //integracion logaritmica en m1
  return(tmp);
}

__device__ float T2h(float *r, float *x){
  float tmp;
  float lgm1, lgm2;
  float m1, m2;
  float nu1, nu2;
  float bc1, bc2;
  float ab1, ab2;
  float costheta1, phi1;
  //float costheta2, phi2;
  float xc, yc, zc;
  float y[3];

  /***********
  x[0]       = M1
  x[2,1]     = bc1,ab1 //bugy
  x[3,4]     = e11,e12
  x[5]       = M2
  x[6,7]     = bc2,ab2 //bugy
  x[8,9]     = e21,e22
  x[10,11,12]   = x,y,z
  **********/
  
  lgm1      = x[0];
  ab1       = x[1];
  bc1       = x[2];
  costheta1 = x[3];
  phi1      = x[4];

  lgm2      = x[5];
  ab2       = x[6];
  bc2       = x[7];
  //costheta2 = x[8];
  //phi2      = x[9];
  y[0]      = x[10];
  y[1]      = x[11];
  y[2]      = x[12];

  m1  = exp10f(lgm1);
  nu1 = Nu_M(lgm1);
  nu1 = exp10f(nu1);
  m2  = exp10f(lgm2);
  nu2 = Nu_M(lgm2);
  nu2 = exp10f(nu2);

  xc = r[0];
  yc = r[1];
  zc = r[2];

  cambio_de_coordenadas(ab1,bc1,costheta1,phi1,&xc,&yc,&zc);

  xc += y[0];
  yc += y[1];
  zc += y[2];

  tmp  = profile_1h(y,lgm2,1.0f,1.0f)/m2;
	tmp *= RHOMEDIO/m1;
  #ifdef CG
  tmp *= RHOMEDIO/m2;
  tmp *= (ncen(lgm2) + nsat(lgm2));
  #endif
  tmp *= f_nu(nu2)*nu2*dlogNu_dlogM(lgm2)*FMNORM;
  tmp *= forma(bc2,ab2);
  tmp *= f_nu(nu1)*nu1*dlogNu_dlogM(lgm1)*FMNORM;
  tmp *= forma(bc1,ab1);
  tmp *= bias(nu1);
  tmp *= bias(nu2);
  tmp *= xilin(xc,yc,zc);
  tmp *= align_masa(xc,yc,zc);
  tmp *= LOGE10;  //integracion logaritmica en m2
  tmp *= LOGE10;  //integracion logaritmica en m1

  return(tmp);
}

/*f(nu) para Sheth & Tormen. Ec. (4) Seljak 2000. Sin normalizar.*/
__device__ float f_nu_log(float x){
  float nu ;
  nu = powf(10.0f,x);
  return(f_nu(nu));
}

__device__ float bias_nu_log(float x){
  float nu ;
  nu = powf(10.0f,x);
  return(bias(nu));
}

/*Integrador para calcular el nc_medio*/
__global__ void integra_nc(float Numin, float Numax, curandState *state){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int bt = blockIdx.x;
  const unsigned int it = threadIdx.x;
  int i;
  float x;
  float value, sigma, tmp;
  __shared__ float s_value[THREADS_PER_BLOCK];
  __shared__ float s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState xx = state[tid];

  /*Tira un numero random x entre NuMin y NuMax*/
  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    x = (Numax - Numin)*curand_uniform(&xx) + Numin;

    tmp = f_nu(powf(10.0f,x))*powf(10.0f,x)/powf(10.0f,M_Nu(x));

    value += tmp;
    sigma += tmp*tmp;
  }

  /*Suma a la memoria global*/
  s_value[it] = value;
  s_sigma[it] = sigma;
  __syncthreads();

  if(it == 0){
    value = 0.0f; sigma = 0.0f;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }

    /*Suma a la memoria global*/
    d_int[bt] = value;
    d_sig[bt] = sigma;
  }
}

__inline__ float nc_medio(float Numin, float Numax, dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int k;
  float r, s;
  float volumen;

  integra_nc<<<dimGrid,dimBlock>>>(Numin,Numax,devStates);
  cudaThreadSynchronize();

  volumen = Numax - Numin;

  /*Copia sumatorias al host*/
  cudaMemcpyFromSymbol(h_int,d_int,RNGS/THREADS_PER_BLOCK*sizeof(float));
  cudaMemcpyFromSymbol(h_sig,d_sig,RNGS/THREADS_PER_BLOCK*sizeof(float));

  /*Termina de reducir en el host*/
  r = 0.0f; s = 0.0f;
  for(k = 0; k < RNGS/THREADS_PER_BLOCK; k++){
    r += h_int[k];
    s += h_sig[k];
  }

  /*Estima la integral y el sigma*/
  r /= (float)((long)RNGS*(long)LAZOSPLUS);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s -= (r*r);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s  = sqrt(s);

  r *= volumen;
  s *= volumen;

  return(r*FMNORM*RHOMEDIO*LOGE10);
}

/*Integrador de la funcion de masa para calcular la normalizacion*/
__global__ void integra_f_nu_log(curandState *state){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int it = threadIdx.x;
  const unsigned int bt = blockIdx.x;
  int    i;
  float value;
  float sigma;
  float tmp;
  __shared__ float s_value[THREADS_PER_BLOCK];
  __shared__ float s_sigma[THREADS_PER_BLOCK];
  float  x;

  /*Estado del RNG*/
  curandState xx = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    /*Tira un numero random x entre NuMin y NuMax*/
    x = (NuMax - NuMin)*(float)curand_uniform(&xx) + NuMin;
    tmp = powf(10.0f,x)*f_nu_log(x)*LOGE10;

    value += tmp;
    sigma += tmp*tmp;
  }
  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();
  if(it == 0){
    value = 0.0f; sigma = 0.0f;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }

    /*Suma a la memoria global*/
    d_int[bt] = value;
    d_sig[bt] = sigma;
  }
  
  /*Guarda el estado del RNG*/
  state[tid] = xx;
}

/*Computa la normalizacion de la funcion de masa*/
__inline__ void normalizacion_func_masa(dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int k;
  float r, s;
  float volumen;

  integra_f_nu_log<<<dimGrid,dimBlock>>>(devStates);
  cudaThreadSynchronize();

  volumen = NuMax - NuMin;

  /*Copia sumatorias al host*/
  HANDLE_ERROR(cudaMemcpyFromSymbol(h_int,d_int,RNGS/THREADS_PER_BLOCK*sizeof(float)));
  HANDLE_ERROR(cudaMemcpyFromSymbol(h_sig,d_sig,RNGS/THREADS_PER_BLOCK*sizeof(float)));

  /*Termina de reducir en el host*/
  r = 0.0f; s = 0.0f;
  for(k = 0; k < RNGS/THREADS_PER_BLOCK ; k++){
    r += h_int[k];
    s += h_sig[k];
  }

  /*Estima la integral y el sigma*/
  r /= (float)((long)RNGS*(long)LAZOSPLUS);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s -= (r*r);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s  = sqrtf(s);
  s *= volumen;
  r *= volumen;

  norma_funmasa = 1.0f/r;

  /*Imprime en file de salida*/
  fprintf(stdout,"Normalizacion de la Funcion de Masa: %e +/- %e\n",norma_funmasa,fabs(1.0/r/r*s));
}

/*Kernel: que inicializa los RNGs*/
__global__ void setup_kernel(curandState *state, int entero){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  curand_init(entero, tid, 0, &state[tid]);
}

/*Kernel: Inicializa los contadores a cero*/
__global__ void init_counters(void){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  d_int[tid] = 0.0f;
  d_sig[tid] = 0.0f;
}

/*Integrador de la funcion de masa para calcular la normalizacion*/
__global__ void integra_bias_f(curandState *state){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int it = threadIdx.x;
  const unsigned int bt = blockIdx.x;

  float value, sigma, tmp;
  float x;
  int i;
  __shared__ float s_value[THREADS_PER_BLOCK], s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState xx = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    /*Tira un numero random x entre NuMin y NuMax*/
    x = (NuMax - NuMin)*curand_uniform(&xx) + NuMin;
  
    tmp = powf(10.0f,x)*bias_nu_log(x)*f_nu_log(x)*FMNORM*LOGE10;
    value += tmp;
    sigma += tmp*tmp;
  }
  
  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();

  if(it == 0){
    value = 0.0f; sigma = 0.0f;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }

    /*Suma a la memoria global*/
    d_int[bt] = value;
    d_sig[bt] = sigma;
  }

  /*Guarda el estado del RNG*/
  state[tid] = xx;
}

/*Computa la normalizacion de la funcion de masa*/
__inline__ void prueba_bias_f(dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int k;
  float r, s;
  float volumen;

  volumen = NuMax - NuMin;

  integra_bias_f<<<dimGrid,dimBlock>>>(devStates);
  cudaThreadSynchronize();

  /*Copia sumatorias al host*/
  cudaMemcpyFromSymbol(h_int,d_int,RNGS/THREADS_PER_BLOCK*sizeof(float));
  cudaMemcpyFromSymbol(h_sig,d_sig,RNGS/THREADS_PER_BLOCK*sizeof(float));

  /*Termina de reducir en el host*/
  r = 0.0f; s = 0.0f;
  for(k = 0; k < RNGS/THREADS_PER_BLOCK ; k++){
    r += h_int[k];
    s += h_sig[k];
  }

  /*Estima la integral y el sigma*/
  r /= (float)((long)RNGS*(long)LAZOSPLUS);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s -= (r*r);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s  = sqrt(s);
  s *= volumen;
  r *= volumen;

  /*Imprime en file de salida*/
  fprintf(stdout,"Prueba Bias*F: %e +/- %e\n",r,s);
}

/*Integrador de la funcion de masa para calcular la normalizacion*/
__global__ void kernel_integra_merchan(curandState *state)
{
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int it = threadIdx.x;
  const unsigned int bt = blockIdx.x;

  float value, sigma, tmp;
  float x,y;
  int i;
  __shared__ float s_value[THREADS_PER_BLOCK], s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState xx = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    x = curand_uniform(&xx);
    y = curand_uniform(&xx);
  
    tmp = forma(x,y);
    value += tmp;
    sigma += tmp*tmp;
  }
  
  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();

  if(it == 0){
    value = 0.0f; sigma = 0.0f;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }

    /*Suma a la memoria global*/
    d_int[bt] = value;
    d_sig[bt] = sigma;
  }

  /*Guarda el estado del RNG*/
  state[tid] = xx;
}

/*Computa la normalizacion de la funcion de masa*/
__inline__ float integra_merchan(dim3 dimGrid, dim3 dimBlock, curandState *devStates)
{
  int k;
  float r, s;

  kernel_integra_merchan<<<dimGrid,dimBlock>>>(devStates);
  cudaThreadSynchronize();

  /*Copia sumatorias al host*/
  cudaMemcpyFromSymbol(h_int,d_int,RNGS/THREADS_PER_BLOCK*sizeof(float));
  cudaMemcpyFromSymbol(h_sig,d_sig,RNGS/THREADS_PER_BLOCK*sizeof(float));

  /*Termina de reducir en el host*/
  r = 0.0f; s = 0.0f;
  for(k = 0; k < RNGS/THREADS_PER_BLOCK ; k++){
    r += h_int[k];
    s += h_sig[k];
  }

  /*Estima la integral y el sigma*/
  r /= (float)((long)RNGS*(long)LAZOSPLUS);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s -= (r*r);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s  = sqrt(s);

  /*Imprime en file de salida*/
  fprintf(stdout,"Integra Merchan: %e +/- %e\n",r,s);

  return(r);
}

/*Integrador de la funcion de masa para calcular la normalizacion*/
__global__ void kernel_integra_align(curandState *state){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int it = threadIdx.x;
  const unsigned int bt = blockIdx.x;

  float value, sigma, tmp;
  float theta, phi;
  int i;
  __shared__ float s_value[THREADS_PER_BLOCK], s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState xx = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    phi = curand_uniform(&xx)*PI_CUDA*0.5f;
    theta = curand_uniform(&xx)*PI_CUDA*0.5f;
    tmp = align(theta,phi);

    value += tmp;
    sigma += tmp*tmp;
  }

  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();

  if(it == 0){
    value = 0.0f; sigma = 0.0f;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }

    /*Suma a la memoria global*/
    d_int[bt] = value;
    d_sig[bt] = sigma;
  }

  /*Guarda el estado del RNG*/
  state[tid] = xx;
}

/*Computa la normalizacion de la funcion de masa*/
__inline__ float integra_align(dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int k;
  float r, s;
  float volumen;

  kernel_integra_align<<<dimGrid,dimBlock>>>(devStates);
  cudaThreadSynchronize();

  /*Copia sumatorias al host*/
  cudaMemcpyFromSymbol(h_int,d_int,RNGS/THREADS_PER_BLOCK*sizeof(float));
  cudaMemcpyFromSymbol(h_sig,d_sig,RNGS/THREADS_PER_BLOCK*sizeof(float));

  /*Termina de reducir en el host*/
  r = 0.0f; s = 0.0f;
  for(k = 0; k < RNGS/THREADS_PER_BLOCK ; k++){
    r += h_int[k];
    s += h_sig[k];
  }

  /*Estima la integral y el sigma*/
  r /= (float)((long)RNGS*(long)LAZOSPLUS);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s -= (r*r);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s  = sqrt(s);

  volumen = (float)M_PI*0.5f*(float)M_PI*0.5f;

  r *= volumen;
  s *= volumen;

  /*Imprime en file de salida*/
  fprintf(stdout,"Integra Align: %e +/- %e\n",r,s);

  return(r);
}

/*Integrador de la funcion de masa para calcular la normalizacion*/
__global__ void integra_perfil(float rmin, float rmax, float m, curandState *state){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int it = threadIdx.x;
  const unsigned int bt = blockIdx.x;

  float value, sigma, tmp;
  int i;
  __shared__ float s_value[THREADS_PER_BLOCK], s_sigma[THREADS_PER_BLOCK];

  float r[3];
  float rr, phi, costheta;

  /*Estado del RNG*/
  curandState xx = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    rr       = curand_uniform(&xx)*(rmax - rmin) + rmin;
    rr       = powf(10.0f,rr);
    phi      = PI_CUDA*0.5f*curand_uniform(&xx);
    costheta = curand_uniform(&xx);

    r[0] = rr*sqrt(1.0f - costheta*costheta)*cos(phi);
    r[1] = rr*sqrt(1.0f - costheta*costheta)*sin(phi);
    r[2] = rr*costheta;

    tmp = profile_1h(r,m,1.0f,1.0f)*powf(rr,3.0f)/powf(10.0f,m)*LOGE10;

    value += tmp;
    sigma += tmp*tmp;
  }

  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();

  if(it == 0){
    value = 0.0f; sigma = 0.0f;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }

    /*Suma a la memoria global*/
    d_int[bt] = value;
    d_sig[bt] = sigma;
  }

  /*Guarda el estado del RNG*/
  state[tid] = xx;
}

/*Computa la normalizacion de la funcion de masa*/
__inline__ void normalizacion_perfil(dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int k;
  float r, s;
  float lgm;
  float rvir;
  float volumen;
  float rmin, rmax;

  lgm = (CENTROS_MASA_MAX + CENTROS_MASA_MIN)*0.5f;

  /*** radio virial esferico ***/
  rvir = log10(4.0f/3.0f*PI_CUDA*DELTAV*RHOCRIT);
  rvir = (lgm - rvir)/3.0f;

  rmin = -4.0f;
  rmax = rvir;

  rvir = powf(10.0f,rvir);
  printf("Log10(masa): %f Rvir: %f\n",lgm,rvir);

  integra_perfil<<<dimGrid,dimBlock>>>(rmin,rmax,lgm,devStates);
  cudaThreadSynchronize();

  /*Copia sumatorias al host*/
  cudaMemcpyFromSymbol(h_int,d_int,RNGS/THREADS_PER_BLOCK*sizeof(float));
  cudaMemcpyFromSymbol(h_sig,d_sig,RNGS/THREADS_PER_BLOCK*sizeof(float));

  /*Termina de reducir en el host*/
  r = 0.0f; s = 0.0f;
  for(k = 0; k < RNGS/THREADS_PER_BLOCK ; k++){
    r += h_int[k];
    s += h_sig[k];
  }

  /*Estima la integral y el sigma*/
  r /= (float)((long)RNGS*(long)LAZOSPLUS);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s -= (r*r);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s  = sqrt(s);

  volumen = 0.5f*M_PI*(rmax-rmin)*8.0f;

  r *= volumen;
  s *= volumen;

  /*Imprime en file de salida*/
  fprintf(stdout,"Normalizacion perfil: %e +/- %e\n",r,s);
}

/*Integrador para calcular el ng_medio*/
__global__ void integra_ng(float Numin, float Numax, curandState *state){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int bt = blockIdx.x;
  const unsigned int it = threadIdx.x;
  int i;
  float x;
  float value, sigma, tmp;
  __shared__ float s_value[THREADS_PER_BLOCK];
  __shared__ float s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState xx = state[tid];

  /*Tira un numero random x entre NuMin y NuMax*/
  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    x = (Numax - Numin)*curand_uniform(&xx) + Numin;

    tmp = M_Nu(x);
    tmp = f_nu(powf(10.0f,x))*powf(10.0f,x)/powf(10.0f,tmp)*(ncen(tmp) + nsat(tmp));

    value += tmp;
    sigma += tmp*tmp;
  }

  /*Suma a la memoria global*/
  s_value[it] = value;
  s_sigma[it] = sigma;
  __syncthreads();

  if(it == 0){
    value = 0.0f; sigma = 0.0f;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }

    /*Suma a la memoria global*/
    d_int[bt] = value;
    d_sig[bt] = sigma;
  }
}

__inline__ float ng_medio(float Numin, float Numax, dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int k;
  float r, s;
  float volumen;

  integra_ng<<<dimGrid,dimBlock>>>(Numin,Numax,devStates);
  cudaThreadSynchronize();

  volumen = Numax - Numin;

  /*Copia sumatorias al host*/
  cudaMemcpyFromSymbol(h_int,d_int,RNGS/THREADS_PER_BLOCK*sizeof(float));
  cudaMemcpyFromSymbol(h_sig,d_sig,RNGS/THREADS_PER_BLOCK*sizeof(float));

  /*Termina de reducir en el host*/
  r = 0.0f; s = 0.0f;
  for(k = 0; k < RNGS/THREADS_PER_BLOCK; k++){
    r += h_int[k];
    s += h_sig[k];
  }

  /*Estima la integral y el sigma*/
  r /= (float)((long)RNGS*(long)LAZOSPLUS);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s -= (r*r);
  s /= (float)((long)RNGS*(long)LAZOSPLUS);
  s  = sqrt(s);

  r *= volumen;
  s *= volumen;

  return(r*FMNORM*RHOMEDIO*LOGE10);
}
