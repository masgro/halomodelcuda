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

__device__ void rotate(float *y,float phi,float costheta,float psi,float *x){
  float c[3][3];
  float sintheta = sqrt(1.0f - costheta*costheta);

  /*Define la matriz de rotacion*/
  c[0][0] =  cos(psi)*cos(phi) - costheta*sin(phi)*sin(psi);
  c[0][1] =  cos(psi)*sin(phi) + costheta*cos(phi)*sin(psi);
  c[0][2] =  sin(psi)*sintheta;
  c[1][0] = -sin(psi)*cos(phi) - costheta*sin(phi)*cos(psi);
  c[1][1] = -sin(psi)*sin(phi) + costheta*cos(phi)*cos(psi);
  c[1][2] =  cos(psi)*sintheta;
  c[2][0] =  sintheta*sin(phi);
  c[2][1] = -sintheta*cos(phi);
  c[2][2] =  costheta;
  
  /*Computa la rotacion*/
  x[0] = c[0][0]*y[0] + c[0][1]*y[1] + c[0][2]*y[2];
  x[1] = c[1][0]*y[0] + c[1][1]*y[1] + c[1][2]*y[2];
  x[2] = c[2][0]*y[0] + c[2][1]*y[1] + c[2][2]*y[2];
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
  fi  = exp(-(f-m1)*(f-m1)*0.5f/(m2*m2));
  fi /= SQRT_TWOPI_CUDA;
  fi /= m2;
  return(fi);
}

__device__ float forma(float bc, float ab){
  float prob;
  float dbc,dab;
  float m1bc,m2bc;
  float m1ab,m2ab;

  m1bc = BCMEDIO; m2bc = 0.1f; //sigma
  //m1bc = d_bcmedio; m2bc = 0.1f; //sigma
  dbc  = db_gauss(bc,m1bc,m2bc)*(1.0f-bc)*bc;

  m1ab = ABMEDIO; m2ab = 0.1f;  //sigma
  //m1ab = d_abmedio; m2ab = 0.1f;  //sigma
  dab  = db_gauss(ab,m1ab,m2ab)*(1.0f-ab)*ab;

  prob = dab*dbc;
  if(prob < 0.0f) prob = 0.0f;

  return(prob);
}

  /*perfil de jing y suto a z=0* a[0]>a[1]>a[2]**/
__host__ __device__ float u(float *r, float bc, float ab, float lgm, float rlim){
  float cvir,rvir,ce,Deltae,dc,rs,rr,rho;
  float prol;

  if(bc < 0.01f) bc = 0.01f;
  if(ab < 0.01f) ab = 0.01f;

  /** parametro de concentracion **/
  cvir = 1.020f - 0.109f*(lgm - 12.0f); 
  cvir = powf(10.0f,cvir);
//log c = 1.020(+/- 0.015) - 0.109(+/- 0.005)*(log Mvir - 12.0); 
//http://arxiv.org/pdf/astro-ph/0608157v2
  /********************************/

  /** radio virial esferico **/
  rvir = log10(4.0f/3.0f*PI_CUDA*DELTAV*RHOCRIT);
  rvir = (lgm - rvir)/3.0f;
  rvir = powf(10.0f,rvir);
#ifdef DOSH
  //if(rvir > rlim)return(0.0f);
#endif
  /***************************/

  /** ajuste de JS **/
  ce = FACTOR_JS*cvir;

  /** como en JS ecuacion 13  smith 59 **/
  prol = 1.0f/((bc*bc)*ab);
  Deltae = 5.0f*DELTAV*powf(prol ,0.75f);
  dc = (Deltae*ce*ce*ce/3.0f)/(log(1.0f+ce)-ce/(1.0f+ce));

  /** ecuacion 59 Smith **/
  rs = rvir/cvir;
  /******************/

  /** ecuacion 4 paper smith **/
  rr = r[2]*r[2] + r[1]*r[1]/(bc*bc) + r[0]*r[0]/(ab*ab*bc*bc);
  rr = sqrtf(rr);

#ifdef DOSH
  if(rr > 1.5f*rvir) return(0.0f);
#else
  //if(rr > 3.0f*rvir) return(0.0f);
#endif

  rr /= rs;
  rho  = RHOCRIT/powf(10.0f,lgm);
  rho *= dc;
  rho /= rr;
  rho /= ((rr + 1.0f)*(rr + 1.0f));

  if(rho < 1.E-33) rho = 0.0f;
  return rho;
}

__host__ __device__ float u_esferico(float *r, float lgm, float rlim){
  float cvir,rvir,dc,rs,rr,y,rho;

  cvir = 1.020f - 0.109f*(lgm - 12.0f); 
  cvir = powf(10.0f,cvir);

/*** radio virial esferico ***/
  rvir = log10(4.0f/3.0f*PI_CUDA*DELTAV*RHOCRIT);
  rvir = (lgm - rvir)/3.0f;
  rvir = powf(10.0f,rvir);
#ifdef DOSH
  //if(rvir > rlim*0.5f)return(0.0f);
#endif
/*****************************/

  dc = (DELTAV*cvir*cvir*cvir/3.0f)/(log(1.0f+cvir)-cvir/(1.0f+cvir));

/*** ecuacion 59 Smith ***/
  rs = rvir/cvir;

/*ecuacion 4 paper smith*/
  rr = r[2]*r[2] + r[1]*r[1] + r[0]*r[0];
  rr = sqrtf(rr);

  if(rr > 2.0f*rvir) return(0.0f);
  
  y = rr/rs;
  rho  = RHOCRIT/powf(10.0f,lgm);
  rho *= dc;
  rho /= y;
  rho /= ((y + 1.0f)*(y + 1.0f));

  if(rho < 1.E-33) rho = 0.0f;
  return rho;
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

__host__ __device__ float n(float lgm){
  //float coeff[6] = 
  //{ 6.69172970e+02,
  //  -2.88844180e+02,
  //  4.96380268e+01,
  //  -4.28038461e+00,
  //  1.84104417e-01,
  //  -3.16125044e-03};
  //  {2.24744465e+03,
  //  -9.25455609e+02,
  //  1.52058721e+02,
  //  -1.24959391e+01,
  //  5.12676413e-01,
  //  -8.40294028e-03};

  float coeff[5] = 
  { 5.69393747e+02,  
   -1.72677633e+02,   
    1.92362083e+01,
   -9.26525034e-01,   
   1.62235706e-02};

  float q;
  q = polim(lgm,coeff,4);
  q = powf(10.0f,q-20.75321);
  return(q);
}

__device__ float align(float theta, float phi){
  float p;
  p = ALIGN_C*expf(-(theta - PI_CUDA*0.5f)*(theta - PI_CUDA*0.5f)*0.5f/(ALIGN_B*ALIGN_B));
  p *= expf(-phi*phi*0.5f/(ALIGN_B*ALIGN_B));
  p += expf(-theta*theta*0.5f/(ALIGN_B*ALIGN_B));
  return(p);
}

__device__ float b_r(float r){
  float lr = log10f(r);
  float tmp;
  tmp = (lr > 0.0f)? -0.08333*lr + 0.4f : 0.4f;
  if(lr >= 3.0f) tmp = 0.0f;
  return(tmp);
}

__device__ float c_r(float r){
  float lr = log10f(r);
  float tmp;
  float a = 0.0438600f ;
  float b = 0.1494210f ;
  float c = 1.0801300f ;
  tmp = (lr > 0.0f)? 0.05f*a*lr*lr + b*lr + c : c;
  return(tmp);
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

  /*TINKER*/
  //q = 0.5f;
  //q *= (1.0f + (float)erf((lgm - 11.5f)/NCEN_SIGMA));

  return(q);
}

__host__ __device__ float nsat(float lgm){
  float q;
  q = (powf(10.0f,lgm) - NSAT_LOGM0);
  q /= NSAT_LOGM1;
  q = (q > 0.0f)? powf(q,NSAT_ALPHA) : 0.0f;

  /*Parametrizacion 4-param*/
  //q = powf(10.0f,lgm);
  //q /= powf(10.0f,NSAT_LOGM0);
  //q = (q >= 0.0f)? powf(q,NSAT_ALPHA) : 0.0f;

  /*Parametrizacion Tinker 2005*/
  //if(lgm < 11.5f){
  //  return(0.0f);
  //} else {
  //  q = powf(10.0f,12.7f)/(powf(10.0f,lgm)-powf(10.0f,11.5f));
  //  q = expf(-1.0f*q);
  //  q *= powf(10.0f,lgm-12.8f);
  //}

  return(q);
}

__host__ __device__ float mom1(float lgm){
  float q,qcen,qsat;
  qcen = ncen(lgm);
  qsat = nsat(lgm);
  if(lgm >= NCEN_LOGMMIN)
    q = qcen + qsat;
  else
    q = qcen;

  //q = nsat(lgm);
  return(q);
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
  curandState seed = state[tid];

  /*Tira un numero random x entre NuMin y NuMax*/
  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    x = (Numax - Numin)*curand_uniform(&seed) + Numin;

    tmp = f_nu(powf(10.0f,x))*powf(10.0f,x)/powf(10.0f,M_Nu(x));

    /*FIT*/
    //tmp = n(x)*powf(10.0,x);

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

  /*Guarda el estado del RNG*/
  state[tid] = seed;
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
  /*FIT*/
  //return(r*LOGE10);
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
  curandState seed = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    /*Tira un numero random x entre NuMin y NuMax*/
    x = (NuMax - NuMin)*(float)curand_uniform(&seed) + NuMin;
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
  state[tid] = seed;
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
  curandState seed = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    /*Tira un numero random x entre NuMin y NuMax*/
    x = (NuMax - NuMin)*curand_uniform(&seed) + NuMin;
  
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
  state[tid] = seed;
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

  double value, sigma;
  float tmp;
  float x,y;
  int i;
  __shared__ double s_value[THREADS_PER_BLOCK], s_sigma[THREADS_PER_BLOCK];

  const float fmin = 0.0f, fmax = 1.0f;

  /*Estado del RNG*/
  curandState seed = state[tid];

  value = 0.0; sigma = 0.0;
  for(i = 0; i < LAZOSPLUS; i++){

    x = (fmax-fmin)*curand_uniform(&seed)+fmin;
    y = (fmax-fmin)*curand_uniform(&seed)+fmin;
  
    tmp = forma(x,y);

    value += tmp;
    sigma += tmp*tmp;
  }
  
  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();

  if(it == 0){
    value = 0.0; sigma = 0.0;
    for(i = 0; i < THREADS_PER_BLOCK; i++){
      value += s_value[i];
      sigma += s_sigma[i];
    }

    /*Suma a la memoria global*/
    d_int[bt] = value;
    d_sig[bt] = sigma;
  }

  /*Guarda el estado del RNG*/
  state[tid] = seed;
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

  //float volumen = (1.0f-0.1f)*(1.0f-0.1f);

  //r *= volumen;
  //s *= volumen;

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
  float x,y,z;
  float costheta, phi;
  int i;
  __shared__ float s_value[THREADS_PER_BLOCK], s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState seed = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    phi = curand_uniform(&seed)*PI_CUDA*0.5f;
    costheta = curand_uniform(&seed);
    x = costheta;
    y = sqrt(1.0f - costheta*costheta)*cos(phi);
    z = sqrt(1.0f - costheta*costheta)*sin(phi);
    tmp = align_masa(x,y,z);

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
  state[tid] = seed;
}

/*Computa la normalizacion de la funcion de masa*/
__inline__ float integra_align(dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int k;
  float r, s;

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

  //volumen = (float)M_PI*0.5f*(float)M_PI*0.5f;

  //s *= volumen;
  //r *= volumen;

  /*Imprime en file de salida*/
  fprintf(stdout,"Integra Align: %e +/- %e\n",r,s);
  fflush(stdout);
  return(r);
}

/*Integrador de la funcion de masa para calcular la normalizacion*/
__global__ void integra_perfil(float rmin, float rmax, float lgm, curandState *state){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int it = threadIdx.x;
  const unsigned int bt = blockIdx.x;

  float value, sigma, tmp;
  float r[3],rr,phi,costheta,sintheta;
  int i;
  __shared__ float s_value[THREADS_PER_BLOCK], s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState seed = state[tid];

  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){

    rr       = curand_uniform(&seed)*(rmax - rmin) + rmin;
    rr       = powf(10.0f,rr);
    phi      = PI_CUDA*0.5f*curand_uniform(&seed);
    costheta = curand_uniform(&seed);
		sintheta = sqrtf(1.0f - costheta*costheta);

    r[0] = rr*sintheta*cosf(phi);
    r[1] = rr*sintheta*sinf(phi);
    r[2] = rr*costheta;

    tmp = u(r,BCMEDIO,ABMEDIO,lgm,1.0f)*powf(rr,3.0f)*LOGE10;

    //r[0] = (curand_uniform(&seed)*(rmax - rmin) + rmin);
    //r[0] = powf(10.0f,r[0]);
    //tmp = r[0]*LOGE10;
    //r[0] /= (ABMEDIO*BCMEDIO);

    //r[1] = (curand_uniform(&seed)*(rmax - rmin) + rmin);
    //r[1] = powf(10.0f,r[1]);
    //tmp *= r[1]*LOGE10;
    //r[1] /= BCMEDIO;

    //r[2] = (curand_uniform(&seed)*(rmax - rmin) + rmin);
    //r[2] = powf(10.0f,r[2]);
    //tmp *= r[2]*LOGE10;

    //tmp *= u(r,BCMEDIO,ABMEDIO,lgm,1.0f);

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
  state[tid] = seed;
}

/*Computa la normalizacion de la funcion de masa*/
__inline__ float normalizacion_perfil(dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int k;
  float r, s;
  float lgm = 15.0f;
  float volumen;
  float rvir,rmin,rmax;

  lgm = (CENTROS_MASA_MAX + CENTROS_MASA_MIN)*0.5f;

  /*** radio virial esferico ***/
  rvir = log10(4.0f/3.0f*PI_CUDA*DELTAV*RHOCRIT);
  rvir = (lgm - rvir)/3.0f;

	/*{
	FILE *pfout;
	float rr[3],phi,costheta,sintheta;
	pfout = fopen("salida.dat","w");
	for(k = 0; k < 1000; k++){
		r = (RMAX - RMIN)*(float)k/1000. + RMIN;
    r = pow(10.0f,r);
    phi      = 2.0f*PI_CUDA*drand48();
    costheta = 2.0f*drand48()-1.0f;
		sintheta = sqrtf(1.0f - costheta*costheta);

    rr[0] = r*sintheta*cosf(phi);
    rr[1] = r*sintheta*sinf(phi);
    rr[2] = r*costheta;

		fprintf(pfout,"%f %f\n",r,u(rr,1.0f,1.0f,lgm,100.0f));
	}
	fclose(pfout);
	}*/

  rmin = -4.0f;
  rmax = rvir - 0.2f;
  printf("rvir: %f rmin: %f rmax: %f\n",rvir,rmin,rmax);

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
  return(r);
}

/*Integrador para calcular el ng_medio*/
__global__ void integra_ng(float Numin, float Numax, curandState *state){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int bt = blockIdx.x;
  const unsigned int it = threadIdx.x;
  int i;
  float x;
  float value, sigma, tmp;
  float lgm, nu;
  float m;
  __shared__ float s_value[THREADS_PER_BLOCK];
  __shared__ float s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState seed = state[tid];

  /*Tira un numero random x entre NuMin y NuMax*/
  value = 0.0f; sigma = 0.0f;
  for(i = 0; i < LAZOSPLUS; i++){
    x = (Numax - Numin)*curand_uniform(&seed) + Numin;

    lgm = M_Nu(x);
    nu = powf(10.0f,x);
    m = powf(10.0f,lgm);
    tmp = f_nu(nu)*nu/m*mom1(lgm);

    /*FUNCION DE MASA FIT*/
    //m = powf(10.0f,x);
    //tmp = n(x)*mom1(x)*m;

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

  /*Guarda el estado del RNG*/
  state[tid] = seed;
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

  /*FUNCION DE MASA FIT*/
  //return(r*LOGE10);
}

float integrando_ng_medio(float m){
  float lgm,nu;
  lgm = log10(m);
  nu = pow(10.0,Nu_M(lgm));
  return(f_nu(nu)*nu/m*mom1(lgm)*dlogNu_dlogM(lgm)*RHOMEDIO*FMNORM*LOGE10);
  //return(n(lgm)*m*mom1(lgm)*LOGE10);
}

float massfunc(float lgm) {
  float q;
  float nu,m;
  nu = Nu_M(lgm);
  nu = pow(10.0,nu);
  m = pow(10.0,lgm);
  q = f_nu(nu)*RHOMEDIO*nu/m/m*dlogNu_dlogM(lgm)*FMNORM;

  //q = n(lgm);
  return(q);
}

__inline__ void prueba_ng_medio(){
  FILE *pfout;
  int k;
  float r, s;
  float Mmin = 8.0;
  float Mmax = 16.0;

  pfout = fopen("integrando_ng_medio.dat","w");
  for(k = 0; k < LAZOSPLUS; k++) {
    r = (Mmax - Mmin)/(float)LAZOSPLUS*(float)k + Mmin;
    s = integrando_ng_medio(pow(10.0,r));
    fprintf(pfout,"%e %e %e %e %e\n",r,s,mom1(r),ncen(r),nsat(r));
  }
  fclose(pfout);

  pfout = fopen("massfunc.dat","w");
  for(k = 0; k < LAZOSPLUS; k++) {
    r = (Mmax - Mmin)/(float)LAZOSPLUS*(float)k + Mmin;
    fprintf(pfout,"%.7e %.7e\n",r,massfunc(r));
  }
  fclose(pfout);

  pfout = fopen("tabla.dat","w");
  for(k = 0; k < LAZOSPLUS; k++) {
    r = (Mmax - Mmin)/(float)LAZOSPLUS*(float)k + Mmin;
    fprintf(pfout,"%.7e %.7e %.7e %.7e\n",r,massfunc(r),mom1(r),bias(powf(10.0f,Nu_M(r))));
  }
  fclose(pfout);
}
