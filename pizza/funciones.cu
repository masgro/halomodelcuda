__host__ __device__ float chingada(float lgm){
  float p1 = 0.80467523f*expf(-1.0f*(lgm-13.87587371f)*(lgm-13.87587371f)/(2.0f*0.45226793f*0.45226793f));
  float p2 = 0.24100481f*expf(-1.0f*(lgm-14.65655543f)*(lgm-14.65655543f)/(2.0f*0.15555877f*0.15555877f));
  if(lgm <= 12.23 || lgm >= 15.33)return(0.0);
  return(p1+p2);
}

__host__ __device__ float cambio_de_coordenadas(float ab, float bc, float costheta, 
                                               float phi, float *r){
  float psi,cp,sp;
  float ct,st,cf,sf,st2,ct2,sf2,cf2;
  float xt,yt,zt;
  //float theta;
  float f;

  //theta = acosf(costheta);
  ct = costheta;
  ct2 = ct*ct;
  st2 = 1.0f-ct2;
  st = sqrt(st2);

  cf = cosf(phi);
  sf = sinf(phi);
  sf2 = sf*sf;
  cf2 = cf*cf;

  ///************/
  //float A,B,C;
  //float ca2,cb2;
  //ca2 = 1./(ab*ab*bc*bc);
  //cb2 = 1./(bc*bc);

  //A = ct2*(ca2*sf2 + cb2*cf2) + ca2*cb2*st2;
  //B = ct*sinf(2.0f*phi)*(ca2 - cb2);
  //C = (cb2*sf2 + ca2*cf2);

  //psi = 0.5f*atanf(B/(A-C));

  //f = st2*(ca2*cf2+cb2*sf2) + ct2;
  ///************/


  /**** NEW: siguiendo STARK********/
  float j,k,l;
  float t,u;
  t = 1./(ab*bc);
  u = 1./(bc);

  f = st2*(t*t*sf2 + u*u*cf2) + ct2;

  j = st2*t*t*u*u + t*t*cf2*ct2 + u*u*sf2*ct2;
  k = (u*u - t*t)*sf*cf*ct;
  l = t*t*sf2 + u*u*cf2;

  psi = 0.5f*atanf(2.0*k/(j-l));
  /************/

  //if(((A - C)*cosf(2.0f*psi) + B*sinf(2.0f*psi)) > 0.0f)
  //  psi = psi - PI_CUDA*0.5f;

  cp = cosf(psi);
  sp = sinf(psi);
  
  /**** NEW: siguiendo STARK********/
  /** X = A0^{-1} * A1^{-1} * A2^{-1} * X" **/
  xt = r[0];
  yt = r[1];
  zt = r[2];

  r[0] = cp*xt-sp*yt;
  r[1] = sp*xt+cp*yt;
  r[2] = zt;

  xt = r[0];
  yt = r[1];
  zt = r[2];

  r[0] = xt;
  r[1] = ct*yt-st*zt;
  r[2] = st*yt+ct*zt;

  xt = r[0];
  yt = r[1];
  zt = r[2];

  r[0] = cf*xt-sf*yt;
  r[1] = sf*xt+cf*yt;
  r[2] = zt;

  //r[0] = xt*(cf*cp - sf*ct*sp) - yt*(cf*sp + sf*ct*cp) + sf*st*zt;
  //r[1] = xt*(sf*cp + cf*ct*sp) + yt*(cf*ct*cp - sf*sp) - cf*st*zt;
  //r[2] = st*sp*xt + st*cp*yt + ct*zt;

  //r[0] = cp*xt-sp*yt;
  //r[1] = sp*xt+cp*yt;
  //r[2] = zt;

  //xt = r[0];
  //yt = r[1];
  //zt = r[2];

  //r[0] = cf*xt-ct*sf*yt+st*sf*zt;
  //r[1] = sf*xt+ct*cf*yt-st*cf*zt;
  //r[2] = st*yt+ct*zt;
  /*********************************/

  /*** Siguiendo Binney 1985 **/
  /** X = A0^{-1} * A1^{-1} * A2^{-1} * X" **/
  //xt = r[0];
  //yt = r[1];
  //zt = r[2];

  //r[0] = cp*xt-sp*yt;
  //r[1] = sp*xt+cp*yt;
  //r[2] = zt;

  //xt = r[0];
  //yt = r[1];
  //zt = r[2];

  //r[0] = xt;
  //r[1] = ct*yt-st*zt;
  //r[2] = st*yt+ct*zt;

  //xt = r[0];
  //yt = r[1];
  //zt = r[2];

  //r[0] = -sf*xt-cf*yt;
  //r[1] = cf*xt-sf*yt;
  //r[2] = zt;
  /****************************/


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
  tmp  = u(r,bc,ab,lgm,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));//sqrt(f);
#ifdef CG
	tmp *= RHOMEDIO/m;
  tmp *= mom1(lgm);
#endif
  tmp *= f_nu(nu)*nu*dlogNu_dlogM(lgm)*FMNORM;
  //tmp *= chingada(lgm);
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
  float ab2, bc2;
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
  ab2       = x[6];
  bc2       = x[7];
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

  //tmp  = u_esferico(y,lgm2,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]))/sqrt(f);
  tmp  = u(y,bc2,ab2,lgm2,sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));//sqrt(f);
  tmp *= (RHOMEDIO/m1);
  #ifdef CG
  tmp *= (RHOMEDIO/m2);
  tmp *= mom1(lgm2);
  #endif
  tmp *= f_nu(nu2)*nu2*dlogNu_dlogM(lgm2)*FMNORM;
  //tmp *= chingada(lgm2);
  /*FUNCION DE MASA FIT*/
  //tmp *= n(lgm2)*m2*m2/RHOMEDIO;
  tmp *= forma(bc2,ab2);
  //tmp *= NORMA_ORIENTACION; //factor para normalizar la func. de probabilidad de
                            //orientacion del halo vecino.
  tmp *= f_nu(nu1)*nu1*dlogNu_dlogM(lgm1)*FMNORM;
  //tmp *= chingada(lgm1);
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
  float theta;

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
      theta = 2.0*PI_CUDA*curand_uniform(&xx);
      p[0] = r*cosf(theta);
      p[1] = r*sinf(theta);

      //p[2] = curand_normal(&xx); /*Linea de la visual*/
      p[2] = RANGOZ*(2.0f*curand_uniform(&xx)-1.0f); /*Linea de la visual*/

      //tmp  = SQRT_TWOPI_CUDA*RANGOZ/exp(-p[2]*p[2]*0.5f);
      //p[2] *= RANGOZ;

      tmp = xilin(p[0],p[1],p[2]);
      tmp *= (2.0*RANGOZ);

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

__inline__ void proyectador_de_xilin(dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int i,j,k;
  FILE *pfout;
  char filename[80];
  float r,s;
  float dpaso,h_radio;
  dpaso = (PASOMAX - PASOMIN)/(float)NPASOS;

  sprintf(filename,"xilin_proyectada.dat");
  pfout = fopen(filename,"w");

  //En cada direccion hace NPASOS pasos//
  for(i = 0; i < NPASOS; i++){

    //Setea posicion en la direccion dada//
    h_radio = dpaso*(float)(i) + PASOMIN;
    h_radio = exp10(h_radio);

    //Lanza kernel//
    proyecta_xilin<<<dimGrid,dimBlock>>>(devStates,h_radio);
    cudaThreadSynchronize();
    //Termina kernel//

    //Copia sumatorias al host//
    HANDLE_ERROR(cudaMemcpyFromSymbol(h_int,d_int,(RNGS/THREADS_PER_BLOCK)*sizeof(float)));
    HANDLE_ERROR(cudaMemcpyFromSymbol(h_sig,d_sig,(RNGS/THREADS_PER_BLOCK)*sizeof(float)));

    //Termina de reducir en el host//
    r = 0.0f; s = 0.0f;
    for(k = 0; k < RNGS/THREADS_PER_BLOCK; k++){
      r += h_int[k];
      s += h_sig[k];
    }

    //Estima la integral y el sigma//
    r /= (float)((long)RNGS*(long)LAZOS);
    s /= (float)((long)RNGS*(long)LAZOS);
    s -= (r*r);
    s /= (float)((long)RNGS*(long)LAZOS);
    s  = sqrt(s);

    //Imprime en file de salida//
    fprintf(pfout,"%e %e %e\n",h_radio,r,s);
  }
  //Cierra archivo de salida//
  fclose(pfout);
}

