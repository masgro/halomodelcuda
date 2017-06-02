#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <omp.h>
#include <vector_types.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <time.h>

#include "constantes.h"
#include "HandleError.h"

//#include "HandleError.cu"

#define SIGMA   0.50000f
#define SIGMA2  0.25000f
#define SIGMA3  0.12500f

#ifndef ANGULO
#define ANGULO 45
#endif

//#define TRACERS_MASA_MIN 11.82f
#define TRACERS_MASA_MIN 10.82f
#define TRACERS_MASA_MAX 15.0f

/*Cantidad total de hilos (RNG) que se van a tirar*/
#define RNGS 65536

/*Cantidad de veces que se lanza el Kernel de integracion*/
#define LAZOS 500
#define LAZOSPLUS 10000

/*Cantidad de Threads por Block*/
#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 128
#endif

/*Tarjeta sobre la cual corre*/
#ifndef DEVICE
#define DEVICE 0
#endif

/*Intervalo espacial para la correlacion lineal*/
#define RANGO 10.0f
/*Grosor de la pizza*/
#define RANGOZ 5.0f

/*Numero de dimensiones de la integral*/
#ifdef DOSH
#define NDIM 15
#else
#define NDIM 7
#endif

/*Cantidad de pasos en cada direccion*/
#define NPASOS   50
#define PASOMIN -1.0
#define PASOMAX  2.0

#define NDIR 2

/*Vectores integral y sigma, host version*/
float h_int[RNGS/THREADS_PER_BLOCK];
float h_sig[RNGS/THREADS_PER_BLOCK];

/*Vectores integral y sigma, device version*/
__device__ float d_int[RNGS/THREADS_PER_BLOCK];
__device__ float d_sig[RNGS/THREADS_PER_BLOCK];
__device__ float d_integral[1],d_sigma[1];

/*Coeficientes de la forma y normalizacion, host version*/
float h_bc[4][3];
float h_ab[4][3];
float h_norm[4];

/*Coeficientes de la forma y normalizacion, device version*/
__device__ float d_bc[4][3];
__device__ float d_ab[4][3];
__device__ float d_norm[4];

/*Coeficientes del alineamiento, host version*/
float h_alig[10][30][5];
float h_alig_norm[10][30];
float h_alig_m[11];
float h_alig_rmin;
float h_alig_rmax;
float h_alig_dr;

/*Coeficientes del alineamiento, device version*/
__device__ float d_alig[10][30][5];
__device__ float d_alig_norm[10][30];
__device__ float d_alig_m[11];
__device__ float d_alig_rmin;
__device__ float d_alig_rmax;
__device__ float d_alig_dr;

/*Vectores de limites de las integrales, host version*/
float h_xmin[NDIM];
float h_xmax[NDIM];
float norma_funmasa;

/*Vectores de limites de las integrales, device version*/
__constant__ float d_xmin[NDIM];
__constant__ float d_xmax[NDIM];

#include "lecturas.cu"

#include "chrono.c"

/*Incluye archivo con todas las funciones necesarias*/
#include "common_functions.cu"
#include "funciones.cu"

/*Kernel: toma un punto aleatorio en el espacio N-Dimensional
          y evalua la funcion integrando (T1h, T2h) en dicho punto.
          La evaluacion se guarda en d_int y el cuadrado en d_sig.
          Al final guarda el estado del RNG en state*/
__global__ void integra(curandState *state, float r, int eje){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;

  /*Estado del RNG*/
  curandState seed = state[tid];

  const unsigned int it = threadIdx.x;
  int   i,j;
  float x[NDIM];
  float p[3];
  double value, sigma, tmp;

  __shared__ float xmin[NDIM], dx[NDIM];
  __shared__ float s_value[THREADS_PER_BLOCK];
  __shared__ float s_sigma[THREADS_PER_BLOCK];

  /*Inicializan variables*/
  if(it < NDIM){
    xmin[it] = d_xmin[it];
      dx[it] = d_xmax[it] - xmin[it];
  }

  /*Esperan hasta que todos terminen*/
  __syncthreads();

  /*Tira un numero random x de dimension NDIM*/
  value = 0.0; sigma = 0.0;
  for(j = 0; j < LAZOS; j++){
    do{
      //x[NDIM-1] = dx[NDIM-1] * curand_uniform(&sedd) + xmin[NDIM-1];

      switch(eje){
        case 2 :
          x[NDIM-2] = curand_uniform(&seed)*2.0f*PI_CUDA;
          p[0] = r*cosf(x[NDIM-2]);
          p[1] = r*sinf(x[NDIM-2]);
          break;

        case 1 :
          x[NDIM-2] = dx[NDIM-2]*curand_uniform(&seed)+xmin[NDIM-2];
          x[NDIM-2] = x[NDIM-2]*GRAD2RAD;
          p[0] = r*cosf(x[NDIM-2]);
          p[1] = r*sinf(x[NDIM-2]);
          break;

        case 0 :
          x[NDIM-2] = dx[NDIM-2]*curand_uniform(&seed)+xmin[NDIM-2];
          x[NDIM-2] = x[NDIM-2]*GRAD2RAD;
          p[1] = r*cosf(x[NDIM-2]);
          p[0] = r*sinf(x[NDIM-2]);
          break;
      }

      p[2] = curand_normal(&seed); /*Linea de la visual*/
      //p[2] = 2.0f*RANGOZ*curand_uniform(&seed) - RANGOZ; /*Linea de la visual*/

      //if(curand_uniform(&seed) < 0.5)
      //  p[0] = -p[0];
      //if(curand_uniform(&seed) < 0.5)
      //  p[1] = -p[1];

#ifdef DOSH
      for(i = 0; i <= 9; i++)
        x[i] = dx[i] * curand_uniform(&seed) + xmin[i];
      
      x[10] = curand_normal(&seed);
      x[11] = curand_normal(&seed);
      x[12] = curand_normal(&seed);

      tmp = x[10]*x[10] + x[11]*x[11] + x[12]*x[12];
      if(tmp < 1.E-32){
        tmp = 1.0E-16/sqrt(tmp);
        x[10] *= tmp;
        x[11] *= tmp;
        x[12] *= tmp;
        tmp = 1.0E-32;
      }

      /*sqrt(2·pi)^3 sigma^3 / exp(-tmp/2)*/
      tmp  = SQRT_TWOPI_CUBO_CUDA*SIGMA3*exp(tmp*0.5f);
      x[10] *= SIGMA;
      x[11] *= SIGMA;
      x[12] *= SIGMA;

      /*sqrt(2·pi) RANGOZ / exp(-z*z/2)*/
      //tmp  *= 2.0f*RANGOZ;
      tmp  *= SQRT_TWOPI_CUDA*RANGOZ*exp(p[2]*p[2]*0.5f);
      p[2] *= RANGOZ;

      tmp *= T2h(p,x);
#else 
      for(i = 0; i <= 4; i++)
        x[i] = dx[i] * curand_uniform(&seed) + xmin[i];

      /*sqrt(2·pi) RANGOZ / exp(-tmp/2)*/
      //tmp  = 2.0f*RANGOZ;
      tmp  = SQRT_TWOPI_CUDA*RANGOZ*exp(p[2]*p[2]*0.5f);
      p[2] *= RANGOZ;

      tmp *= T1h(p,x);
#endif

      if(isfinite(tmp))break;
    }while(1);

    value += tmp;
    sigma += tmp*tmp;
  }

  /*Guarda en la memoria compartida*/
  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();

  if(it == 0){
    value = 0.0; sigma = 0.0;
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
  state[tid] = seed;
}

__global__ void suma(void){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  int inext;

  inext = RNGS/THREADS_PER_BLOCK/2;
  while(inext >= 1){
    if(tid < inext){
      d_int[tid] += d_int[tid+inext];
      d_sig[tid] += d_sig[tid+inext];
    }
    inext = inext/2;
    __syncthreads();
  }

  if(tid == 0){
    d_integral[0] = d_int[0];
    d_sigma[0] = d_sig[0];
  }
}

/*Funcion que imprime las propiedades de la placa*/
__inline__ void printDevProp(cudaDeviceProp devProp){
  printf("#############################################\n");
  printf(" Running on device:         %s\n",  devProp.name);
  printf(" Total global memory:       %zu\n", devProp.totalGlobalMem);
  printf(" Total shared memory/block: %zu\n", devProp.sharedMemPerBlock);
  printf(" Total registers/block:     %d\n",  devProp.regsPerBlock);
  printf(" Maximum threads/block:     %d\n",  devProp.maxThreadsPerBlock);
  printf(" Total constant memory:     %zu\n", devProp.totalConstMem);
  printf(" Kernel execution timeout:  %s\n", (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
  printf(" WarpSize:                  %d\n",devProp.warpSize);
  printf(" Compute Capability:        %d%d\n",devProp.major,devProp.minor);
  printf("#############################################\n");
}

int setseed(void){
  int seed;
  FILE *pf;
  pf = fopen("/dev/urandom","r");
  fread(&seed,sizeof(int),1,pf);
  fclose(pf);
  return seed;
}

int main(int argc, char **argv){
  double time,elapsed;
  double time1;
  FILE  *pfout;
  char  filename[200],term[200];
  int   i,j,l;
  int   run;
  float r,s;
  float volumen;
  float h_radio;
  //float *d_radio;
  //int   *d_eje;
  curandState *devStates;

  elapsed = 0.0f;
  chrono(START,&time);

  run = (argc > 1)? atoi(argv[1]) : 0;

#ifdef DOSH
  sprintf(term,"2h");
#else
  sprintf(term,"1h");
#endif

  /*Setea el device a utilizar*/
  cudaSetDevice(DEVICE);

  /*Lee e imprime las propiedades del device*/
  cudaDeviceProp devProp;
  cudaGetDeviceProperties(&devProp, DEVICE);
  printDevProp(devProp);

  //HANDLE_ERROR(cudaMalloc((void **)&d_radio,sizeof(float)));
  //HANDLE_ERROR(cudaMalloc((void **)&d_eje,sizeof(int)));
  
  /*Chequea Cantidad de Threads y de Blocks*/
  assert(THREADS_PER_BLOCK <= 1024);
  assert(RNGS%THREADS_PER_BLOCK == 0); // should be divisible by blocks

  /*Setea Cantidad de Threads y de Blocks*/
  dim3 dimBlock(THREADS_PER_BLOCK,1,1);
  dim3 dimGrid(RNGS/THREADS_PER_BLOCK,1,1);

  fprintf(stdout,"Corriendo %d Blocks con %d threads cada uno\n",
                        RNGS/THREADS_PER_BLOCK,THREADS_PER_BLOCK);

  /*Allocatea memoria para el RNG*/
  HANDLE_ERROR(cudaMalloc((void **)&devStates,RNGS*sizeof(curandState)));

  /*Setea la semilla*/
  int seed = setseed();

  /*Setea las semillas de los RNG en el device*/
  setup_kernel<<<dimGrid,dimBlock>>>(devStates,seed);
  cudaThreadSynchronize();

  /*lee los coeficientes de los ajustes*/
  read_coefficients();

  float Numin   = Nu_M(CENTROS_MASA_MIN);
  float Numax   = Nu_M(CENTROS_MASA_MAX);
  float ncmedio = nc_medio(Numin,Numax,dimGrid,dimBlock,devStates);
  /*FUNCION DE MASA FIT*/
  //float ncmedio = nc_medio(CENTROS_MASA_MIN,CENTROS_MASA_MAX,dimGrid,dimBlock,devStates);

#ifdef CG 
  Numin   = Nu_M(TRACERS_MASA_MIN);
  Numax   = Nu_M(TRACERS_MASA_MAX);
  float ngmedio = ng_medio(Numin,Numax,dimGrid,dimBlock,devStates);
  /*FUNCION DE MASA FIT*/
  //float ngmedio = ng_medio(TRACERS_MASA_MIN,TRACERS_MASA_MAX,dimGrid,dimBlock,devStates);
#endif

  /*Chequea si la integral del perfil hasta el radio virial da 1*/
  float norma_perfil = normalizacion_perfil(dimGrid,dimBlock,devStates);
  
  prueba_ng_medio();

  /*Verifica que la integral de bias(nu)*f(nu) de 1*/
  //prueba_bias_f(dimGrid,dimBlock,devStates);

  /*Calcula la normalizacion de la funcion de masa*/
  //normalizacion_func_masa(dimGrid,dimBlock,devStates);

  /*Setea los limites de integracion*/
  /*Halo Centro*/
  h_xmin[0] = (float)CENTROS_MASA_MIN; /*Masa minima*/
  h_xmax[0] = (float)CENTROS_MASA_MAX; /*Masa maxima*/
  /*Forma*/
  h_xmin[1] = 0.00f; /* ab minimo */
  h_xmax[1] = 1.00f; /* ab maximo */
  h_xmin[2] = 0.00f; /* bc minimo */
  h_xmax[2] = 1.00f; /* bc maximo */

  //h_xmin[1] = ABMEDIO; /* ab minimo */
  //h_xmax[1] = ABMEDIO; /* ab maximo */
  //h_xmin[2] = BCMEDIO; /* bc minimo */
  //h_xmax[2] = BCMEDIO; /* bc maximo */

  h_xmin[3] =  0.0f;     /*Orientacion del Halo Centro*/ 
  h_xmax[3] =  1.0f;      
  h_xmin[4] =  0.0f;      
  h_xmax[4] =  0.5f*M_PI; 

#ifdef DOSH
  /*Halo Vecino*/
  h_xmin[5] = (float)TRACERS_MASA_MIN; /*Masa minima*/
  h_xmax[5] = (float)TRACERS_MASA_MAX; /*Masa maxima*/
  h_xmin[6] = 0.0f; /* ab minimo */
  h_xmax[6] = 0.0f; /* ab maximo */
  h_xmin[7] = 0.0f; /* bc minimo */
  h_xmax[7] = 0.0f; /* bc maximo */

  h_xmin[8] =  0.0f;
  h_xmax[8] =  0.0f;
  h_xmin[9] =  0.0f;
  h_xmax[9] =  0.0f;

  h_xmin[10] =  0.0f;
  h_xmax[10] =  0.0f;
  h_xmin[11] =  0.0f;
  h_xmax[11] =  0.0f;
  h_xmin[12] =  0.0f;
  h_xmax[12] =  0.0f;
#endif

  /** Orientacion **/
  h_xmin[NDIM-2] = 0.0f;
  h_xmax[NDIM-2] = (float)ANGULO;
  h_xmin[NDIM-1] = 0.0f;
  h_xmax[NDIM-1] = 0.0f;

  clock_t cuenta;
  //cuenta = clock();
  //float norma_merchan = integra_merchan(dimGrid,dimBlock,devStates);
  //time1 = ((double)(clock()-cuenta))/((double)CLOCKS_PER_SEC);
  //printf("  NormaForma: %E time %.15E\n",norma_merchan,time1);

  cuenta = clock();
  double i1,i2,a,b;
  a = ABMEDIO; b = 0.1;
  i1  = (2.0*(a-1.0)*b*exp(-a*a*0.5/b/b));
  i1 -= (sqrt(2.0*M_PI)*((a-1.0)*a+b*b)*erf((a-1.0)/sqrt(2.0)/b));
  i1 += (sqrt(2.0*M_PI)*((a-1.0)*a+b*b)*erf(a/sqrt(2.0)/b));
  i1 -= (2.*a*b*exp(-(a-1.0)*(a-1.0)*0.5/b/b));
  i1 *= (-0.5*b);
  i1 /= (sqrt(2.0*M_PI)*b);

  a = BCMEDIO; b = 0.1;
  i2  = (2.0*(a-1.0)*b*exp(-a*a*0.5/b/b));
  i2 -= (sqrt(2.0*M_PI)*((a-1.0)*a+b*b)*erf((a-1.0)/sqrt(2.0)/b));
  i2 += (sqrt(2.0*M_PI)*((a-1.0)*a+b*b)*erf(a/sqrt(2.0)/b));
  i2 -= (2.*a*b*exp(-(a-1.0)*(a-1.0)*0.5/b/b));
  i2 *= (-0.5*b);
  i2 /= (sqrt(2.0*M_PI)*b);

  float norma_merchan = i1*i2;
  time1 = ((double)(clock()-cuenta))/((double)CLOCKS_PER_SEC);
  printf("  NormaForma: %E time %.15E\n",norma_merchan,time1);
  


#ifdef DOSH
  float norma_align = integra_align(dimGrid,dimBlock,devStates);
#endif

  /*Copia los limites de integracion al device*/
  HANDLE_ERROR(cudaMemcpyToSymbol(d_xmin,h_xmin,NDIM*sizeof(float)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_xmax,h_xmax,NDIM*sizeof(float)));

  /*Calcula el hipervolumen de integracion*/
  volumen = 1.0f;
#ifdef DOSH
  for(i = 0; i <= 4; i++) volumen *= (h_xmax[i] - h_xmin[i]);
  //for(i = 3; i <= 4; i++) volumen *= (h_xmax[i] - h_xmin[i]);
  for(i = 5; i <= 5; i++) volumen *= (h_xmax[i] - h_xmin[i]);
#else
  for(i = 0; i <= 4; i++) volumen *= (h_xmax[i] - h_xmin[i]);
  //for(i = 3; i <= 4; i++) volumen *= (h_xmax[i] - h_xmin[i]);
#endif

  /*Calcula la memoria total en el device*/
  size_t memfree, memtot;
  HANDLE_ERROR(cudaMemGetInfo(&memfree,&memtot));
  printf("Memoria CUDA Total: %8.3lf Mb Used: %8.3lf Mb free: %8.3lf Mb \n",
         (float)memtot/1024.0f/1024.0f,(float)(memtot-memfree)/1024.0f/1024.0f,
         (float)memfree/1024.0f/1024.0f);

  /*Imprime alguna informacion*/
  printf("--------------------------\n");
  printf("  Volumen:  %E\n",volumen);
  printf("  RHOMEDIO: %E\n",RHOMEDIO);
  printf("  NCMEDIO:  %E\n",ncmedio);
#ifdef CG
  printf("  NGMEDIO: %E\n",ngmedio);
#endif
  printf("  NormaForma: %E\n",norma_merchan);
#ifdef DOSH
  printf("  NormaAlign: %E\n",norma_align);
#endif
  printf("  RNGs:  %d\n",RNGS);
  printf("  LAZOS: %d\n",LAZOS);
  printf("--------------------------\n");

  printf("Computando termino de %s....\n",term);

  float dpaso;
  dpaso = (PASOMAX - PASOMIN)/(float)NPASOS;

  /*Recorre las 3 direcciones j=0(parallel),1(perpendicular),2(iso)*/
  for(j = 0; j < 3; j++){
    /*Abre archivo de salida*/
    sprintf(filename,"funcorr_%1d_%s.%02d",j,term,run);
    pfout = fopen(filename,"w");

    /*En cada direccion hace NPASOS pasos*/
    for(i = 0; i < NPASOS; i++){
      /*Setea posicion en la direccion dada*/
      h_radio = dpaso*(float)(i) + PASOMIN;
      h_radio = powf(10.0f,h_radio);

      /*Copia posicion al device*/
      //HANDLE_ERROR(cudaMemcpy(d_radio,&h_radio,sizeof(float),cudaMemcpyHostToDevice));
      //HANDLE_ERROR(cudaMemcpy(d_eje,&j,sizeof(int),cudaMemcpyHostToDevice));

      /*Lanza kernel*/
      integra<<<dimGrid,dimBlock>>>(devStates,h_radio,j);
      cudaThreadSynchronize();

      //suma<<<dimGrid,dimBlock>>>();
      //cudaThreadSynchronize();
      /*Termina kernel*/

      CHECK_KERNEL_SUCCESS();

      ////*Copia sumatorias al host y termina de reducir en el host*/
      HANDLE_ERROR(cudaMemcpyFromSymbol(h_int,d_int,(RNGS/THREADS_PER_BLOCK)*sizeof(float)));
      HANDLE_ERROR(cudaMemcpyFromSymbol(h_sig,d_sig,(RNGS/THREADS_PER_BLOCK)*sizeof(float)));

      /*Termina de reducir en el host*/
      r = 0.0; s = 0.0;
      for(l = 0; l < RNGS/THREADS_PER_BLOCK; l++){
        if(isnan(h_int[l]))printf("%f %d %d\n",h_int[l],i,j);
        r += h_int[l];
        s += h_sig[l];
      }

      //HANDLE_ERROR(cudaMemcpyFromSymbol(&r,d_integral,sizeof(float)));
      //HANDLE_ERROR(cudaMemcpyFromSymbol(&s,d_sigma,sizeof(float)));

      /*Estima la integral y el sigma*/
      r /= (float)((long)RNGS*(long)LAZOS);
      s /= (float)((long)RNGS*(long)LAZOS);
      s -= (r*r);
      s /= (float)((long)RNGS*(long)LAZOS);
      s  = sqrt(s);

      s *= volumen;
      r *= volumen;

      r /= ncmedio;
      s /= ncmedio;

      r /= norma_merchan;
      s /= norma_merchan;

#ifndef DOSH
      r /= norma_perfil;
      s /= norma_perfil;
#endif

#ifdef DOSH
      r /= norma_align;
      s /= norma_align;
      //r /= norma_merchan;
      //s /= norma_merchan;
#endif

#ifdef CG
			r /= ngmedio;
			s /= ngmedio;
#endif

      /*Imprime en file de salida*/
      fprintf(pfout,"%e %e %e\n",h_radio,r,s);
    }
    /*Cierra archivo de salida*/
    fclose(pfout);
  }

  /*Libera memoria allocateada en el device*/
  HANDLE_ERROR(cudaFree(devStates));

  /*Computa el tiempo total utilizado en el device*/
  chrono(STOP,&time);
  elapsed += time;
  printf("Tiempo: %lf [seg] \n", elapsed);

  /*Fin del programa*/
  return(EXIT_SUCCESS);
}
