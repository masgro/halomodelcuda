#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <omp.h>
#include "cutil.h"
#include "constantes.h"
#include <cuda.h>
#include <curand_kernel.h>

#define TRACERS_MASA_MIN 9.0f
#define TRACERS_MASA_MAX 17.0f

#ifndef CENTROS_MASA_MIN
#define CENTROS_MASA_MIN 13.0f
#endif
#ifndef CENTROS_MASA_MAX
#define CENTROS_MASA_MAX 13.5f
#endif

#define SIGMA    0.5f
#define SIGMA3 0.125f

#define BMAX 0.5f     
#define BMAX_YZ 0.3f  

#ifndef ANGULO
#define ANGULO 1
#endif

/*Cantidad total de hilos (RNG) que se van a tirar*/
#define RNGS 65536

/*Cantidad de veces que se lanza el Kernel de integracion*/
#define LAZOS 100
#define LAZOSPLUS 100000

/*Cantidad de Threads por Block*/
#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 512
#endif

/*Tarjeta sobre la cual corre*/
#ifndef DEVICE
#define DEVICE 1   /*en quadro, DEVICE 0: gtx, DEVICE 1: quadro*/
#endif

/*Intervalo espacial para la correlacion lineal*/
#ifdef DOSH
#define RANGO 10.0f
#endif

/*Numero de dimensiones de la integral*/
#ifdef DOSH
#define NDIM 12
#else
#define NDIM 3
#endif

/*Cantidad de pasos en cada direccion*/
#define NPASOS   100
#define PASOMIN -1.0
#define PASOMAX  2.0

/*Vectores integral y sigma, host version*/
float h_int[RNGS/THREADS_PER_BLOCK];
float h_sig[RNGS/THREADS_PER_BLOCK];

/*Vectores integral y sigma, device version*/
__device__ float d_int[RNGS/THREADS_PER_BLOCK];
__device__ float d_sig[RNGS/THREADS_PER_BLOCK];

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
__device__ float d_xmin[NDIM];
__device__ float d_xmax[NDIM];


#include "lecturas.cu"

#include "chrono.c"

/*Incluye archivo con todas las funciones necesarias*/
#include "funciones.cu"

/*Kernel: toma un punto aleatorio en el espacio N-Dimensional
          y evalua la funcion integrando (T1h, T2h) en dicho punto.
          La evaluacion se guarda en d_int y el cuadrado en d_sig.
          Al final guarda el estado del RNG en state*/
__global__ void integra(curandState *state, float *r)
{
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;

  /*Estado del RNG*/
  curandState xx = state[tid];

  const unsigned int it = threadIdx.x;
  int   i,j;
  float x[NDIM];

  float value, sigma, tmp;

  __shared__ float xmin[NDIM], dx[NDIM];
  __shared__ float s_value[THREADS_PER_BLOCK];
  __shared__ float s_sigma[THREADS_PER_BLOCK];
  __shared__ float p[3];

  /*Inicializan variables*/
  if(it < 3) p[it] = r[it]; /*Setea las posicion*/

  if(it < NDIM)
  {
    xmin[it] = d_xmin[it];
      dx[it] = d_xmax[it] - xmin[it];
  }

  /*Esperan hasta que todos terminen*/
  __syncthreads();

  /*Tira un numero random x de dimension NDIM*/
  value = 0.0; sigma = 0.0;
  for(j = 0; j < LAZOS; j++)
  {
#ifdef DOSH
    for(i = 0; i < 6; i++)
      x[i] = dx[i] * curand_uniform(&xx) + xmin[i];
    
    x[6] = curand_normal(&xx);
    x[7] = curand_normal(&xx);
    x[8] = curand_normal(&xx);

    tmp = x[6]*x[6] + x[7]*x[7] + x[8]*x[8];
    if(tmp < 1.E-4)
    {
      tmp = 1.0E-2/sqrt(tmp);
      x[6] *= tmp;
      x[7] *= tmp;
      x[8] *= tmp;
      tmp = 1.0E-4;
    }

    for(i = 9; i < NDIM; i++)
      x[i] = dx[i] * curand_uniform(&xx) + xmin[i];

    tmp  = SQRT_TWOPI_CUBO_CUDA*SIGMA3/exp(-tmp*0.5f);

    x[6] *= SIGMA;
    x[7] *= SIGMA;
    x[8] *= SIGMA;

    x[6] += p[0];
    x[7] += p[1];
    x[8] += p[2];

    tmp *= T2h(p,x);
#else 
    for(i = 0; i < NDIM; i++)
      x[i] = dx[i] * curand_uniform(&xx) + xmin[i];

    tmp = T1h(p,x);
#endif

    value += tmp;
    sigma += tmp*tmp;
  }

  /*Guarda en la memoria compartida*/
  s_value[it] = value;
  s_sigma[it] = sigma;

  __syncthreads();

  if(it == 0)
  {
    value = 0.0; sigma = 0.0;
    for(i = 0; i < THREADS_PER_BLOCK; i++)
    {
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

/*Funcion que imprime las propiedades de la placa*/
__inline__ void printDevProp(cudaDeviceProp devProp)
{
  printf("Running on device:             %s\n",  devProp.name);
  printf("Total global memory:           %zu\n",  devProp.totalGlobalMem);
  printf("Total shared memory per block: %zu\n",  devProp.sharedMemPerBlock);
  printf("Total registers per block:     %d\n",  devProp.regsPerBlock);
  printf("Maximum threads per block:     %d\n",  devProp.maxThreadsPerBlock);
  printf("Total constant memory:         %zu\n",  devProp.totalConstMem);
  printf("Kernel execution timeout:      %s\n", (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
  return;
}

int main(int argc, char **argv)
{
  FILE *pfout;
  char filename[200];
  int i,j,l;
  float r, s;
  float volumen;
  float *h_pos;
  float *d_pos;
  curandState *devStates;
  


  int run;
  if(argc > 1)
    run = atoi(argv[1]);
  else
    run = 0;

  /*Setea el device a utilizar*/
  cudaSetDevice(DEVICE);

  /*Lee e imprime las propiedades del device*/
  cudaDeviceProp devProp;
  cudaGetDeviceProperties(&devProp, DEVICE);
  printDevProp(devProp);




#ifdef GRAMCHARLIER
  /*Lee los parametros de los ajustes*/
  lee_parametros();
#endif

  /*Chequea Cantidad de Threads y de Blocks*/
  assert(THREADS_PER_BLOCK <= 1024);
  assert(RNGS%THREADS_PER_BLOCK == 0); // should be divisible by blocks

  /*Setea Cantidad de Threads y de Blocks*/
  dim3 dimBlock(THREADS_PER_BLOCK,1,1);
  dim3 dimGrid(RNGS/THREADS_PER_BLOCK,1,1);

  fprintf(stdout,"Corriendo %d Blocks con %d threads cada uno\n",RNGS/THREADS_PER_BLOCK,THREADS_PER_BLOCK);

  /*Allocatea memoria para el RNG*/
  CUDA_SAFE_CALL(cudaMalloc((void **)&devStates,RNGS*sizeof(curandState)));

  /*Setea las semillas de los RNG*/
  int entero;
  FILE *pfrandom;
  pfrandom = fopen("/dev/urandom","r");
  fread(&entero,sizeof(int),1,pfrandom);

  setup_kernel<<<dimGrid,dimBlock>>>(devStates,entero);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("run kernel failed");

  /*lee los coeficientes de los ajustes*/
  read_coefficients();

  float ncmedio;
  float Numin, Numax;
  Numin   = Nu_M(CENTROS_MASA_MIN);
  Numax   = Nu_M(CENTROS_MASA_MAX);
  ncmedio = nc_medio(Numin,Numax,dimGrid,dimBlock,devStates);

#ifdef CG 
  Numin   = Nu_M(TRACERS_MASA_MIN);
  Numax   = Nu_M(TRACERS_MASA_MAX);
  float ngmedio = ng_medio(Numin,Numax,dimGrid,dimBlock,d_states);
#endif

  /*Calcula la normalizacion de la funcion de masa*/
  //normalizacion_func_masa(dimGrid,dimBlock,devStates);

  /*Setea los limites de integracion*/
  /*Halo Centro*/
  h_xmin[0] = (float)CENTROS_MASA_MIN; /*Masa minima*/
  h_xmax[0] = (float)CENTROS_MASA_MAX; /*Masa maxima*/
  h_xmin[1] = 0.00f; /* ab minimo */
  h_xmax[1] = 1.00f; /* ab maximo */
  h_xmin[2] = 0.00f; /* bc minimo */
  h_xmax[2] = 1.00f; /* bc maximo */

#ifdef DELTAS
  h_xmin[1] = ABMEDIO - ANCHO; /* ab minimo */
  h_xmax[1] = ABMEDIO + ANCHO; /* ab maximo */
  h_xmin[2] = BCMEDIO - ANCHO; /* bc minimo */
  h_xmax[2] = BCMEDIO + ANCHO; /* bc maximo */
#endif

#ifdef DOSH
  /*Halo Vecino*/
  h_xmin[3] = TRACERS_MASA_MIN; /*Masa minima*/
  h_xmax[3] = TRACERS_MASA_MAX; /*Masa maxima*/
  h_xmin[4] = 0.00f; /* ab minimo */
  h_xmax[4] = 1.00f; /* ab maximo */
  h_xmin[5] = 0.00f; /* bc minimo */
  h_xmax[5] = 1.00f; /* bc maximo */

#ifdef DELTAS
  h_xmin[4] = ABMEDIO - ANCHO; /* ab minimo */
  h_xmax[4] = ABMEDIO + ANCHO; /* ab maximo */
  h_xmin[5] = BCMEDIO - ANCHO; /* bc minimo */
  h_xmax[5] = BCMEDIO + ANCHO; /* bc maximo */
#endif

  /*Volumen*/
  h_xmin[6] = -RANGO; 
  h_xmax[6] = +RANGO;
  h_xmin[7] = -RANGO; 
  h_xmax[7] = +RANGO;
  h_xmin[8] = -RANGO; 
  h_xmax[8] = +RANGO;

  /*Orientacion del halo vecino*/
  h_xmin[9]  =  0.0f;
  h_xmax[9]  =  M_PI;
  h_xmin[10] =  0.0f;
  h_xmax[10] =  1.0f;
  h_xmin[11] =  0.0f;
  h_xmax[11] =  M_PI;
#endif

#ifdef MERCHAN
  float norma_merchan = integra_merchan(dimGrid,dimBlock,devStates);
#endif
  //float norma_align = integra_align(dimGrid,dimBlock,devStates);
  //float norma_align = 7.309210e-01; //0.91 0.35
  //float norma_align = 1.023280e+00; //0.91 0.1
  //float norma_align = 1.226160e+00; //0.91 1.0
  //float norma_align = 1.174561e-01; //0.30 0.35
  //float norma_align = 1.265457e+00; //3.00 0.35
  float norma_align = 3.062452e-01; //3.00 0.35
  /*Copia los limites de integracion al device*/
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_xmin,h_xmin,NDIM*sizeof(float)));
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_xmax,h_xmax,NDIM*sizeof(float)));

  /*Allocatea memoria para el vector posicion en host y device*/
  h_pos = (float *) calloc(3,sizeof(float));
  CUDA_SAFE_CALL(cudaMalloc((float **)&d_pos,3*sizeof(float)));

  /*Calcula el hipervolumen de integracion*/
  volumen = 1.0;
#ifdef DOSH
  for(i = 0; i < 6 ; i++) volumen *= (h_xmax[i] - h_xmin[i]);
  for(i = 9; i < NDIM; i++) volumen *= (h_xmax[i] - h_xmin[i]);
#ifdef MERCHAN
  volumen /= (norma_merchan*norma_merchan);
#endif
#else
  for(i = 0; i < NDIM; i++) volumen *= (h_xmax[i] - h_xmin[i]);
#ifdef MERCHAN
  volumen /= norma_merchan;
#endif
#endif

  /*Calcula la memoria total en el device*/
  size_t memfree, memtot;
  CUDA_SAFE_CALL(cudaMemGetInfo(&memfree,&memtot));
  printf("Memoria CUDA Total: %8.3lf Mb Used: %8.3lf Mb free: %8.3lf Mb \n",
         (float)memtot/1024.0/1024.0,(float)(memtot-memfree)/1024.0/1024.0,
         (float)memfree/1024.0/1024.0);
  
  /*Imprime alguna informacion*/
  printf("--------------------------\n");
  printf("  Volumen:  %E\n",volumen);
  printf("  RHOMEDIO: %E\n",RHOMEDIO);
  printf("  NCMEDIO:  %E\n",ncmedio);
#ifdef CG
  printf("  NGMEDIO: %E\n",ngmedio);
#endif
  printf("  RNGs:  %d\n",RNGS);
  printf("  LAZOS: %d\n",LAZOS);
  printf("--------------------------\n");

#ifdef DOSH
  printf("Computando termino de 2-Halos....\n");
#else
  printf("Computando termino de 1-Halo.....\n");
#endif

  float dpaso;
  dpaso = (PASOMAX - PASOMIN)/(float)NPASOS;

  /*Recorre las 3 direcciones j=0(x),1(y),2(z)*/
  for(j = 0; j < 3; j++)
  {
    /*Abre archivo de salida*/
#ifdef DOSH
    //sprintf(filename,"funcorr_%1d_2h_%4.1f-%4.1f_%2s_%02d.dat",j,CENTROS_MASA_MIN,CENTROS_MASA_MAX,argv[2],run);
    sprintf(filename,"funcorr_%1d_2h.%02d",j,run);
#else
    //sprintf(filename,"funcorr_%1d_1h_%4.1f-%4.1f_%2s_%02d.dat",j,CENTROS_MASA_MIN,CENTROS_MASA_MAX,argv[2],run);
    sprintf(filename,"funcorr_%1d_1h.%02d",j,run);
#endif
    pfout = fopen(filename,"w");

    h_pos[0] = 0.0f;
    h_pos[1] = 0.0f;
    h_pos[2] = 0.0f;

#ifdef DOSH
    h_xmin[6] = h_pos[0] - RANGO; h_xmax[6] = h_pos[0] + RANGO;
    h_xmin[7] = h_pos[1] - RANGO; h_xmax[7] = h_pos[1] + RANGO;
    h_xmin[8] = h_pos[2] - RANGO; h_xmax[8] = h_pos[2] + RANGO;
#endif

    /*En cada direccion hace NPASOS pasos*/
    for(i = 0; i < NPASOS; i++)
    {
      /*Setea posicion en la direccion dada*/
      h_pos[j] = dpaso*(float)(i) + PASOMIN;
      h_pos[j] = exp10(h_pos[j]);

      /*Copia posicion al device*/
      CUDA_SAFE_CALL(cudaMemcpy(d_pos,h_pos,3*sizeof(float),cudaMemcpyHostToDevice));

      /*Lanza kernel*/
      integra<<<dimGrid,dimBlock>>>(devStates, d_pos);
      cudaThreadSynchronize();
      CUT_CHECK_ERROR("run kernel failed");
      /*Termina kernel*/

      /*Copia sumatorias al host*/
      CUDA_SAFE_CALL(cudaMemcpyFromSymbol(h_int,d_int,(RNGS/THREADS_PER_BLOCK)*sizeof(float)));
      CUDA_SAFE_CALL(cudaMemcpyFromSymbol(h_sig,d_sig,(RNGS/THREADS_PER_BLOCK)*sizeof(float)));

      /*Termina de reducir en el host*/
      r = 0.0; s = 0.0;
      for(l = 0; l < RNGS/THREADS_PER_BLOCK; l++)
      {
        r += h_int[l];
        s += h_sig[l];
      }

      /*Estima la integral y el sigma*/
#ifdef DOSH
      r *= RHOMEDIO;
      s *= RHOMEDIO*RHOMEDIO;
#endif
      r /= (float)((long)RNGS*(long)LAZOS);
      s /= (float)((long)RNGS*(long)LAZOS);
      s -= (r*r);
      s /= (float)((long)RNGS*(long)LAZOS);
      s  = sqrt(s);
      s *= (volumen/ncmedio);
      r *= (volumen/ncmedio);
#ifdef DOSH
      s /= norma_align;
      r /= norma_align;
#endif

#ifdef CG
			r /= ngmedio;
			s /= ngmedio;
#endif
      /*Imprime en file de salida*/
      fprintf(pfout,"%e %e %e\n",h_pos[j],r,s);
    }
    /*Cierra archivo de salida*/
    fclose(pfout);
  }

  /*Libera memoria allocateada en el device*/
  CUDA_SAFE_CALL(cudaFree(d_pos));
  CUDA_SAFE_CALL(cudaFree(devStates));

  /*Fin del programa*/
  return(EXIT_SUCCESS);
}
