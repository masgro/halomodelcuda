#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <omp.h>
#include "constantes.h"
#include "colores.h"
#include <vector_types.h>
#include <cuda.h>
#include <curand_kernel.h>

#include "HandleError.cu"

#define TRACERS_MASA_MIN 9.0f
#define TRACERS_MASA_MAX 18.0f

/*Cantidad total de hilos (RNG) que se van a tirar*/
#define RNGS 65536

/*Cantidad de veces que se lanza el Kernel de integracion*/
#define LAZOS 1000
#define LAZOSPLUS 100000

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

/*Numero de dimensiones de la integral*/
#ifdef DOSH
#define NDIM 14
#else
#define NDIM 5
#endif

/*Cantidad de pasos en cada direccion*/
#define NPASOS   50
#define PASOMIN -1.0
#define PASOMAX  2.0

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
__device__ float d_xmin[NDIM];
__device__ float d_xmax[NDIM];

#include "lecturas.cu"

#include "chrono.c"

/*Incluye archivo con todas las funciones necesarias*/
#include "funciones.cu"

/*Funcion que imprime las propiedades de la placa*/
void printDevProp(cudaDeviceProp devProp){
  printf("Running on device:             %s\n",devProp.name);
  printf("Total global memory:           %zu\n",devProp.totalGlobalMem);
  printf("Total shared memory per block: %zu\n",devProp.sharedMemPerBlock);
  printf("Total registers per block:     %d\n",devProp.regsPerBlock);
  printf("Maximum threads per block:     %d\n",devProp.maxThreadsPerBlock);
  printf("Total constant memory:         %zu\n",devProp.totalConstMem);
  printf("Kernel execution timeout:      %s\n",(devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
  printf("WarpSize:                      %d\n",devProp.warpSize);
  printf("Compute Capability:            %d%d\n",devProp.major,devProp.minor);
}

int main(int argc, char **argv){
  float time, elapsed;
  FILE  *pfout;
  char  filename[200];
  int   i,j,l;
  float r, s;
  float volumen;
  float h_radio;
  curandState *devStates;
	size_t ierr;

  elapsed = 0.0f;
  chrono(START,&time);

  int run = 0;
  if(argc > 1) run = atoi(argv[1]);

  /*Setea el device a utilizar*/
  cudaSetDevice(DEVICE);

  /*Lee e imprime las propiedades del device*/
  cudaDeviceProp devProp;
  cudaGetDeviceProperties(&devProp, DEVICE);
  printDevProp(devProp);

  /*Chequea Cantidad de Threads y de Blocks*/
  assert(THREADS_PER_BLOCK <= 1024);
  assert(RNGS%THREADS_PER_BLOCK == 0); // should be divisible by blocks

  /*Setea Cantidad de Threads y de Blocks*/
  dim3 dimBlock(THREADS_PER_BLOCK,1,1);
  dim3 dimGrid(RNGS/THREADS_PER_BLOCK,1,1);

  fprintf(stdout,"Corriendo %d Blocks con %d threads cada uno\n",RNGS/THREADS_PER_BLOCK,THREADS_PER_BLOCK);

  /*Allocatea memoria para el RNG*/
  HANDLE_ERROR(cudaMalloc((void **)&devStates,RNGS*sizeof(curandState)));

  /*Setea las semillas de los RNG*/
  int entero;
  FILE *pfrandom;
  pfrandom = fopen("/dev/urandom","r");
  ierr = fread(&entero,sizeof(int),1,pfrandom);
	assert(ierr == 1);

  setup_kernel<<<dimGrid,dimBlock>>>(devStates,entero);
  cudaThreadSynchronize();

  /*lee los coeficientes de los ajustes*/
  read_coefficients();

  float ncmedio;
  float Numin, Numax;
  Numin   = Nu_M(CENTROS_MASA_MIN);
  Numax   = Nu_M(CENTROS_MASA_MAX);
  ncmedio = nc_medio(Numin,Numax,dimGrid,dimBlock,devStates);

  float ngmedio;
  Numin   = Nu_M(TRACERS_MASA_MIN);
  Numax   = Nu_M(TRACERS_MASA_MAX);
  ngmedio = ng_medio(Numin,Numax,dimGrid,dimBlock,devStates);

  printf("--------------------------\n");
  printf("  ngmedio: %E\n",ngmedio);
  printf("  ncmedio: %E\n",ncmedio);
  printf("--------------------------\n");

  prueba_ng_medio();

  /*Libera memoria allocateada en el device*/
  cudaFree(devStates);

  /*Computa el tiempo total utilizado en el device*/
  chrono(STOP,&time);
  elapsed += time;
  printf("Tiempo: %lf [seg] \n", elapsed);

  /*Fin del programa*/
  return(EXIT_SUCCESS);
}
