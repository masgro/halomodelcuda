#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <omp.h>
#include <vector_types.h>
#include <cuda.h>
#include <curand_kernel.h>

#include "constantes.h"
#include "colores.h"

#include "HandleError.cu"

#define SIGMA   2.000f
#define SIGMA2  4.000f
#define SIGMA3  8.000f

#ifndef ANGULO
#define ANGULO 45
#endif

#define TRACERS_MASA_MIN 11.81f
#define TRACERS_MASA_MAX 16.0f

/*Cantidad total de hilos (RNG) que se van a tirar*/
#define RNGS 65536

/*Cantidad de veces que se lanza el Kernel de integracion*/
#define LAZOS 500
#define LAZOSPLUS 1000

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
#define RANGOZ 20.0f

/*Numero de dimensiones de la integral*/
#ifdef DOSH
#define NDIM 15
#else
#define NDIM 7
#endif

/*Cantidad de pasos en cada direccion*/
#define NPASOS   50
#define PASOMIN -1.0
#define PASOMAX  1.0

#define NDIR 3

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

__device__ float d_abmedio;
__device__ float d_bcmedio;
__device__ float d_costheta;
__device__ float d_phi;

float h_abmedio;
float h_bcmedio;
float h_costheta;
float h_phi;
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

float norma_merchan;

__inline__ void proyectador_del_perfil(dim3 dimGrid, dim3 dimBlock, curandState *devStates);
__global__ void proyecta_perfil(curandState *state, float r, int eje);

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
  float time,elapsed;
  FILE  *pfout;
  char  filename[200],term[200];
  int   i,j,l;
  int   run;
  float r,s;
  float volumen;
  float h_radio;
  curandState *devStates;

  elapsed = 0.0f;
  chrono(START,&time);

  /*Setea el device a utilizar*/
  cudaSetDevice(DEVICE);

  if(argc < 5){
    printf("Usage: %s abmedio bcmedio costheta phi\n",argv[0]);
    return(EXIT_FAILURE);
  }

  h_abmedio = atof(argv[1]);
  h_bcmedio = atof(argv[2]);
  h_costheta = atof(argv[3]);
  h_phi = atof(argv[4]);
  
  HANDLE_ERROR(cudaMemcpyToSymbol(d_abmedio,&h_abmedio,sizeof(float)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_bcmedio,&h_bcmedio,sizeof(float)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_costheta,&h_costheta,sizeof(float)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_phi,&h_phi,sizeof(float)));

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

  norma_merchan = integra_merchan(dimGrid,dimBlock,devStates);

  proyectador_del_perfil(dimGrid,dimBlock,devStates);

  /*Computa el tiempo total utilizado en el device*/
  chrono(STOP,&time);
  elapsed += time;
  printf("Tiempo: %lf [seg] \n", elapsed);

  /*Fin del programa*/
  return(EXIT_SUCCESS);
}

__global__ void proyecta_perfil(curandState *state, float r, int eje){
  const unsigned int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  const unsigned int it = threadIdx.x;
  int   i,j,k;
  float value, sigma;
  float p[3], tmp;
  float bc,ab;
  float costheta, phi;
  float f;

  //__shared__ float xmin[NDIM], dx[NDIM];
  __shared__ float s_value[THREADS_PER_BLOCK];
  __shared__ float s_sigma[THREADS_PER_BLOCK];

  /*Estado del RNG*/
  curandState seed = state[tid];

  /*Esperan hasta que todos terminen*/
  //__syncthreads();

  value = 0.0f; sigma = 0.0f;
  for(j = 0; j < LAZOS; j++){
    do{
      switch(eje){
        case 2 :
          tmp = 2.0f*PI_CUDA*curand_uniform(&seed);
          p[1] = r*cosf(tmp);
          p[0] = r*sinf(tmp);
          break;

        case 1 :
          tmp = ANGULO*curand_uniform(&seed);
          tmp = tmp*GRAD2RAD;
          p[1] = r*cosf(tmp);
          p[0] = r*sinf(tmp);
          break;

        case 0 :
          tmp = ANGULO*curand_uniform(&seed);
          tmp = tmp*GRAD2RAD;
          p[0] = r*cosf(tmp);
          p[1] = r*sinf(tmp);
          break;
      }

      p[2] = curand_normal(&seed); /*Linea de la visual*/

      tmp   = SQRT_TWOPI_CUDA*RANGO/exp(-p[2]*p[2]*0.5f);
      p[2] *= RANGO;

      //bc = 0.9f*curand_uniform(&seed) + 0.1;
      //ab = 0.9f*curand_uniform(&seed) + 0.1;

      //costheta = 2.0f*curand_uniform(&seed) - 1.0f;
      //phi = 2.0f*PI_CUDA*curand_uniform(&seed);

      bc = d_bcmedio;
      ab = d_abmedio;

      costheta = d_costheta;
      phi = d_phi;

      f = cambio_de_coordenadas(bc,ab,costheta,phi,p);

      tmp *= (u(p,bc,ab,14.0f,1.0f)/exp10f(14.0f)/sqrtf(f));
      //tmp *= forma(bc,ab);

      if(isfinite(tmp))break;
    }while(1);

    value += tmp;
    sigma += tmp*tmp;
  }

  /*Guarda el estado del RNG*/
  state[tid] = seed;

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
}

__inline__ void proyectador_del_perfil(dim3 dimGrid, dim3 dimBlock, curandState *devStates){
  int i,j,k;
  FILE *pfout;
  char filename[80];
  float r,s;
  float dpaso,h_radio;
  dpaso = (PASOMAX - PASOMIN)/(float)NPASOS;

  for(j = 0; j < NDIR; j++){
    sprintf(filename,"perfil_%.2f_%.2f_%.2f_%.2f.%1d",h_abmedio,h_bcmedio,h_costheta,h_phi,j);
    pfout = fopen(filename,"w");

    //En cada direccion hace NPASOS pasos//
    for(i = 0; i < NPASOS; i++){

      //Setea posicion en la direccion dada//
      h_radio = dpaso*(float)(i) + PASOMIN;
      h_radio = exp10(h_radio);

      //Lanza kernel//
      proyecta_perfil<<<dimGrid,dimBlock>>>(devStates,h_radio,j);
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

      //r /= norma_merchan;
      //s /= norma_merchan;

      //Imprime en file de salida//
      fprintf(pfout,"%e %e %e\n",h_radio,r,s);
    }
    //Cierra archivo de salida//
    fclose(pfout);
  }
}

