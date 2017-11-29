#ifndef _FUNCIONES_H_
#define _FUNCIONES_H_
__host__ __device__ void cambio_de_coordenadas(float,float,float,float, float*);
__device__ float T1h(float*,float*);
__device__ float T2h(float*,float*);
__global__ void proyecta_xilin(curandState*,float);
#endif
