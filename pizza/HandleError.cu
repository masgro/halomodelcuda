#include <stdio.h>
#include "colores.h"
#include "HandleError.h"


void HandleError(cudaError_t err, const char *file, int line){
  if(err != cudaSuccess){
    printf("%s in %s at line %d\n",cudaGetErrorString(err),file,line);
    exit(EXIT_FAILURE);
  }
}

void CheckKernelSuccess(const char *file, int line){
  cudaError_t err;
  err = cudaGetLastError();
	if (err != cudaSuccess){
	  sprintf(message,"Error: %s\n", cudaGetErrorString(err));
		RED(message);
	}
}

