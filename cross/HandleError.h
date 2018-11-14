#ifndef _HANDLE_ERROR_H_
#define _HANDLE_ERROR_H_

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define CHECK_KERNEL_SUCCESS() (CheckKernelSuccess(__FILE__, __LINE__ ))

void HandleError(cudaError_t, const char *, int);
void CheckKernelSuccess(const char *, int);
#endif
