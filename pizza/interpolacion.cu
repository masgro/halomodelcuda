#define N_inter 100
__device__ float dev_r[N_inter];
__device__ float dev_f[N_inter];
__device__ int dev_n;

float host_r[N_inter];
float host_f[N_inter];
int host_n;

__device__ float interpolador(float *p){
  int i;
  float r,y;
  r = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
  r = sqrtf(r);
  r = log10f(r);
  i = 0;
  do{
    if(dev_r[i] > r)break;
    i++;
  }while(i < N_inter);
  if(i >= N_inter)return(0.0f);

  y = dev_f[i-1] + (dev_f[i] - dev_f[i-1])/(dev_r[i-1] - dev_r[i])*(dev_r[i-1] - r);
  return(powf(10.0,y));
}
