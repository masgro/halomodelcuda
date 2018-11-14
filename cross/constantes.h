#define OMEGA_M  0.258f    //Omega Matter
#define OMEGA_L  0.742f    //Omega Lambda
#define OMEGA_B  0.0441f   //Omega Barionico
#define OMEGA_N  0.0f      //Omega Neutrinos 
#define FB       1.0f      //0.17093f  //Barion Fraction
#define NUMDEGEN 0.0f      //Densidad numerica de neutrinos degenerados
#define REDSHIFT 0.00f     //Redshift
#define H        0.719f    //Parametro de Hubble
#define TILT     1.0f      //Pendiente del espectro de potencia inicial
#define SIGMA8   0.796f    //rms de densidad en una esfera de 8 mpc

#define RHOCRIT               2.77536627E11
#define TRESCUARTOSSOBREPI    0.23873241f
#define CUATROTERCIOSPI       4.18879020f
#define LOGE10                2.30258509f  //logaritmo de 10 en base e
#define LOG10E                0.43429f      //logaritmo de e en base 10
#define PI_CUDA               3.14159265f
#define TWOPI_CUDA            6.2831853071f
#define SQRT_TWOPI_CUDA       2.50662827f
#define SQRT_TWOPI_CUBO_CUDA 15.74960995f

/*All definitions that should be set at compilation time but they was not*/
#ifndef MAG
#define MAG 19 
#endif
#ifndef BCMEDIO
#define BCMEDIO 0.5 
#endif
#ifndef ABMEDIO
#define ABMEDIO 0.5 
#endif
#ifndef ALIGN_C
#define ALIGN_C 1.4 
#endif
#ifndef ALIGN_B
#define ALIGN_B 1.4 
#endif
#ifndef CENTROS_MASA_MIN
#define CENTROS_MASA_MIN 11.0 
#endif
#ifndef CENTROS_MASA_MAX
#define CENTROS_MASA_MAX 12.0 
#endif
/*************************************************************************/

/*todas estas cosas deberiamos calcularlas adentro*/
/*Scorpio*/
//#define RHOMEDIO 7.1608e+10 //RHOCRIT*omegam
//#define DELTAV   95.4f //Delta_c del paper http://arxiv.org/pdf/1005.0411v1.pdf ec. (9-11)
//#define DELTAC   1.673950f //delta critico de colapso esferico

/*Millenium*/
//#define RHOMEDIO 6.9388e+10 //Millenium
//#define DELTAV   94.215f    //Delta_c del paper http://arxiv.org/pdf/1005.0411v1.pdf ec. (9-11)
//#define DELTAC   1.673660f  //delta critico de colapso esferico

/*TAO*/
#define RHOMEDIO 7.4934889290e+10 //TAO
#define DELTAV   97.010f    
#define DELTAC   1.674369f 

#define FMNORM   0.3232045f 

/*Constantes de la funcion de masa*/
#define ALPHA      0.707f
#define SQRT_ALPHA 0.84083f
#define POTEN      0.3f

/*Constantes de la funcion de bias*/
#define BIAS_ALPHA      0.75f
#define BIAS_POTEN      0.3f

//#define BIAS_ALPHA      0.707f
//#define BIAS_SQRT_ALPHA 0.84083f
//#define BIAS_POTEN      0.3f
//#define BIAS_BETA       0.5f
//#define BIAS_GAMMA      0.6f

#define NuMin -5.704587f 
#define NuMax  5.438179f

#define NORMA_ORIENTACION 0.0795775f //factor para normalizar la func. de probabilidad de
                                     //orientacion del halo vecino (4*pi).
//#define NORMA_ORIENTACION_CROSS 0.10132118364233778 //factor para normalizar la func. de probabilidad de
                                   //orientacion del halo vecino (pi*pi).
//#define NORMA_ORIENTACION_CROSS 0.00806288360829987 //factor para normalizar la func. de probabilidad de
                                   //orientacion del halo vecino (2*pi*pi*2*pi).
#define NORMA_ORIENTACION_CROSS 0.0126651479552922 //(1/(2*pi*2*2*pi)
//#define NORMA_ORIENTACION_PIZZA 0.31831f //1/pi
//#define NORMA_ORIENTACION_PIZZA 0.63662f //1/(pi*0.5) 
//#define NORMA_ORIENTACION_PIZZA 0.4052847f //1/(pi*0.5)**2
//#define NORMA_ORIENTACION_PIZZA 0.0506605918211689 //1/(pi*2*pi)
//#define NORMA_ORIENTACION_PIZZA 0.202642367284676 //1/(pi*0.5*pi)
//#define NORMA_ORIENTACION_PIZZA 0.15915f //1/(2*pi)
//#define NORMA_ORIENTACION_PIZZA  0.07957747154594767f//1/(4*pi)

#define FACTOR_JS 0.50f
#define GRAD2RAD 1.74532925e-2f

#if(MAG==16)
/*HOD Manuel -16*/
#define NSAT_LOGM0   1.591427e+11
#define NSAT_LOGM1   6.555390e+11
#define NSAT_ALPHA   0.9771117
#define NCEN_LOGMMIN 10.64826
#define NCEN_SIGMA   0.2073095
#endif

#if(MAG==17)
/*HOD Manuel -17*/
//#define NSAT_LOGM0   1.847395e+11
//#define NSAT_LOGM1   9.565590e+11
//#define NSAT_ALPHA   0.9722029
//#define NCEN_LOGMMIN 10.91039
//#define NCEN_SIGMA   0.2761359

//SCORPIO//
//#define NSAT_LOGM0   9.440609e+11
//#define NSAT_LOGM1   1.446105e+12
//#define NSAT_ALPHA   0.9546865
//#define NCEN_LOGMMIN 12.08870
//#define NCEN_SIGMA   1.559783

#define NSAT_LOGM0   1.407400e+11
#define NSAT_LOGM1   1.545452e+12
#define NSAT_ALPHA   0.9648788
#define NCEN_LOGMMIN 1.000000e+00
#define NCEN_SIGMA   1.000000e+00

#endif

#if(MAG==18)
/*HOD Manuel -18*/
#define NSAT_LOGM0    8.413953e+10
#define NSAT_LOGM1    1.482277e+12
#define NSAT_ALPHA    0.9673009
#define NCEN_LOGMMIN  11.21975
#define NCEN_SIGMA    0.3569215
#endif

#if(MAG==19)
/*HOD Manuel -19*/
//#define NSAT_LOGM0   8.41395500e+10
//#define NSAT_LOGM1   2.07576800e+12
//#define NSAT_ALPHA   0.918499700
//#define NCEN_LOGMMIN 11.5368800
//#define NCEN_SIGMA   0.399705800

//SCORPIO//
#define NSAT_LOGM0   1.471672e+11
#define NSAT_LOGM1   4.461058e+12
#define NSAT_ALPHA   0.9641127
#define NCEN_LOGMMIN 10.42298
#define NCEN_SIGMA   0.07802054
#endif

#if(MAG==20)
/*HOD Manuel -20*/
#define NSAT_LOGM0    1.883649e+11
#define NSAT_LOGM1    1.477901e+12
#define NSAT_ALPHA    0.6864812
#define NCEN_LOGMMIN  18.87520
#define NCEN_SIGMA    4.850250
#endif

#if(MAG==21)
/*HOD Manuel -21*/
//#define NSAT_LOGM0    1.883650e+11
//#define NSAT_LOGM1    9.489366e+12
//#define NSAT_ALPHA    0.6717666
//#define NCEN_LOGMMIN  20.28347
//#define NCEN_SIGMA    3.936721

//SCORPIO//
#define NSAT_LOGM0   1.478100e+11
#define NSAT_LOGM1   6.855145e+13
#define NSAT_ALPHA   0.9607115
#define NCEN_LOGMMIN 10.42298
#define NCEN_SIGMA   0.07802054
#endif


#ifdef MILLENIUM

#undef NSAT_LOGM0   
#undef NSAT_LOGM1   
#undef NSAT_ALPHA   
#undef NCEN_LOGMMIN 
#undef NCEN_SIGMA   

#if(MAG==17)
#define NSAT_LOGM0   11.227f
#define NSAT_LOGM1   11.83343913f
#define NSAT_ALPHA   0.99538603f
#define NCEN_LOGMMIN 9.84716443f
#define NCEN_SIGMA   0.1f
#endif

#if(MAG==18)
#define NSAT_LOGM0   
#define NSAT_LOGM1   
#define NSAT_ALPHA   
#define NCEN_LOGMMIN 
#define NCEN_SIGMA   
#endif

#if(MAG==19)
#define NSAT_LOGM0   
#define NSAT_LOGM1   
#define NSAT_ALPHA   
#define NCEN_LOGMMIN 
#define NCEN_SIGMA   
#endif

#if(MAG==20)
#define NSAT_LOGM0   
#define NSAT_LOGM1   
#define NSAT_ALPHA   
#define NCEN_LOGMMIN 
#define NCEN_SIGMA   
#endif

#if(MAG==21)
#define NSAT_LOGM0    12.621f
#define NSAT_LOGM1    13.41568715f
#define NSAT_ALPHA    0.91777866f
#define NCEN_LOGMMIN  10.76157022f
#define NCEN_SIGMA    0.1f
#endif
#endif

#ifdef GG
#undef NSAT_LOGM0   
#undef NSAT_LOGM1   
#undef NSAT_ALPHA   
#undef NCEN_LOGMMIN 
#undef NCEN_SIGMA   

#if(MAG==18)
//#define NSAT_LOGM0   8.413953e+10
//#define NSAT_LOGM1   1.482277e+12
//#define NSAT_ALPHA   0.6673009
//#define NCEN_LOGMMIN 11.21975
//#define NCEN_SIGMA   0.3569215
#define NSAT_LOGM0   1.88e+11
#define NSAT_LOGM1   3.0200e+12
#define NSAT_ALPHA   1.01
#define NCEN_LOGMMIN 0.0
#define NCEN_SIGMA   0.0
#endif

#if(MAG==19)
#define NSAT_LOGM0   1.88e+11
#define NSAT_LOGM1   5.8884e+12
#define NSAT_ALPHA   1.06
#define NCEN_LOGMMIN 0.0
#define NCEN_SIGMA   0.0
#endif

#if(MAG==20)
//#define NSAT_LOGM0   1.883649e+11
//#define NSAT_LOGM1   1.477901e+12
//#define NSAT_ALPHA   0.6864812
//#define NCEN_LOGMMIN 18.87520
//#define NCEN_SIGMA   4.850250
#define NSAT_LOGM0   1.88e+11
#define NSAT_LOGM1   1.4125e+13
#define NSAT_ALPHA   1.09
#define NCEN_LOGMMIN 0.0
#define NCEN_SIGMA   0.0
#endif

#if(MAG==21)
//#define NSAT_LOGM0   1.883650e+11
//#define NSAT_LOGM1   9.489366e+12
//#define NSAT_ALPHA   0.6717666
//#define NCEN_LOGMMIN 20.28347
//#define NCEN_SIGMA   3.936721
//
//#define NSAT_LOGM0   1.88e+11
//#define NSAT_LOGM1   6.6069e+12
//#define NSAT_ALPHA   1.13
//#define NCEN_LOGMMIN 0.0
//#define NCEN_SIGMA   0.0
//
#define NSAT_LOGM0   1.478100e+11
#define NSAT_LOGM1   6.855145e+13
#define NSAT_ALPHA   1.13
#define NCEN_LOGMMIN 12.47
#define NCEN_SIGMA   0.07802054
#endif
#endif

#ifdef TAO
#undef NSAT_LOGM0   
#undef NSAT_LOGM1   
#undef NSAT_ALPHA   
#undef NCEN_LOGMMIN 
#undef NCEN_SIGMA   

#if(MAG==16)
#define NSAT_LOGM0   1.9482e+10
#define NSAT_LOGM1   6.8412e+11
#define NSAT_ALPHA   0.97332639f
#define NCEN_LOGMMIN 10.53567766f
#define NCEN_SIGMA   0.29825291f
#endif

#if(MAG==17)
#define NSAT_LOGM0   2.7144e+10
#define NSAT_LOGM1   1.0964e+12
#define NSAT_ALPHA   0.9742907f
#define NCEN_LOGMMIN 10.90751949f
#define NCEN_SIGMA   0.32871516f
#endif

#if(MAG==18)
#define NSAT_LOGM0   1.638406e+11
#define NSAT_LOGM1   1.64967e+12
#define NSAT_ALPHA   0.97023153f
#define NCEN_LOGMMIN 11.07579064f
#define NCEN_SIGMA   0.25656407f
#endif               

#if(MAG==19)
#define NSAT_LOGM0   3.335838e+11
#define NSAT_LOGM1   2.513996e+12 
#define NSAT_ALPHA   0.96177889927 
#define NCEN_LOGMMIN 11.2242686451 
#define NCEN_SIGMA   0.189693397486
#endif               

#if(MAG==20)
#define NSAT_LOGM0   7.1776474e+11
#define NSAT_LOGM1   4.3660788e+12
#define NSAT_ALPHA   0.9551466    
#define NCEN_LOGMMIN 11.3845761   
#define NCEN_SIGMA   0.1741197    
#endif               

#if(MAG==21)
#define NSAT_LOGM0   2.5677576e+12
#define NSAT_LOGM1   1.0731786e+13
#define NSAT_ALPHA   0.9547510    
#define NCEN_LOGMMIN 11.5617118   
#define NCEN_SIGMA   0.1557338    
#endif               
#endif


#ifdef TAO_II
#undef NSAT_LOGM0   
#undef NSAT_LOGM1   
#undef NSAT_ALPHA   
#undef NCEN_LOGMMIN 
#undef NCEN_SIGMA   

#if(MAG==16)
#define NSAT_LOGM0   2.37138e+10 
#define NSAT_LOGM1   9.96539e+11 
#define NSAT_ALPHA   0.97570 
#define NCEN_LOGMMIN  10.82384 
#define NCEN_SIGMA   0.32570
#endif

#if(MAG==17)
#define NSAT_LOGM0   1.55315e+11 
#define NSAT_LOGM1   1.49959e+12 
#define NSAT_ALPHA   0.97063 
#define NCEN_LOGMMIN  11.03274 
#define NCEN_SIGMA   0.26694
#endif

#if(MAG==18)
#define NSAT_LOGM0   2.86721e+11 
#define NSAT_LOGM1   2.26727e+12 
#define NSAT_ALPHA   0.96360 
#define NCEN_LOGMMIN  11.18408 
#define NCEN_SIGMA   0.19388
#endif

#if(MAG==19)
#define NSAT_LOGM0   5.83486e+11 
#define NSAT_LOGM1   3.77472e+12 
#define NSAT_ALPHA   0.95479 
#define NCEN_LOGMMIN  11.25349 
#define NCEN_SIGMA   0.13262
#endif

#if(MAG==20)
#define NSAT_LOGM0   1.69806e+12 
#define NSAT_LOGM1   8.31366e+12 
#define NSAT_ALPHA   0.95235 
#define NCEN_LOGMMIN  11.69348 
#define NCEN_SIGMA   0.24693
#endif

#if(MAG==21)
#define NSAT_LOGM0   4.47362e+12 
#define NSAT_LOGM1   3.16228e+13 
#define NSAT_ALPHA   0.93860 
#define NCEN_LOGMMIN  11.74757 
#define NCEN_SIGMA   0.12878
#endif

#endif

#ifdef SLOAN
#undef NSAT_LOGM0   
#undef NSAT_LOGM1   
#undef NSAT_ALPHA   
#undef NCEN_LOGMMIN 
#undef NCEN_SIGMA   

#if(MAG==16)
#endif

#if(MAG==17)
#define NSAT_LOGM0   1.13559e+13
#define NSAT_LOGM1   5.56929e+12
#define NSAT_ALPHA   1.10000
#define NCEN_LOGMMIN 12.39508
#define NCEN_SIGMA   0.52903
#define NCEN_PLATEAU 2.89390
#endif

#if(MAG==18)
#endif

#if(MAG==19)
#define NSAT_LOGM0   6.05922e+12
#define NSAT_LOGM1   8.79901e+12
#define NSAT_ALPHA   1.09975
#define NCEN_LOGMMIN 12.44040
#define NCEN_SIGMA   0.54055
#define NCEN_PLATEAU 2.30997
#endif

#if(MAG==20)
#endif

#if(MAG==21)
#endif
#endif

#define RMAX 1.0f
#define RMIN -4.0f
