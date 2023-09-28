#include "mnblas.h"
#include "complexe.h"

void mncblas_sswap(const int N, float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX, 
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_cswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY){

    save =((complexe_float_t*)Y)[j].real; 
    ((complexe_float_t*)Y)[j].real = ((complexe_float_t*)X)[i].real; 
    ((complexe_float_t*)X) [i].real = save ;

    save =((complexe_float_t*)Y)[j].imaginary; 
    ((complexe_float_t*)Y)[j].imaginary = ((complexe_float_t*)X)[i].imaginary; 
    ((complexe_float_t*)X) [i].imaginary = save ;

  }
  return ;
}

void mncblas_zswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY){

    save =((complexe_double_t*)Y)[j].real; 
    ((complexe_double_t*)Y)[j].real = ((complexe_double_t*)X)[i].real; 
    ((complexe_float_t*)X) [i].real = save ;

    save =((complexe_double_t*)Y)[j].imaginary; 
    ((complexe_double_t*)Y)[j].imaginary = ((complexe_double_t*)X)[i].imaginary; 
    ((complexe_double_t*)X) [i].imaginary = save ;

  }
  return ;
}

