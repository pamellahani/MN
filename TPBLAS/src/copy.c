#include "mnblas.h"
#include "complexe.h"
#include <omp.h> 

void mncblas_scopy(const int N, const float *X, const int incX, float *Y, const int incY)
{
  int i = 0;

  #pragma omp parallel for private (i)
  for (i = 0; i < N ; i += incX)
  {
    register unsigned int j = 0; 
    Y[j] = X[i];
    j+=incY;
  }
}

void mncblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY)
{
  int i = 0;

   #pragma omp parallel for private (i)
  for (i = 0; i < N ; i += incX)
  {
    register unsigned int j = 0; 
    Y[j] = X[i];
    j+=incY; 
  }
}

void mncblas_ccopy(const int N, const void *X, const int incX, void *Y, const int incY)
{
  int i = 0;

   #pragma omp parallel for private (i)
  for (i = 0 ; i < N ; i += incX)
  {
    register unsigned int j = 0; 
    ((complexe_float_t*)Y)[j].real = ((complexe_float_t*)X)[i].real;
    ((complexe_float_t*)Y)[j].imaginary = ((complexe_float_t*)X)[i].imaginary;
    j+= incY; 
  }
}

void mncblas_zcopy(const int N, const void *X, const int incX, void *Y, const int incY)
{
  int i = 0;

   #pragma omp parallel for private (i)
  for (i = 0 ; i < N ; i += incX)
  {
    register unsigned int j = 0; 
    ((complexe_double_t*)Y)[j].real = ((complexe_double_t*)X)[i].real;
    ((complexe_double_t*)Y)[j].imaginary = ((complexe_double_t*)X)[i].imaginary;
    j+= incY; 

  }
}
