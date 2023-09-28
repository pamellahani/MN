#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>

float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  float dot = 0.0 ;

  for (i = 0 ; i < N ; i += incX)
  {
    dot += X[i] * Y[j] ;
    j += incY ;
  }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX, 
                 const double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  double dot = 0.0 ;

  for (i = 0 ; i < N ; i += incX)
  {
    dot += X[i] * Y[j] ;
    j += incY ;
  }

  return dot ;
}

// unconjugated
void mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  complexe_float_t *dot=(complexe_float_t*)dotu;
  dot->real=0;
  dot->imaginary=0;

  for (i = 0; i < N; i += incX)
  {
    *dot=add_complexe_float(*dot, mult_complexe_float(((complexe_float_t*)X)[i], ((complexe_float_t*)Y)[j]));
    j += incY;
  }
}

// conjugated
void mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  complexe_float_t *dot=(complexe_float_t*)dotc;
  dot->real=0;
  dot->imaginary=0;
  complexe_float_t tmp;

  for (i = 0; i < N; i += incX)
  {
    tmp.real = ((complexe_float_t*)X)[i].real;
    tmp.imaginary = -(((complexe_float_t*)X)[i].imaginary);

    *dot = add_complexe_float(*dot, mult_complexe_float(tmp, ((complexe_float_t*)Y)[j]));
    j += incY;
  }
}

void mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  complexe_double_t *dot=(complexe_double_t*)dotu;
  dot->real=0.0;
  dot->imaginary=0.0;

  for (i = 0; i < N; i += incX)
  {
    *dot=add_complexe_double(*dot, mult_complexe_double(((complexe_double_t*)X)[i], ((complexe_double_t*)Y)[j]));
    j += incY;
  }
}

void mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  complexe_double_t *dot=(complexe_double_t*)dotc;
  dot->real=0.0;
  dot->imaginary=0.0;
  complexe_double_t tmp;

  for (i = 0; i < N; i += incX)
  {
    tmp.real = ((complexe_double_t*)X)[i].real;
    tmp.imaginary = -(((complexe_double_t*)X)[i].imaginary);

    *dot = add_complexe_double(*dot, mult_complexe_double(tmp, ((complexe_double_t*)Y)[j]));
    j += incY;
  }
}