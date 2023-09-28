#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>

void mnblas_saxpy(const int n, const float a, const float *x, const int incx, float *y, const int incy)
{
    register unsigned int i = 0;
    register unsigned int j = 0;

    for (i = 0; i < n; i += incx)
    {
        y[j] += a * x[i] + y[j];
        j += incy;
    }
}

void mnblas_daxpy(const int n, const double a, const double *x, const int incx, double *y, const int incy)
{
    register unsigned int i = 0;
    register unsigned int j = 0;

    for (i = 0; i < n; i += incx)
    {
        y[j] += a * x[i] + y[j];
        j += incy;
    }
}

void mnblas_caxpy(const int n, const void *a, const void *x, const int incx, void *y, const int incy)
{
    register unsigned int i = 0;
    register unsigned int j = 0;
    const float *a_ptr = (const float *)a;
    const complexe_float_t *x_ptr = (const complexe_float_t *)x;
    complexe_float_t *y_ptr = (complexe_float_t *)y;

    for (i = 0; i < n; i += incx)
    {
        y_ptr[i].real += a_ptr[0] * x_ptr[i].real - a_ptr[1] * x_ptr[i].imaginary;
        y_ptr[i].imaginary += a_ptr[0] * x_ptr[i].imaginary + a_ptr[1] * x_ptr[i].real;
    }
}

void mnblas_zaxpy(const int n, const void *a, const void *x, const int incx, void *y, const int incy){
    register unsigned int i = 0;
    register unsigned int j = 0;
    const float *a_ptr = (const float *)a;
    const complexe_double_t *x_ptr = (const complexe_double_t *)x;
    complexe_double_t *y_ptr = (complexe_double_t *)y;

    for (i = 0; i < n; i += incx)
    {
        y_ptr[i].real += a_ptr[0] * x_ptr[i].real - a_ptr[1] * x_ptr[i].imaginary;
        y_ptr[i].imaginary += a_ptr[0] * x_ptr[i].imaginary + a_ptr[1] * x_ptr[i].real;
    }
}