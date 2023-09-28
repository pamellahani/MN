#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>

float fabsf(float x)
{
    if (x < 0.0f)
    {
        return -x;
    }
    else
    {
        return x;
    }
}

double fabs(double x)
{
    if (x < 0.0)
    {
        return -x;
    }
    else
    {
        return x;
    }
}

float mnblas_sasum(const int n, const float *x, const int incx)
{
    
    float sum = 0.0;
    int i;
    for (i = 0; i < n; i+=incx)
    {
        sum += fabsf(x[i]);
    }
    return sum;
}

float mnblas_scasum(const int n, const void *x, const int incx)
{
    float sum = 0.0;
    int i;
    const complexe_float_t *cx = (const complexe_float_t *)x;
    for (i = 0; i < n; i+=incx)
    {
        sum += fabsf(cx[i].real) + fabsf(cx[i].imaginary);
    }
    return sum;
}

double mnblas_dasum(const int n, const double *x, const int incx){
    double sum = 0.0;
    int i;
    for (i = 0; i < n; i+=incx)
    {
        sum += fabs(x[i]);
    }
    return sum;
}

double mnblas_dzasum(const int n, const void *x, const int incx){
    double sum = 0.0;
    int i;
    const complexe_double_t *cx = (const complexe_double_t *)x;
    for (i = 0; i < n; i+=incx)
    {
        sum += fabs(cx[i].real) + fabs(cx[i].imaginary);
    }
    return sum;
}