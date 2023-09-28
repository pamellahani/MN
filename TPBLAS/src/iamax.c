#include "complexe.h"
#include <stdio.h>
#include "mnblas.h"

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

// Finds the index of the element with the largest absolute value.

CBLAS_INDEX mnblas_ismax(const int n, const float *x, const int incx)
{
    CBLAS_INDEX index = 0; 
    float max = 0; 
    for (int i = 0; i < n; i+=incx){
        if (fabsf(x[i]) > max) {
            index = i;
            max = fabsf(x[i]); 
        }
    }
    return index;
}

CBLAS_INDEX mnblas_idamax(const int n, const double *x, const int incx)
{
    CBLAS_INDEX index = 0 ; 
    double max = 0; 
    for (int i = 0; i < n; i+=incx){
        if (fabs(x[i]) > max) {
            index = i;
            max = fabs(x[i]); 
        }
    }
    return index;

}

CBLAS_INDEX mnblas_icamax(const int n, const void *x, const int incx)
{
    CBLAS_INDEX index = 0; 
    float max = 0; 
    for (int i = 0; i < n; i+=incx){
        float abs_value = fabsf(((complexe_float_t *)x)[i].real)+fabsf(((complexe_float_t *)x)[i].imaginary);
        if (abs_value > max) {
            index = i;
            max = abs_value;  
        }
    }
    return index;
}

CBLAS_INDEX mnblas_izamax(const int n, const void *x, const int incx)
{
    CBLAS_INDEX index = 0; 
    double max = 0; 
    for (int i = 0; i < n; i+=incx){
        double abs_value = fabs(((complexe_double_t *)x)[i].real)+fabs(((complexe_double_t *)x)[i].imaginary);
        if (abs_value > max) {
            index = i;
            max = abs_value;  
        }
    }
    return index;
}